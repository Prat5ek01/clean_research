"""Audit CLEAN metrics, false positives, sequence overlap and label drift (E1)."""
import argparse
import time
import warnings
from collections import defaultdict

warnings.filterwarnings("ignore")

from cleanx import config

config.enter_app_dir()

import numpy as np
import pandas as pd

from cleanx import labels as L
from cleanx.index import TrainIndex
from cleanx.metrics import bootstrap, ci, evaluate, group_ids, paired_delta
from cleanx.pipeline import EvalSet, fmt, load_fmap, write_rows
from cleanx.scoring import centroid_distances, maxsep_predict


def label_drift(ev, table, fmap):
    df = pd.read_csv(table, sep="\t", dtype=str, keep_default_na=False, compression="infer")
    df.columns = ["accession", "ec", "date_created", "xref_pfam", "length"][:len(df.columns)]
    now = {a: set(L.complete_ecs(L.split_ecs(e))) for a, e in zip(df.accession, df.ec)}
    if sum(i in now for i in ev.ids) < 0.5 * len(ev.ids):
        return []
    rows = []
    for i in ev.ids:
        bench = set(ev.truth[i])
        cur = now.get(i)
        if cur is None:
            cat = "no_longer_has_EC_or_absent"
        elif cur == bench:
            cat = "same"
        else:
            moved = {s for e in bench for s in (fmap.get(e, [e]) if fmap else [e])}
            if cur == moved:
                cat = "IUBMB_transfer_only"
            elif bench < cur:
                cat = "EC_added"
            elif cur < bench:
                cat = "EC_removed"
            else:
                cat = "EC_changed"
        rows.append({"dataset": ev.name, "Entry": i, "benchmark_EC": ";".join(sorted(bench)),
                     "current_complete_EC": ";".join(sorted(cur)) if cur is not None else "", "category": cat})
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", nargs="+", default=["new392", "price149"])
    ap.add_argument("--backbone", default="split100", choices=list(config.BACKBONES))
    ap.add_argument("--bootstrap", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=config.SEED)
    ap.add_argument("--smoke", action="store_true", help="New-392 only, 200 bootstrap resamples")
    args = ap.parse_args()
    if args.smoke:
        args.datasets, args.bootstrap = ["new392"], 200
    config.ensure_dirs()
    out = config.CX_RESULTS / ("E1_smoke" if args.smoke else "E1")
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    idx = TrainIndex(args.backbone)
    fmap = load_fmap()
    train_seq = idx.train_sequences()
    seq2train = defaultdict(list)
    for tid, s in train_seq.items():
        seq2train[s].append(tid)

    rows, ci_rows, blind_rows, dup_rows, drift_rows, over_rows, summary = [], [], [], [], [], [], []
    for name in args.datasets:
        ev = EvalSet(name, fmap)
        Q = idx.embed(ev.ids)
        D = centroid_distances(Q, idx)
        preds = maxsep_predict(D, idx.ecs, ev.ids)
        order = np.argsort(D, axis=1, kind="stable")
        for k in (1, 2, 3, 5, 10):
            pk = ev.score_labels({q: [idx.ecs[j] for j in order[i, :k]] for i, q in enumerate(ev.ids)})
            r = evaluate(ev.truth, pk, ev.ids, 4)
            over_rows.append({"dataset": name, "ECs_predicted_per_protein": k,
                              "weighted_P(reported)": r["weighted_P"], "weighted_F1(reported)": r["weighted_F1"],
                              "micro_P": r["micro_P"], "micro_F1": r["micro_F1"], "macro_true_F1": r["macro_true_F1"],
                              "macro_union_F1": r["macro_union_F1"], "recall": r["micro_R"],
                              "FP_outside_true_labelset": r["fp_outside_true_labelset"], "FP": r["fp"]})
        L.write_pred_csv(config.CX_PREDS / f"{name}_{args.backbone}_centroid_maxsep{'_smoke' if args.smoke else ''}.csv",
                         preds)
        methods = {"clean_maxsep": ev.score_labels(preds)}
        pv = config.APP_DIR / "results" / "new_pvalue.csv"
        if name == "new392" and args.backbone == "split100" and pv.is_file():
            methods["clean_pvalue(file)"] = ev.score_labels(L.read_pred_csv(pv))
        tiers = ev.tiers(idx)
        for m, P in methods.items():
            for depth in (1, 2, 3, 4):
                r = evaluate(ev.truth, P, ev.ids, depth)
                rows.append({"dataset": name, "backbone": args.backbone, "method": m, "subset_type": "all",
                             "subset": "all", **r})
                for t, sub in group_ids(ev.ids, tiers).items():
                    rows.append({"dataset": name, "backbone": args.backbone, "method": m, "subset_type": "true_tier",
                                 "subset": t, **evaluate(ev.truth, P, sub, depth)})
            r4 = evaluate(ev.truth, P, ev.ids, 4)
            blind_rows.append({"dataset": name, "method": m, "predictions": r4["n_pred_labels"], "TP": r4["tp"],
                               "FP": r4["fp"], "FP_outside_true_labelset": r4["fp_outside_true_labelset"],
                               "share_invisible_to_weighted_P": r4["fp_outside_true_labelset"] / max(r4["fp"], 1),
                               "weighted_P(reported)": r4["weighted_P"], "micro_P": r4["micro_P"],
                               "macro_true_F1": r4["macro_true_F1"], "micro_F1": r4["micro_F1"],
                               "weighted_F1(reported)": r4["weighted_F1"]})
        bs = bootstrap(ev.truth, methods, ev.ids, depth=4, B=args.bootstrap, seed=args.seed)
        for m in methods:
            for key, arr in bs[m].items():
                lo, hi = ci(arr)
                ci_rows.append({"dataset": name, "method": m, "metric": key, "mean": float(arr.mean()),
                                "ci_lo": lo, "ci_hi": hi, "B": args.bootstrap})
        if "clean_pvalue(file)" in methods:
            for key in bs["clean_maxsep"]:
                d, lo, hi, p = paired_delta(bs, "clean_pvalue(file)", "clean_maxsep", key)
                ci_rows.append({"dataset": name, "method": "pvalue_minus_maxsep", "metric": key, "mean": d,
                                "ci_lo": lo, "ci_hi": hi, "B": args.bootstrap, "p_two_sided": p})
        for q in ev.ids:
            hits = seq2train.get(ev.seqs.get(q, ""), [])
            if hits:
                pl = [e for e, _ in preds[q]]
                dup_rows.append({"dataset": name, "Entry": q, "test_EC": ";".join(ev.truth[q]),
                                 "train_ids": ";".join(hits),
                                 "train_EC": " | ".join(";".join(idx.id_ec[h]) for h in hits),
                                 "clean_prediction": ";".join(pl),
                                 "correct": bool(set(ev.score_labels({q: preds[q]})[q]) & set(ev.truth[q]))})
        drift_table = config.CX_DOWNLOADS / "uniprot_sp_ec_all.tsv.gz"
        if not ev.temporal and drift_table.is_file():
            dr_new = label_drift(ev, drift_table, fmap)
            if not dr_new:
                summary.append(f"- label drift {name}: not applicable (ids are not Swiss-Prot accessions)")
            drift_rows += dr_new
        a4 = next(r for r in rows if r["dataset"] == name and r["method"] == "clean_maxsep"
                  and r["subset_type"] == "all" and r["depth"] == 4)
        summary.append(f"- **{name}** (n={a4['n']}): reported (support-weighted) P/R/F1 = "
                       f"{fmt(a4['weighted_P'])}/{fmt(a4['weighted_R'])}/{fmt(a4['weighted_F1'])}; "
                       f"micro P/R/F1 = {fmt(a4['micro_P'])}/{fmt(a4['micro_R'])}/{fmt(a4['micro_F1'])}; "
                       f"macro(true) F1 = {fmt(a4['macro_true_F1'])}; "
                       f"{a4['fp_outside_true_labelset']}/{a4['fp']} FPs invisible to the reported precision; "
                       f"exact training duplicates: {sum(1 for d in dup_rows if d['dataset'] == name)}")

    df = write_rows(rows, out / "E1_metrics_long.csv")
    write_rows(ci_rows, out / "E1_bootstrap_ci.csv")
    write_rows(blind_rows, out / "E1_blind_spot.csv")
    write_rows(dup_rows, out / "E1_exact_duplicates.csv")
    ov = write_rows(over_rows, out / "E1_overprediction.csv")
    for name, g in ov.groupby("dataset"):
        a, b = g.iloc[0], g.iloc[-1]
        summary.append(f"- over-prediction stress test {name}: top-{int(a.ECs_predicted_per_protein)} -> "
                       f"top-{int(b.ECs_predicted_per_protein)} ECs per protein changes reported P "
                       f"{a['weighted_P(reported)']:.3f} -> {b['weighted_P(reported)']:.3f} but micro P "
                       f"{a.micro_P:.3f} -> {b.micro_P:.3f}")
    if drift_rows:
        dr = write_rows(drift_rows, out / "E1_label_drift.csv")
        for name, g in dr.groupby("dataset"):
            summary.append(f"- label drift {name}: " + ", ".join(f"{k}={v}" for k, v in g.category.value_counts().items()))

    td = df[(df.subset_type == "true_tier") & (df.method == "clean_maxsep")]
    piv = td.pivot_table(index=["dataset", "subset"], columns="depth", values="micro_F1")
    piv.to_csv(out / "E1_tier_x_depth_microF1.csv")
    old = config.APP_DIR / "audit_out" / "maxsep_tier_x_depth.csv"
    if "new392" in args.datasets and old.is_file():
        o = pd.read_csv(old).set_index("tier")
        mine = piv.loc["new392"]
        diffs = [abs(mine.loc[t, d] - o.loc[t, f"f1_d{d}"]) for t in o.index for d in (1, 2, 3, 4) if t in mine.index]
        summary.append(f"- consistency with the earlier audit (audit_out/maxsep_tier_x_depth.csv): "
                       f"max |F1 diff| = {max(diffs):.2e} over {len(diffs)} cells")
    with open(out / "E1_summary.md", "w", encoding="utf-8") as fh:
        fh.write(f"# E1 summary ({args.backbone}, bootstrap B={args.bootstrap})\n\n" + "\n".join(summary) + "\n")
    print("\n".join(summary))
    print(f"\nE1 DONE in {time.time() - t0:.0f}s -> {out}")


if __name__ == "__main__":
    main()
