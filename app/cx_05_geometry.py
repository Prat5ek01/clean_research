"""E2: does the single-centroid-per-EC assumption limit CLEAN? Same frozen embeddings, different
inference geometry, all with CLEAN's max-separation rule:
  centroid   CLEAN (one mean per EC)
  exemplar   per-EC 1-NN
  kmean3     local mean of the 3 nearest members (LMKNN)
  multiproto silhouette-gated k-means sub-prototypes
  balanced1/balanced3  exemplar / kmean3 x n_c**lambda (equal-prior kNN density), lambda tuned on
             the calibration split only (--tune-on), never on a test set

Pre-registered endpoints (fixed 2026-09-27, before any sp_temporal_test result existed). Only
metrics that count every false positive are used (E1 shows support-weighted and true-label macro
scores rise when a method over-predicts):
  primary   macro-F1 over the union label space (EC depth 4), balanced3 vs centroid, sp_temporal_test,
            proteins whose ECs are all in CLEAN's vocabulary ('closed'), paired bootstrap, two-sided
  secondary micro-F1 for the same comparison; per-true-EC-tier micro-F1 differences"""
import argparse
import time
import warnings

warnings.filterwarnings("ignore")

from cleanx import config

config.enter_app_dir()

import numpy as np
import pandas as pd

from cleanx import labels as L
from cleanx.index import TrainIndex
from cleanx.metrics import bootstrap, ci, evaluate, group_ids, paired_delta
from cleanx.pipeline import EvalSet, fmt, load_fmap, write_rows
from cleanx.scoring import (balanced_distances, build_multiproto, exemplar_distances, kmean_distances,
                            maxsep_predict, scores_for)

REF = "centroid"
LAMBDAS = [0.0, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2]


def tune_lambda(idx, ev, grid):
    """lambda per k maximising union-label-space macro-F1 (depth 4) on the closed part of the tuning set."""
    Q = idx.embed(ev.ids)
    tiers = ev.tiers(idx)
    closed = [i for i in ev.ids if tiers[i] != "unseen"]
    rows, best = [], {}
    for k in (1, 3):
        base = exemplar_distances(Q, idx) if k == 1 else kmean_distances(Q, idx, k=k)
        for lam in grid:
            preds = maxsep_predict(balanced_distances(Q, idx, k, lam, base=base), idx.ecs, ev.ids)
            r = evaluate(ev.truth, ev.score_labels(preds), closed, 4)
            rows.append({"tuning_set": ev.name, "k": k, "lambda": lam, "n_closed": len(closed),
                         "macro_union_F1": r["macro_union_F1"], "micro_F1": r["micro_F1"],
                         "macro_true_F1": r["macro_true_F1"], "weighted_F1": r["weighted_F1"],
                         "pred_labels_per_protein": r["n_pred_labels"] / max(r["n"], 1)})
        cand = [x for x in rows if x["k"] == k]
        top = max(x["macro_union_F1"] for x in cand)
        best[k] = min(x["lambda"] for x in cand if x["macro_union_F1"] >= top - 1e-12)
    return best, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", nargs="+", default=["new392", "price149", "sp_temporal_test"])
    ap.add_argument("--backbone", default="split100", choices=list(config.BACKBONES))
    ap.add_argument("--methods", nargs="+",
                    default=["centroid", "exemplar", "kmean3", "multiproto", "balanced1", "balanced3"])
    ap.add_argument("--tune-on", default="sp_temporal_cal",
                    help="set used ONLY to choose lambda for balanced1/balanced3")
    ap.add_argument("--lambdas", nargs="+", type=float, default=LAMBDAS)
    ap.add_argument("--bootstrap", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=config.SEED)
    ap.add_argument("--smoke", action="store_true",
                    help="code-path test: New-392 only, lambda tuned on New-392 itself, 200 resamples")
    args = ap.parse_args()
    if args.smoke:
        args.datasets, args.bootstrap, args.tune_on = ["new392"], 200, "new392"
        args.methods = ["centroid", "exemplar", "balanced3"]
    config.ensure_dirs()
    tag = f"E2_{args.backbone}" + ("_smoke" if args.smoke else "")
    out = config.CX_RESULTS / tag
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    idx = TrainIndex(args.backbone)
    fmap = load_fmap()
    lam = {}
    if any(m.startswith("balanced") for m in args.methods):
        try:
            tune_ev = EvalSet(args.tune_on, fmap)
            best, trows = tune_lambda(idx, tune_ev, args.lambdas)
            write_rows(trows, out / "E2_tuning.csv")
            L.save_json({"tuning_set": args.tune_on, "criterion": "macro_union_F1, depth 4, closed subset",
                         "lambda_k1": best[1], "lambda_k3": best[3], "grid": args.lambdas},
                        out / "E2_tuning.json")
            lam = {"balanced1": best[1], "balanced3": best[3]}
            print(f"lambda tuned on {args.tune_on}: balanced1 {best[1]}, balanced3 {best[3]}", flush=True)
            if args.tune_on in args.datasets:
                print(f"  WARNING: {args.tune_on} is also evaluated; its balanced* results are optimistic")
        except FileNotFoundError as e:
            print(f"  skipping balanced1/balanced3: tuning set unavailable ({str(e)[:160]})")
            args.methods = [m for m in args.methods if not m.startswith("balanced")]
    mp = None
    if "multiproto" in args.methods:
        t1 = time.time()
        mp = build_multiproto(idx)
        print(f"multi-prototypes: {int((mp['nk'] > 1).sum())} of {len(idx.ecs)} ECs split "
              f"({time.time() - t1:.0f}s, cached)")
        pd.DataFrame({"ec": idx.ecs, "n": idx.sizes, "k": mp["nk"], "silhouette": mp["sil"]}).to_csv(
            out / "multiproto_ecs.csv", index=False)

    rows, bs_rows, d_rows, per_q, summary = [], [], [], [], []
    for name in args.datasets:
        ev = EvalSet(name, fmap)
        Q = idx.embed(ev.ids)
        tiers = ev.tiers(idx)
        closed = [i for i in ev.ids if tiers[i] != "unseen"]
        idbin = ev.identity() if args.backbone == "split100" else None
        methods, dnn = {}, None
        for m in args.methods:
            S = scores_for(m, Q, idx, mp, lam.get(m))
            if m == "exemplar":
                dnn = S.min(1)
            preds = maxsep_predict(S, idx.ecs, ev.ids)
            L.write_pred_csv(config.CX_PREDS / f"{name}_{args.backbone}_{m}_maxsep{'_smoke' if args.smoke else ''}.csv",
                             preds)
            methods[m] = ev.score_labels(preds)
        hom = config.CX_PREDS / f"{name}_homology_tophit.csv"
        if hom.is_file() and args.backbone == "split100":
            hp = L.read_pred_csv(hom)
            if all(i in hp for i in ev.ids):
                methods["homology_tophit"] = ev.score_labels(hp)
        subsets = {("all", "all"): ev.ids}
        if len(closed) < len(ev.ids):
            subsets[("closed", "seen_ECs_only")] = closed
        for t, sub in group_ids(ev.ids, tiers).items():
            subsets[("true_tier", t)] = sub
        if idbin:
            for b, sub in group_ids(closed, idbin).items():
                subsets[("identity", b)] = sub
        if dnn is not None:
            pos = {i: k for k, i in enumerate(ev.ids)}
            qs = np.quantile(dnn[[pos[i] for i in closed]], [0.25, 0.5, 0.75])
            qlab = {i: f"Q{1 + int(np.searchsorted(qs, dnn[k]))}" for k, i in enumerate(ev.ids)}
            for b, sub in group_ids(closed, qlab).items():
                subsets[("dnn_quartile", b)] = sub
        for (stype, sname), sub in subsets.items():
            for m, P in methods.items():
                for depth in (1, 2, 3, 4):
                    rows.append({"dataset": name, "backbone": args.backbone, "method": m, "subset_type": stype,
                                 "subset": sname, **evaluate(ev.truth, P, sub, depth)})
        primary = closed
        bs = bootstrap(ev.truth, methods, primary, depth=4, B=args.bootstrap, seed=args.seed)
        for m in methods:
            for key, arr in bs[m].items():
                lo, hi = ci(arr)
                bs_rows.append({"dataset": name, "subset": "closed", "method": m, "metric": key,
                                "mean": float(arr.mean()), "ci_lo": lo, "ci_hi": hi, "B": args.bootstrap})
            if m != REF and REF in methods:
                for key in bs[m]:
                    d, lo, hi, p = paired_delta(bs, m, REF, key)
                    d_rows.append({"dataset": name, "subset": "closed", "comparison": f"{m} - {REF}",
                                   "metric": key, "delta": d, "ci_lo": lo, "ci_hi": hi, "p_two_sided": p})
        per_tier = [m for m in ("exemplar", "kmean3", "balanced3") if m in methods]
        if per_tier and REF in methods:
            for t, sub in group_ids(primary, tiers).items():
                if len(sub) < 10:
                    continue
                bt = bootstrap(ev.truth, {k: methods[k] for k in per_tier + [REF]}, sub, depth=4,
                               B=args.bootstrap, seed=args.seed, keys=("micro_F1",))
                for m in per_tier:
                    d, lo, hi, p = paired_delta(bt, m, REF, "micro_F1")
                    d_rows.append({"dataset": name, "subset": f"tier:{t}", "comparison": f"{m} - {REF}",
                                   "metric": "micro_F1", "delta": d, "ci_lo": lo, "ci_hi": hi, "p_two_sided": p})
        for k, i in enumerate(ev.ids):
            r = {"dataset": name, "Entry": i, "tier": tiers[i], "identity_bin": idbin[i] if idbin else "",
                 "d_nn": float(dnn[k]) if dnn is not None else np.nan, "true_EC": ";".join(ev.truth[i])}
            for m, P in methods.items():
                r[f"pred_{m}"] = ";".join(P[i])
                r[f"f1_{m}"] = evaluate(ev.truth, P, [i], 4)["micro_F1"]
            per_q.append(r)
        prim = [r for r in rows if r["dataset"] == name and r["depth"] == 4 and
                r["subset_type"] == ("closed" if len(closed) < len(ev.ids) else "all")]
        line = "; ".join(f"{r['method']} micro-F1 {fmt(r['micro_F1'])} (macro-union {fmt(r['macro_union_F1'])})"
                         for r in prim)
        dl = [d for d in d_rows if d["dataset"] == name and d["subset"] == "closed"
              and d["metric"] in ("micro_F1", "macro_union_F1")]
        dtxt = "; ".join(f"{d['comparison']} {d['metric']}: {d['delta']:+.3f} [{d['ci_lo']:+.3f},{d['ci_hi']:+.3f}] "
                         f"p={d['p_two_sided']:.3g}" for d in dl)
        summary.append(f"- **{name}** (closed n={len(closed)}/{len(ev.ids)}): {line}\n  - paired deltas: {dtxt}")
        print(summary[-1], flush=True)

    write_rows(rows, out / "E2_metrics_long.csv")
    write_rows(bs_rows, out / "E2_bootstrap_ci.csv")
    write_rows(d_rows, out / "E2_deltas.csv")
    write_rows(per_q, out / "E2_per_query.csv")
    with open(out / "E2_summary.md", "w", encoding="utf-8") as fh:
        fh.write(f"# E2 summary ({args.backbone}, B={args.bootstrap})\n\n" + "\n".join(summary) + "\n")
    print(f"\nE2 DONE in {time.time() - t0:.0f}s -> {out}")


if __name__ == "__main__":
    main()
