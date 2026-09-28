"""E5: does the reported EC metric change conclusions across the field? (evaluation card + EC-Bench re-scoring)

PRE-REGISTRATION (written 2026-09-28, before any EC-Bench metric other than the reported one was computed; the
only EC-Bench numbers computed before this text were the reproduction gate below).

Data. EC-Bench (Davoudi et al. 2026, doi:10.1093/bioadv/vbag004, MIT licence) publishes the predictions of 16 EC
predictors, commit 5e98a2c (cleanx/ecbench.py): its test set (468 proteins) and Price-149, for models trained on
its 100%- and 30%-identity clusterings. Labels are parsed exactly as in EC-Bench's own scoring code.
Gate. Our re-implementation must reproduce EC-Bench's printed support-weighted F1 (EC levels 1 and 4, test set,
100% regime, all 16 methods) to 1e-9; otherwise the run stops with FAIL.
Theory. Props. 1-4 in cleanx/evalcard.py: scores computed after projecting predictions onto the test label set,
and support-weighted / true-label-macro scores on any label space, cannot see false positives on EC numbers that
are absent from the evaluated labels; support-weighted recall equals micro recall.

H5 (primary). On the EC-Bench test set (100% regime, EC level 4) at least one pair of methods is ordered one way
by the reported metric (support-weighted F1, EC-Bench's scoring) and significantly the other way by micro-F1.
Primary endpoint: the number of such pairs among the 120, where "significantly" = two-sided paired-bootstrap p
(B = 2000 resamples of proteins, seed 20260927) with Benjamini-Hochberg adjustment across the 120 pairs,
q < 0.05, and the bootstrap median agrees in sign with the observed micro-F1 difference.
Secondary. (S1) the same count with union-label macro-F1 and with per-protein F1; (S2) Kendall tau-b between
the method rankings under the reported metric and each metric that counts every false positive, with a
bootstrap 95% CI; (S3) the share of each method's false positives that the reported metric cannot see;
(S4) Spearman rho between inflation (reported F1 - micro-F1) and mean predicted ECs per protein; (S5) S1-S4 for
the 30% regime, EC levels 1-3, and Price-149 (disclosed: Price-149 was used as the pilot, --pilot); (S6)
complete-EC-only sensitivity (tokens containing '-' dropped); (S7) after the full run, the invisible share of
CLEAN max-separation false positives on SP-temporal TEST, expected below New-392's 0.67 because its label set is
larger (Prop. 4).
Decision rules. >= 1 significant reversal: C1 becomes "the blind spot changes method comparisons in a published
benchmark". No reversal but tau-b < 1: "rankings largely preserved; reported precision inflated by X". tau-b = 1
and no reversal: report the benchmark as robust and keep C1 to CLEAN's own benchmarks and gameability.

AMENDMENT (2026-09-28, after the disclosed Price-149 pilot, before any scoring of the EC-Bench test set beyond
the gate). Pilot: rho(inflation, ECs/protein) = 0.94 (100% regime) and 0.93 (30%); 7 directional but no
significant reversals among 120 pairs (n = 149). Therefore S4 becomes the key secondary endpoint with a
directional hypothesis (rho > 0, one-sided alpha = 0.05, test set, 100% regime, EC level 4). The primary
endpoint is unchanged; the pilot suggests it may be underpowered, and a null result will be reported as such.

Modes: default (confirmatory, downloads the EC-Bench files at the pinned commit), --pilot (EC-Bench Price-149
only; writes E5_pilot), --smoke (synthetic EC-Bench-format data + New-392 CLEAN predictions; writes E5_smoke),
--card TRUTH PREDS (evaluation card for any prediction file; writes E5_card/E5_card_<name>.csv, outside the
confirmatory E5 folder)."""
import argparse
import itertools
import sys
import time
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

from cleanx import config

config.enter_app_dir()

import numpy as np
import pandas as pd
from scipy import stats

from cleanx import ecbench as EB
from cleanx import labels as L
from cleanx.evalcard import card, invariance_check
from cleanx.metrics import bh, bootstrap, ci
from cleanx.nomenclature import forward_counts
from cleanx.pipeline import EvalSet, load_fmap, write_rows

ALL_FP_METRICS = ("micro_F1", "macro_union_F1", "samples_F1")
REPORTED = "weighted_F1"  # union label space = EC-Bench's scoring (and equal to CLEAN's by Prop. 1)
fails = []


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f" -- {detail}" if detail else ""), flush=True)
    if not ok:
        fails.append(name)
    return ok


def complete_only(truth, preds_by_method, ids):
    keep = [i for i in ids if any("-" not in e for e in truth[i])]
    t = {i: [e for e in truth[i] if "-" not in e] for i in keep}
    p = {m: {i: [e for e in P[i] if "-" not in e] for i in keep} for m, P in preds_by_method.items()}
    return keep, t, p


def reproduction_gate(paths, out):
    ids, truth, preds = EB.read_results(paths[("test", "100")])
    rows = []
    for level, printed in EB.PRINTED_TEST100.items():
        for m, v in zip(EB.METHODS, printed):
            mine = card(truth, preds[m], ids, level)["reported_ECBench_F1"]
            rows.append({"method": m, "level": level, "printed": v, "recomputed": mine, "abs_diff": abs(mine - v)})
    df = write_rows(rows, out / "E5_validation.csv")
    return check("EC-Bench reported weighted F1 reproduced (32 values, test set, 100% regime)",
                 df.abs_diff.max() < 1e-9, f"max |diff| {df.abs_diff.max():.1e}")


def rescore(dataset, regime, ids, truth, preds_by_method, depths, B, seed, variant):
    """Cards, pairwise reversals and rank agreement for one prediction table."""
    cards, pairs, agree = [], [], []
    methods = list(preds_by_method)
    for depth in depths:
        cm = {m: card(truth, preds_by_method[m], ids, depth) for m in methods}
        for m in methods:
            cards.append({"source": "ecbench", "dataset": dataset, "regime": regime, "variant": variant,
                          "method": m, **cm[m]})
        bs = bootstrap(truth, preds_by_method, ids, depth, B=B, seed=seed,
                       keys=(REPORTED,) + ALL_FP_METRICS)
        rep_obs = np.array([cm[m]["reported_ECBench_F1"] for m in methods])
        for metric in ALL_FP_METRICS:
            key = "u_" + metric[:-3] + "_F1"
            obs = np.array([cm[m][key] for m in methods])
            rows = []
            for a, b in itertools.combinations(range(len(methods)), 2):
                d = bs[methods[a]][metric] - bs[methods[b]][metric]
                lo, hi = ci(d)
                p = min(1.0, 2 * min(np.mean(d <= 0), np.mean(d >= 0)))
                d_obs = float(obs[a] - obs[b])
                rows.append({"dataset": dataset, "regime": regime, "variant": variant, "depth": depth,
                             "metric": metric, "A": methods[a], "B": methods[b],
                             "reported_diff": float(rep_obs[a] - rep_obs[b]), "diff": d_obs, "ci_lo": lo,
                             "ci_hi": hi, "p_two_sided": float(p),
                             "median_sign_agrees": bool(np.sign(np.median(d)) == np.sign(d_obs))})
            q = bh([r["p_two_sided"] for r in rows])
            for r, qq in zip(rows, q):
                r["q_bh"] = float(qq)
                r["reversed"] = bool(r["reported_diff"] * r["diff"] < 0)
                r["significant_reversal"] = bool(r["reversed"] and qq < 0.05 and r["median_sign_agrees"])
            pairs += rows
            tau = stats.kendalltau(rep_obs, obs)[0]  # tau-b; [0]/[1] indexing works across scipy versions
            R = np.array([bs[m][REPORTED] for m in methods])
            A = np.array([bs[m][metric] for m in methods])
            taus = np.array([stats.kendalltau(R[:, k], A[:, k])[0] for k in range(R.shape[1])])
            taus = taus[np.isfinite(taus)]
            lo, hi = ci(taus) if len(taus) else (float("nan"), float("nan"))
            infl = rep_obs - np.array([cm[m]["u_micro_F1"] for m in methods])
            nlab = np.array([cm[m]["mean_pred_labels"] for m in methods])
            rho = stats.spearmanr(infl, nlab)
            agree.append({"dataset": dataset, "regime": regime, "variant": variant, "depth": depth, "metric": metric,
                          "n_proteins": len(ids), "n_methods": len(methods), "kendall_tau_b": float(tau),
                          "tau_ci_lo": lo, "tau_ci_hi": hi,
                          "n_reversed": int(sum(r["reversed"] for r in rows)),
                          "n_significant_reversals": int(sum(r["significant_reversal"] for r in rows)),
                          "spearman_inflation_vs_labels": float(rho[0]), "spearman_p": float(rho[1]),
                          "spearman_p_one_sided_positive": float(rho[1] / 2 if rho[0] > 0 else 1 - rho[1] / 2)})
    return cards, pairs, agree


def legacy_cards():
    """Cards for our committed CLEAN / geometry / homology predictions (descriptive; no weights needed)."""
    root = config.APP_DIR / "results"
    sets = {"new392": {"clean_maxsep": root / "new_maxsep.csv", "clean_pvalue": root / "new_pvalue.csv"},
            "price149": {"clean_maxsep": root / "inputs" / "price_maxsep.csv"}}
    for name in sets:
        for m in ("centroid", "exemplar", "kmean3", "multiproto"):
            sets[name][f"split100_{m}"] = config.CX_PREDS / f"{name}_split100_{m}_maxsep.csv"
        sets[name]["homology_tophit"] = config.CX_PREDS / f"{name}_homology_tophit.csv"
    cards, inv = [], []
    for name, files in sets.items():
        ids, truth, _ = L.load_named_dataset(name)
        for m, f in files.items():
            if not f.is_file():
                continue
            preds = L.labels_only(L.read_pred_csv(f))
            if not all(i in preds for i in ids):
                continue
            for depth in (1, 2, 3, 4):
                cards.append({"source": "legacy", "dataset": name, "regime": "released_CLEAN", "variant": "as_published",
                              "method": m, **card(truth, preds, ids, depth)})
            inv.append({"dataset": name, "method": m, **invariance_check(truth, preds, ids, 4)})
    return cards, inv


def temporal_cards():
    """S7: cards on SP-temporal TEST once cx_05_geometry.py has written its predictions."""
    cards, inv = [], []
    if not config.dataset_csv("sp_temporal_test").is_file():
        return cards, inv
    fmap = load_fmap()
    ev = EvalSet("sp_temporal_test", fmap)
    for m in ("centroid", "exemplar", "kmean3", "balanced3", "homology_tophit"):
        f = config.CX_PREDS / (f"sp_temporal_test_split100_{m}_maxsep.csv" if m != "homology_tophit"
                               else "sp_temporal_test_homology_tophit.csv")
        if not f.is_file():
            continue
        preds = ev.score_labels(L.read_pred_csv(f))
        if not all(i in preds for i in ev.ids):
            continue
        for depth in (1, 2, 3, 4):
            cards.append({"source": "temporal", "dataset": "sp_temporal_test", "regime": "released_CLEAN",
                          "variant": "harmonized", "method": m, **card(ev.truth, preds, ev.ids, depth)})
        inv.append({"dataset": "sp_temporal_test", "method": m, **invariance_check(ev.truth, preds, ev.ids, 4)})
    return cards, inv


def benchmark_cards():
    """Label-space composition of the temporal sets (C2): unseen ECs split into likely-new EC numbers versus
    existing ECs without a CLEAN training member; pure vs mixed unseen proteins; largest single-EC share."""
    freq = pd.read_csv(config.DATA_DIR / "ec_train_frequency.csv")
    counts = dict(zip(freq.ec, freq.train_count))
    fmap = load_fmap()
    cur = forward_counts(counts, fmap) if fmap else counts
    maxser = {}
    for ec in counts:
        p = ec.split(".")
        if len(p) == 4 and p[3].isdigit():
            k = ".".join(p[:3])
            maxser[k] = max(maxser.get(k, 0), int(p[3]))
    rows = []
    for name in ("sp_temporal_cal", "sp_temporal_test"):
        p = config.dataset_csv(name)
        if not p.is_file():
            continue
        ids, truth, _ = L.read_dataset(p)
        ec_prot = {}
        pure = mixed = new_num = existing = absent_sub = 0
        for i in ids:
            for e in truth[i]:
                ec_prot.setdefault(e, set()).add(i)
            uns = [e for e in truth[i] if cur.get(e, 0) == 0]
            if not uns:
                continue
            if len(uns) == len(truth[i]):
                pure += 1
            else:
                mixed += 1
            cats = set()
            for e in uns:
                parts = e.split(".")
                k = ".".join(parts[:3])
                if k not in maxser:
                    cats.add("absent")
                elif parts[3].isdigit() and int(parts[3]) > maxser[k]:
                    cats.add("new")
                else:
                    cats.add("existing")
            new_num += "new" in cats
            existing += "existing" in cats
            absent_sub += "absent" in cats
        top_ec, top_ids = max(ec_prot.items(), key=lambda kv: len(kv[1]))
        rows.append({"dataset": name, "n_proteins": len(ids), "n_distinct_ECs": len(ec_prot),
                     "unseen_pure": pure, "unseen_mixed": mixed, "unseen_share": (pure + mixed) / len(ids),
                     "unseen_with_EC_numbered_above_vocab": new_num, "unseen_with_existing_EC_no_member": existing,
                     "unseen_with_absent_subsubclass": absent_sub, "largest_EC": top_ec,
                     "largest_EC_proteins": len(top_ids), "largest_EC_share": len(top_ids) / len(ids)})
    return rows


def load_any_preds(path):
    """CLEAN prediction CSV ('id,EC:x/score,...') or plain 'id,EC;EC'."""
    raw = L.labels_only(L.read_pred_csv(path))
    return {q: [t for e in v for t in L.split_ecs(e.replace(",", ";"))] for q, v in raw.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bootstrap", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=config.SEED)
    ap.add_argument("--pilot", action="store_true", help="EC-Bench Price-149 only (disclosed pilot)")
    ap.add_argument("--smoke", action="store_true", help="synthetic EC-Bench-format data; 200 resamples")
    ap.add_argument("--card", nargs=2, metavar=("TRUTH_TSV", "PREDS_CSV"), help="card for any prediction file")
    ap.add_argument("--skip-ecbench", action="store_true",
                    help="cards for our own predictions only (legacy + temporal), no download; writes E5_local")
    ap.add_argument("--force-download", action="store_true")
    args = ap.parse_args()
    config.ensure_dirs()
    t0 = time.time()

    if args.card:
        ids, truth, _ = L.read_dataset(args.card[0])
        preds = load_any_preds(args.card[1])
        missing = [i for i in ids if i not in preds]
        if missing:
            print(f"  note: {len(missing)} proteins without predictions are scored as abstentions")
        rows = [card(truth, preds, ids, d) for d in (1, 2, 3, 4)]
        out = config.CX_RESULTS / "E5_card" / f"E5_card_{Path(args.card[1]).stem}.csv"
        write_rows(rows, out)
        r = rows[-1]
        print(f"EC level 4 (n={r['n']}): reported (support-weighted) P/R/F1 {r['reported_ECBench_P']:.3f}/"
              f"{r['reported_ECBench_R']:.3f}/{r['reported_ECBench_F1']:.3f} | micro {r['u_micro_P']:.3f}/"
              f"{r['u_micro_R']:.3f}/{r['u_micro_F1']:.3f} | union-macro F1 {r['u_macro_union_F1']:.3f} | "
              f"per-protein F1 {r['u_samples_F1']:.3f} | {r['fp_invisible']}/{r['fp']} FPs invisible -> {out}")
        return

    tag = "E5_smoke" if args.smoke else ("E5_pilot" if args.pilot else ("E5_local" if args.skip_ecbench else "E5"))
    out = config.CX_RESULTS / tag
    out.mkdir(parents=True, exist_ok=True)
    B = 200 if args.smoke else args.bootstrap
    cards, pairs, agree, inv = [], [], [], []

    if args.smoke:
        syn = EB.write_synthetic(out / "synthetic_ecbench.csv", n=60, seed=args.seed)
        tables = {("synthetic", "na"): syn}
        cached = {k: EB.DEST / v for k, v in EB.FILES.items()}
        if all(p.is_file() for p in cached.values()):
            reproduction_gate(EB.verify(cached), out)
        else:
            print("  reproduction gate skipped in smoke mode: EC-Bench files not downloaded yet "
                  "(the confirmatory run downloads them)")
    elif args.skip_ecbench:
        tables = {}
    else:
        paths = EB.fetch(force=args.force_download)
        if not reproduction_gate(paths, out):
            print("\nE5 STOPPED: reproduction gate failed")
            sys.exit(1)
        keys = [("price149", "100"), ("price149", "30")] if args.pilot else \
            [("test", "100"), ("test", "30"), ("price149", "100"), ("price149", "30")]
        tables = {k: paths[k] for k in keys}

    for (dataset, regime), path in tables.items():
        ids, truth, preds = EB.read_results(path)
        for variant in ("as_published", "complete_only"):
            if variant == "complete_only":
                ids_v, truth_v, preds_v = complete_only(truth, preds, ids)
            else:
                ids_v, truth_v, preds_v = ids, truth, preds
            depths = (4,) if args.smoke else (1, 2, 3, 4)
            c, p, a = rescore(dataset, regime, ids_v, truth_v, preds_v, depths, B, args.seed, variant)
            cards += c; pairs += p; agree += a
        for m in EB.METHODS:
            inv.append({"dataset": dataset, "method": f"ecbench:{regime}:{m}",
                        **invariance_check(truth, preds[m], ids, 4)})
        print(f"  {dataset}/{regime}: {len(ids)} proteins x {len(EB.METHODS)} methods scored "
              f"({time.time() - t0:.0f}s)", flush=True)

    lc, li = legacy_cards()
    cards += lc; inv += li
    if not args.smoke and not args.pilot:
        tc, ti = temporal_cards()
        cards += tc; inv += ti
    bench = benchmark_cards()

    cdf = write_rows(cards, out / "E5_cards.csv")
    pdf = write_rows(pairs, out / "E5_pairs.csv")
    adf = write_rows(agree, out / "E5_rank_agreement.csv")
    idf = write_rows(inv, out / "E5_invariance.csv")
    write_rows(bench, out / "E5_benchmark_card.csv")
    check("Propositions 1-2 and Lemma 3 hold numerically on every prediction table",
          bool(len(idf)) and bool(idf.ok.all()), f"{int(idf.ok.sum())}/{len(idf)} tables")

    lines = [f"# {tag} summary (B={B}, seed {args.seed})\n"]
    if args.smoke:
        lines.append("Synthetic EC-Bench-format data: code-path test only, no scientific content.\n")
    if not len(adf):  # --skip-ecbench: no EC-Bench tables, so no rank-agreement rows
        adf = pd.DataFrame(columns=["dataset", "regime", "depth", "variant", "metric"])
    prim = adf[(adf.dataset == "test") & (adf.regime == "100") & (adf.depth == 4) &
               (adf.variant == "as_published") & (adf.metric == "micro_F1")]
    if len(prim):
        r = prim.iloc[0]
        lines.append(f"- **PRIMARY (H5)**: EC-Bench test, 100% regime, EC level 4: {int(r.n_significant_reversals)} "
                     f"significantly reversed pairs (reported weighted F1 vs micro-F1, BH q<0.05) of 120; "
                     f"{int(r.n_reversed)} reversed in direction; Kendall tau-b {r.kendall_tau_b:.3f} "
                     f"[{r.tau_ci_lo:.3f}, {r.tau_ci_hi:.3f}]")
        lines.append(f"- **KEY SECONDARY (S4)**: Spearman rho(inflation, ECs per protein) {r.spearman_inflation_vs_labels:.2f}, "
                     f"one-sided p {r.spearman_p_one_sided_positive:.3g}")
    for r in adf[(adf.depth == 4) & (adf.variant == "as_published")].itertuples():
        lines.append(f"- {r.dataset}/{r.regime} EC4 {r.metric}: tau-b {r.kendall_tau_b:.3f} "
                     f"[{r.tau_ci_lo:.3f},{r.tau_ci_hi:.3f}], reversed {r.n_reversed}, significant "
                     f"{r.n_significant_reversals}; Spearman(inflation, ECs/protein) {r.spearman_inflation_vs_labels:.2f} "
                     f"(p={r.spearman_p:.3g})")
    eb4 = cdf[(cdf.source == "ecbench") & (cdf.depth == 4) & (cdf.variant == "as_published")]
    for (ds, rg), g in eb4.groupby(["dataset", "regime"]):
        lines.append(f"- {ds}/{rg} EC4 invisible share of FPs across methods: median "
                     f"{g.invisible_share.median():.2f} (range {g.invisible_share.min():.2f}-{g.invisible_share.max():.2f})")
    for r in cdf[(cdf.source != "ecbench") & (cdf.depth == 4)].itertuples():
        share = "n/a" if r.invisible_share != r.invisible_share else f"{r.invisible_share:.0%}"
        lines.append(f"- {r.source} {r.dataset} {r.method}: reported F1 {r.reported_ECBench_F1:.3f} vs micro "
                     f"{r.u_micro_F1:.3f}; {r.fp_invisible}/{r.fp} FPs invisible ({share})")
    for r in bench:
        lines.append(f"- benchmark card {r['dataset']}: {r['n_proteins']} proteins, {r['n_distinct_ECs']} ECs, unseen "
                     f"{r['unseen_share']:.1%} (pure {r['unseen_pure']}, mixed {r['unseen_mixed']}; EC numbered above "
                     f"CLEAN's vocabulary {r['unseen_with_EC_numbered_above_vocab']}, existing EC without training "
                     f"member {r['unseen_with_existing_EC_no_member']}); largest EC {r['largest_EC']} "
                     f"({r['largest_EC_share']:.1%})")
    with open(out / "E5_summary.md", "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\n{'CX09 SMOKE' if args.smoke else 'E5'} {'FAILED' if fails else 'DONE'} in {time.time() - t0:.0f}s -> {out}")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
