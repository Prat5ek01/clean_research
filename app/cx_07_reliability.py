"""E4: reliability. Is CLEAN's released confidence calibrated across rarity tiers, and can
tier-conditional (Mondrian) risk control give certified abstention per tier, including on
unseen-EC enzymes and non-enzymes?

Top-1 prediction per protein under each inference geometry (centroid = CLEAN; exemplar; balanced3
with the lambda E2 tuned on the calibration split).
- calibration: CLEAN's GMM confidence as released; isotonic maps fitted on the calibration split
  (global, and per tier of the predicted EC); ECE/Brier on the test sets
- risk control: SGR thresholds (Geifman & El-Yaniv 2017) on the raw distance, global vs one
  threshold per predicted-EC tier (delta split by Bonferroni); realized risk/coverage on test
- open set: share of unseen-EC enzymes and non-enzymes abstained at the certified thresholds;
  AUROC of each score for seen-EC enzymes vs non-enzymes
- hierarchical fallback: when EC4 abstains, emit the EC3 prefix if its own certified threshold passes

Changes recorded 2026-09-28, before any SP-temporal result existed:
- unseen-EC enzymes are split into 'pure' (no true EC in CLEAN's vocabulary: abstaining is the right answer) and
  'mixed' (at least one true EC is in the vocabulary, so a correct call exists; e.g. proteins annotated with both
  1.14.99.54 and 1.14.99.56); abstention is reported for both, plus accuracy among accepted mixed calls
- the EC3 fallback threshold is certified on the calibration proteins that the EC4 rule rejects, the population it
  is applied to (it was certified on all calibration proteins before, which does not cover the rejected subset)"""
import argparse
import time
import warnings

warnings.filterwarnings("ignore")

from cleanx import config

config.enter_app_dir()

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from cleanx import labels as L
from cleanx.calibration import (IsoCalibrator, MondrianIso, apply_threshold, brier, gmm_confidence,
                                gmm_confidence_fast, load_gmms, reliability, risk_coverage,
                                selective_summary, sgr_threshold)
from cleanx.index import TrainIndex
from cleanx.nomenclature import forward_counts
from cleanx.pipeline import EvalSet, load_fmap, write_rows
from cleanx.scoring import scores_for

PRED_TIERS = ["common", "rare", "very_rare", "singleton"]


def features(ev, idx, gmms, geoms, lam, ids=None):
    ids = ev.ids if ids is None else ids
    Q = idx.embed(ids)
    true_tier = ev.tiers(idx)
    counts = forward_counts(idx.counts, ev.fmap) if ev.fmap else idx.counts
    kind = {i: ("" if true_tier[i] != "unseen" or not ev.truth[i] else
                ("mixed" if any(counts.get(e, 0) > 0 for e in ev.truth[i]) else "pure")) for i in ids}
    frames = {}
    for geom in geoms:
        S = scores_for(geom, Q, idx, lam=lam.get(geom))
        o2 = np.argpartition(S, 1, axis=1)[:, :2]
        d_o2 = np.take_along_axis(S, o2, 1)
        swap = d_o2[:, 1] < d_o2[:, 0]
        o2[swap] = o2[swap][:, ::-1]
        d1 = S[np.arange(len(ids)), o2[:, 0]]
        d2 = S[np.arange(len(ids)), o2[:, 1]]
        rows = []
        for k, i in enumerate(ids):
            ec = idx.ecs[o2[k, 0]]
            truth = set(ev.truth[i])
            lab = L.harmonize([ec], truth, ev.fmap)[0] if ev.fmap else ec
            rows.append({"Entry": i, "top1": ec, "top1_scored_as": lab, "d1": d1[k], "gap": d2[k] - d1[k],
                         "pred_tier": idx.tiers[o2[k, 0]], "true_tier": true_tier[i] if truth else "non_enzyme",
                         "unseen_kind": kind[i], "correct4": lab in truth,
                         "correct3": L.truncate(lab, 3) in {L.truncate(t, 3) for t in truth}})
        f = pd.DataFrame(rows)
        f["score"] = -f.d1
        if geom == "centroid":
            f["clean_gmm"] = gmm_confidence_fast(f.d1.to_numpy(), gmms)
        frames[geom] = f
    return frames


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cal", default="sp_temporal_cal")
    ap.add_argument("--tests", nargs="+", default=["sp_temporal_test", "price149", "new392"])
    ap.add_argument("--nonenzyme", default="sp_nonenzyme")
    ap.add_argument("--backbone", default="split100", choices=["split100"])
    ap.add_argument("--geometries", nargs="+", default=["centroid", "exemplar", "balanced3"])
    ap.add_argument("--alphas", nargs="+", type=float, default=[0.1, 0.2, 0.3])
    ap.add_argument("--delta", type=float, default=0.1)
    ap.add_argument("--mondrian-min-n", type=int, default=50)
    ap.add_argument("--seed", type=int, default=config.SEED)
    ap.add_argument("--smoke", action="store_true", help="New-392 split in halves as cal/test (code path only)")
    args = ap.parse_args()
    config.ensure_dirs()
    out = config.CX_RESULTS / ("E4_smoke" if args.smoke else "E4")
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    if "centroid" not in args.geometries:
        args.geometries = ["centroid"] + args.geometries
    lam = {}
    if any(g.startswith("balanced") for g in args.geometries):
        tj = config.CX_RESULTS / (f"E2_{args.backbone}" + ("_smoke" if args.smoke else "")) / "E2_tuning.json"
        if tj.is_file():
            t = L.load_json(tj)
            lam = {"balanced1": t["lambda_k1"], "balanced3": t["lambda_k3"]}
            print(f"using lambda from {tj.parent.name} (tuned on {t['tuning_set']}): {lam}")
        else:
            print(f"  skipping balanced geometries: {tj} not found (run cx_05_geometry.py first)")
            args.geometries = [g for g in args.geometries if not g.startswith("balanced")]
    geoms = args.geometries

    idx = TrainIndex(args.backbone)
    fmap = load_fmap()
    gmms = load_gmms()
    probe = np.linspace(0.5, 12, 24)
    if np.max(np.abs(gmm_confidence(probe, gmms) - gmm_confidence_fast(probe, gmms))) > 1e-12:
        raise SystemExit("vectorised GMM confidence differs from CLEAN.evaluate.infer_confidence_gmm")

    feats = {}
    if args.smoke:
        ev = EvalSet("new392", fmap)
        rng = np.random.default_rng(args.seed)
        perm = list(rng.permutation(ev.ids))
        half = len(perm) // 2
        feats["cal"] = features(ev, idx, gmms, geoms, lam, perm[:half])
        feats["new392_half"] = features(ev, idx, gmms, geoms, lam, perm[half:])
        tests, nonenz = ["new392_half"], None
    else:
        feats["cal"] = features(EvalSet(args.cal, fmap), idx, gmms, geoms, lam)
        for name in args.tests:
            feats[name] = features(EvalSet(name, fmap), idx, gmms, geoms, lam)
        tests = list(args.tests)
        nonenz = None
        if args.nonenzyme and config.dataset_csv(args.nonenzyme).is_file():
            nonenz = features(EvalSet(args.nonenzyme, fmap), idx, gmms, geoms, lam)
            feats[args.nonenzyme] = nonenz
    for name, fr in feats.items():
        for geom, f in fr.items():
            f.to_csv(out / f"E4_features_{name}_{geom}.csv", index=False)
    print(f"features ready ({time.time() - t0:.0f}s)", flush=True)

    cal_rows, rc_rows, aurc_rows, os_rows, fb_rows = [], [], [], [], []
    calib = {}
    for geom in geoms:
        c = feats["cal"][geom]
        calib[(geom, "global")] = IsoCalibrator().fit(c.score, c.correct4)
        calib[(geom, "mondrian")] = MondrianIso(args.mondrian_min_n).fit(c.score, c.correct4, c.pred_tier)

    def conf_of(f, geom, kind):
        if kind == "clean_gmm(raw)":
            return f.clean_gmm.to_numpy()
        cal = calib[(geom, kind)]
        return cal.predict(f.score) if kind == "global" else cal.predict(f.score, f.pred_tier)

    variants = [("centroid", "clean_gmm(raw)")] + [(g, k) for g in geoms for k in ("global", "mondrian")]
    for name in tests:
        for geom, kind in variants:
            f = feats[name][geom]
            conf = conf_of(f, geom, kind)
            for grp_col in ("all", "pred_tier", "true_tier"):
                groups = {"all": np.ones(len(f), bool)} if grp_col == "all" else \
                    {g: (f[grp_col] == g).to_numpy() for g in f[grp_col].unique()}
                for g, m in groups.items():
                    if m.sum() < 10:
                        continue
                    _, e = reliability(conf[m], f.correct4[m])
                    cal_rows.append({"dataset": name, "geometry": geom, "confidence": kind, "group_by": grp_col,
                                     "group": g, "n": int(m.sum()), "accuracy": float(f.correct4[m].mean()),
                                     "mean_conf": float(conf[m].mean()), "ECE": e,
                                     "brier": brier(conf[m], f.correct4[m])})
            for subset, m in (("all", np.ones(len(f), bool)), ("seen_ECs_only", (f.true_tier != "unseen").to_numpy())):
                if m.sum():
                    _, _, aurc = risk_coverage(conf[m], f.correct4[m])
                    aurc_rows.append({"dataset": name, "geometry": geom, "confidence": kind, "subset": subset,
                                      "n": int(m.sum()), "AURC": aurc, "top1_accuracy": float(f.correct4[m].mean())})

    for geom in geoms:
        c = feats["cal"][geom]
        for alpha in args.alphas:
            g_theta = sgr_threshold(c.score, c.correct4, alpha, args.delta)[0]
            th = {"global": {t: g_theta for t in PRED_TIERS}, "mondrian": {}}
            for t in PRED_TIERS:
                m = (c.pred_tier == t).to_numpy()
                th["mondrian"][t] = sgr_threshold(c.score[m], c.correct4[m], alpha, args.delta / len(PRED_TIERS))[0]
            cal_rej = ~apply_threshold(c.score, c.pred_tier.map(th["mondrian"]).to_numpy())
            th3 = (sgr_threshold(c.score[cal_rej], c.correct3[cal_rej], alpha, args.delta)[0]
                   if cal_rej.any() else np.inf)
            for rule, tt in th.items():
                for name in tests + ([args.nonenzyme] if nonenz is not None else []):
                    f = feats[name][geom]
                    theta = f.pred_tier.map(tt).to_numpy()
                    acc = apply_threshold(f.score, theta)
                    base = {"geometry": geom, "alpha": alpha, "delta": args.delta, "rule": rule, "dataset": name}
                    if name == args.nonenzyme:
                        os_rows.append({**base, "group": "non_enzyme", "n": len(f),
                                        "abstained": float(1 - acc.mean())})
                        continue
                    rc_rows.append({**base, "group": "all", **selective_summary(acc, f.correct4)})
                    for t in PRED_TIERS:
                        m = (f.pred_tier == t).to_numpy()
                        if m.any():
                            rc_rows.append({**base, "group": f"pred_tier:{t}",
                                            **selective_summary(acc[m], f.correct4[m]),
                                            "theta": float(tt[t])})
                    for t in f.true_tier.unique():
                        m = (f.true_tier == t).to_numpy()
                        rc_rows.append({**base, "group": f"true_tier:{t}", **selective_summary(acc[m], f.correct4[m])})
                    un = (f.true_tier == "unseen").to_numpy()
                    if un.any():
                        os_rows.append({**base, "group": "unseen_EC_enzymes", "n": int(un.sum()),
                                        "abstained": float(1 - acc[un].mean())})
                    for kd in ("pure", "mixed"):
                        mk = (f.unseen_kind == kd).to_numpy()
                        if mk.any():
                            ok = mk & acc
                            os_rows.append({**base, "group": f"unseen_EC_enzymes_{kd}", "n": int(mk.sum()),
                                            "abstained": float(1 - acc[mk].mean()),
                                            "accuracy_if_accepted": float(f.correct4[ok].mean()) if ok.any() else np.nan})
                    if rule == "mondrian":
                        fb = ~acc & apply_threshold(f.score, th3)
                        fb_rows.append({**base, "ec4_accepted": int(acc.sum()), "ec3_fallback_accepted": int(fb.sum()),
                                        "ec3_fallback_risk": float((~f.correct3[fb]).mean()) if fb.any() else np.nan,
                                        "coverage_any_level": float((acc | fb).mean()), "theta_ec3": float(th3)})

    minc_rows = []
    grid = np.round(np.arange(0.05, 0.951, 0.05), 2)
    for geom in geoms:
        c = feats["cal"][geom]
        groups = [("all (global threshold)", np.ones(len(c), bool), args.delta)] + \
                 [(f"pred_tier:{t}", (c.pred_tier == t).to_numpy(), args.delta / len(PRED_TIERS)) for t in PRED_TIERS]
        for gname, m, dlt in groups:
            best = None
            for a in grid:
                theta, bound, cov = sgr_threshold(c.score[m], c.correct4[m], float(a), dlt)
                if np.isfinite(theta):
                    best = (float(a), bound, cov)
                    break
            minc_rows.append({"geometry": geom, "group": gname, "n_cal": int(m.sum()), "delta": dlt,
                              "min_certifiable_risk": best[0] if best else np.nan,
                              "risk_bound": best[1] if best else np.nan,
                              "cal_coverage_at_that_risk": best[2] if best else 0.0})
    write_rows(minc_rows, out / "E4_min_certifiable_risk.csv")

    if nonenz is not None:
        for name in tests:
            for geom in geoms:
                f = feats[name][geom]
                seen = f[f.true_tier != "unseen"]
                y = np.r_[np.ones(len(seen)), np.zeros(len(feats[args.nonenzyme][geom]))]
                for kind in ("score", "clean_gmm") if geom == "centroid" else ("score",):
                    s = np.r_[seen[kind].to_numpy(), feats[args.nonenzyme][geom][kind].to_numpy()]
                    os_rows.append({"geometry": geom, "dataset": name,
                                    "group": f"AUROC seen-EC enzymes vs non-enzymes ({kind})",
                                    "n": len(y), "auroc": float(roc_auc_score(y, s))})
                if (f.true_tier == "unseen").any():
                    y2 = (f.true_tier != "unseen").to_numpy().astype(int)
                    os_rows.append({"geometry": geom, "dataset": name,
                                    "group": "AUROC seen-EC vs unseen-EC enzymes (score)",
                                    "n": len(f), "auroc": float(roc_auc_score(y2, f.score))})

    cal_df = write_rows(cal_rows, out / "E4_calibration.csv")
    rc_df = write_rows(rc_rows, out / "E4_risk_control.csv")
    write_rows(aurc_rows, out / "E4_aurc.csv")
    write_rows(os_rows, out / "E4_openset.csv")
    write_rows(fb_rows, out / "E4_fallback.csv")

    lines = [f"# E4 summary ({'smoke: New-392 halves' if args.smoke else args.cal + ' -> ' + ', '.join(tests)})\n",
             f"geometries: {geoms}; lambda: {lam}\n"]
    main_test = tests[0]
    g = cal_df[(cal_df.dataset == main_test) & (cal_df.group_by == "pred_tier")]
    for geom, kind in [("centroid", "clean_gmm(raw)")] + [(x, "mondrian") for x in geoms]:
        s = g[(g.confidence == kind) & (g.geometry == geom)]
        lines.append(f"- {main_test} ECE by predicted-EC tier, {geom}/{kind}: " +
                     ", ".join(f"{r.group} {r.ECE:.3f} (acc {r.accuracy:.2f}, conf {r.mean_conf:.2f}, n={r.n})"
                               for r in s.itertuples()))
    a_show = 0.2 if 0.2 in args.alphas else args.alphas[0]
    r = rc_df[(rc_df.dataset == main_test) & (rc_df.alpha == a_show)]
    for rule in ("global", "mondrian"):
        for geom in geoms:
            s = r[(r.rule == rule) & (r.geometry == geom) & (r.group.str.startswith("pred_tier") | (r.group == "all"))]
            if len(s):
                lines.append(f"- risk control alpha={a_show} {rule}/{geom}: " + ", ".join(
                    f"{x.group} risk {x.selective_risk:.3f} cov {x.coverage:.2f}" for x in s.itertuples()))
    for geom in geoms:
        lines.append(f"- smallest certifiable selective risk (calibration split), {geom}: " + ", ".join(
            f"{x['group']} {x['min_certifiable_risk']} (cal coverage {x['cal_coverage_at_that_risk']:.2f}, n={x['n_cal']})"
            for x in minc_rows if x["geometry"] == geom))
    with open(out / "E4_summary.md", "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\nE4 DONE in {time.time() - t0:.0f}s -> {out}")


if __name__ == "__main__":
    main()
