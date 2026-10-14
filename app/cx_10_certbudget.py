"""Estimate model-based calibration budgets by predicted rarity tier (E6).

See docs/analysis_protocols.md#cx-10-certbudget."""
import argparse
import math
import sys
import time
import warnings

warnings.filterwarnings("ignore")

from cleanx import config

config.enter_app_dir()

import numpy as np
import pandas as pd

from cleanx import certify as C
from cleanx.labels import load_json
from cleanx.pipeline import write_rows

TIERS = ["common", "rare", "very_rare", "singleton"]
GEOMS = ["centroid", "exemplar", "balanced3"]
fails = []


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f" -- {detail}" if detail else ""), flush=True)
    if not ok:
        fails.append(name)
    return ok


def fmt_n(x):
    return "inf" if not np.isfinite(x) else f"{int(x)}"


def synthetic_features(n=800, seed=0):
    rng = np.random.default_rng(seed)
    tier = rng.choice(TIERS, size=n, p=[0.6, 0.2, 0.1, 0.1])
    base = {"common": 0.75, "rare": 0.55, "very_rare": 0.4, "singleton": 0.35}
    s = rng.uniform(size=n)
    p = np.array([min(0.98, base[t] + 0.5 * (x - 0.5)) for t, x in zip(tier, s)])
    return pd.DataFrame({"score": s, "correct4": rng.uniform(size=n) < p, "pred_tier": tier})


def accrual(cal_tier_share):
    """New Swiss-Prot enzymes per year (benchmark filters) and per-tier yearly counts."""
    meta = pd.read_csv(config.CX_SETS / "sp_temporal_meta.csv", dtype=str, keep_default_na=False)
    enz = meta[meta.tier != "non_enzyme"]
    man = config.CX_SETS / "sp_temporal_manifest.json"
    end = None
    if man.is_file():
        end = (load_json(man).get("uniprot") or {}).get("uniprot_release_date")
    end = pd.to_datetime(end) if end else pd.to_datetime(enz.date_created).max()
    years = (end - pd.to_datetime(config.SNAPSHOT_DATE)).days / 365.25
    total = len(enz) / years
    rows = [{"group": "all_enzymes", "per_year": total, "years_span": years, "basis": "benchmark filters"}]
    for t in TIERS:
        rows.append({"group": f"pred_tier:{t}", "per_year": total * cal_tier_share.get(t, 0.0), "years_span": years,
                     "basis": "all enzymes x tier share of calibration predictions"})
        rows.append({"group": f"true_tier:{t}", "per_year": (enz.tier == t).sum() / years, "years_span": years,
                     "basis": "rarest true EC tier"})
    return pd.DataFrame(rows)


def tier_budgets(f, geom, alphas, delta, B, seed, source):
    rows = []
    rng = np.random.default_rng(seed)
    t0 = time.time()
    for t in TIERS:
        g = f[f.pred_tier == t]
        n_av = len(g)
        print(f"  {geom}/{t}: {n_av} calibration proteins, {B} bootstrap resamples x {len(alphas)} alphas "
              f"({time.time() - t0:.0f}s elapsed)", flush=True)
        for a in alphas:
            if n_av == 0:
                rows.append({"source": source, "geometry": geom, "tier": t, "alpha": a, "n_available": 0})
                continue
            est = C.best_budget(g.score.to_numpy(), g.correct4.to_numpy(), a, delta, groups=len(TIERS))
            boot = []
            for _ in range(B):
                idx = rng.integers(0, n_av, size=n_av)
                boot.append(C.best_budget(g.score.to_numpy()[idx], g.correct4.to_numpy()[idx], a, delta,
                                          groups=len(TIERS))["n_needed"])
            boot = np.array(boot, float)
            lo, hi = C.percentile_ci(boot)  # non-certifiable resamples count as +inf (pre-registered)
            rows.append({"source": source, "geometry": geom, "tier": t, "alpha": a, "n_available": n_av,
                         "top1_accuracy": float(g.correct4.mean()), **est,
                         "deficit": est["n_needed"] - n_av, "deficit_ci_lo": lo - n_av, "deficit_ci_hi": hi - n_av,
                         "share_resamples_certifiable": float(np.isfinite(boot).mean()) if B else np.nan})
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--alphas", nargs="+", type=float, default=[0.1, 0.2, 0.3])
    ap.add_argument("--delta", type=float, default=0.1)
    ap.add_argument("--power", type=float, default=0.8)
    ap.add_argument("--bootstrap", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=config.SEED)
    ap.add_argument("--pilot", action="store_true", help="E4 smoke features (New-392 half) as calibration data")
    ap.add_argument("--smoke", action="store_true", help="synthetic features + known-answer checks")
    args = ap.parse_args()
    config.ensure_dirs()
    t0 = time.time()
    tag = "E6_smoke" if args.smoke else ("E6_pilot" if args.pilot else "E6")
    out = config.CX_RESULTS / tag
    out.mkdir(parents=True, exist_ok=True)

    # known answers (every mode)
    n0 = C.n_star(0.0, 0.2, 0.025)[0]
    check("n* at zero risk equals floor(ln delta / ln(1 - alpha)) + 1", n0 == 17, f"n*(r=0, alpha=0.2, delta=0.025) = {n0}")
    check("no budget certifies a risk at or above alpha", not math.isfinite(C.n_star(0.2, 0.2, 0.025)[1]))
    grid = []
    dp = C.per_test_delta(args.delta, len(TIERS), 256)
    for a in args.alphas:
        for r in np.round(np.arange(0.0, a, 0.01), 2):
            first, stable = C.n_star(float(r), a, dp, args.power)
            grid.append({"alpha": a, "true_selective_risk": float(r), "per_test_delta": dp, "n_first": first,
                         "n_stable": stable})
    gdf = write_rows(grid, out / "E6_nstar_grid.csv")
    mono = all((g.n_stable.diff().dropna() >= 0).all() for _, g in gdf.replace(np.inf, 1e12).groupby("alpha"))
    check("budget n* is non-decreasing in the true selective risk", mono)

    if args.smoke:
        feats = {"centroid": synthetic_features(seed=args.seed)}
        B = 50
    else:
        base = config.CX_RESULTS / ("E4_smoke" if args.pilot else "E4")
        name = "cal"
        feats = {}
        for g in GEOMS:
            p = base / f"E4_features_{name}_{g}.csv"
            if p.is_file():
                feats[g] = pd.read_csv(p)
        if "centroid" not in feats:
            print(f"E6 needs {base / 'E4_features_cal_centroid.csv'} (run cx_07_reliability.py"
                  f"{' --smoke' if args.pilot else ''} first)")
            sys.exit(1)
        B = args.bootstrap
    source = "synthetic" if args.smoke else ("pilot_new392_half" if args.pilot else "sp_temporal_cal")

    rows = []
    for g, f in feats.items():
        rows += tier_budgets(f, g, args.alphas, args.delta, B if g == "centroid" else min(B, 200), args.seed, source)
    bdf = write_rows(rows, out / "E6_tier_budget.csv")
    share = feats["centroid"].pred_tier.value_counts(normalize=True).to_dict()
    adf = accrual(share)
    adf.to_csv(out / "E6_accrual.csv", index=False)
    rate = dict(zip(adf.group, adf.per_year))
    bdf["years_of_accrual_for_deficit"] = [
        (r.deficit / rate.get(f"pred_tier:{r.tier}", np.nan)) if np.isfinite(r.deficit) and r.deficit > 0 else
        (0.0 if np.isfinite(r.deficit) else np.inf) for r in bdf.itertuples()]
    bdf["years_at_deficit_ci_lo"] = [
        (r.deficit_ci_lo / rate.get(f"pred_tier:{r.tier}", np.nan)) if np.isfinite(r.deficit_ci_lo) and
        r.deficit_ci_lo > 0 else (0.0 if np.isfinite(r.deficit_ci_lo) else np.nan) for r in bdf.itertuples()]
    bdf.to_csv(out / "E6_tier_budget.csv", index=False)

    sims = []  # secondary (e): how close SGR's own threshold search comes to the oracle-coverage budget
    for r in bdf[(bdf.geometry == "centroid") & np.isfinite(bdf.n_needed)].itertuples():
        m = int(r.n_needed)
        if m > 5000:
            continue
        reps = 100 if args.smoke else 400
        sims.append({"tier": r.tier, "alpha": r.alpha, "n_budget": m, "coverage": r.best_coverage,
                     "risk": r.risk_at_best, "target_power": args.power,
                     "sgr_rate_at_budget": C.simulate_sgr(m, r.best_coverage, r.risk_at_best, r.alpha, args.delta,
                                                          len(TIERS), reps=reps, seed=args.seed),
                     "sgr_rate_at_2x_budget": C.simulate_sgr(2 * m, r.best_coverage, r.risk_at_best, r.alpha,
                                                             args.delta, len(TIERS), reps=reps, seed=args.seed)})
    sdf = write_rows(sims, out / "E6_sgr_simulation.csv")

    lines = [f"# {tag} summary (delta={args.delta}, target power {args.power}, source {source})\n"]
    if args.smoke:
        lines.append("Synthetic features: code-path test only, no scientific content.\n")
    for r in bdf[(bdf.geometry == "centroid")].itertuples():
        if r.n_available == 0:
            lines.append(f"- {r.tier} alpha={r.alpha}: no calibration proteins")
            continue
        head = f"- centroid {r.tier} alpha={r.alpha}: available {r.n_available}, top-1 acc {r.top1_accuracy:.2f}; "
        if r.status != "estimated":
            what = "not certifiable at any budget" if r.status == "risk_above_alpha" else "not estimable"
            lines.append(head + f"{r.status.replace('_', ' ')} ({what}); deficit CI [{fmt_n(r.deficit_ci_lo)}, "
                                f"{fmt_n(r.deficit_ci_hi)}] (>= {r.years_at_deficit_ci_lo:.1f} years of accrual at the "
                                f"lower bound); bootstrap share of certifiable resamples "
                                f"{r.share_resamples_certifiable:.2f}")
            continue
        lines.append(head + f"needed {fmt_n(r.n_needed)} (coverage {r.best_coverage}, selective risk "
                            f"{r.risk_at_best:.2f}); deficit {fmt_n(r.deficit)} [CI {fmt_n(r.deficit_ci_lo)}, "
                            f"{fmt_n(r.deficit_ci_hi)}]; ~{r.years_of_accrual_for_deficit:.1f} years of Swiss-Prot accrual "
                            f"(>= {r.years_at_deficit_ci_lo:.1f} at the lower bound)")
    lines.append(f"- accrual: {rate['all_enzymes']:.0f} new Swiss-Prot enzymes/year pass the benchmark filters")
    for r in sdf.itertuples():
        lines.append(f"- SGR simulation {r.tier} alpha={r.alpha}: certification rate {r.sgr_rate_at_budget:.2f} at the "
                     f"oracle budget ({r.n_budget}), {r.sgr_rate_at_2x_budget:.2f} at twice it (target {r.target_power})")
    with open(out / "E6_summary.md", "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\n{'CX10 SMOKE' if args.smoke else tag} {'FAILED' if fails else 'DONE'} in {time.time() - t0:.0f}s -> {out}")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
