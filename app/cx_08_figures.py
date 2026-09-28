"""Paper figures from the E1-E6 outputs. Each figure is written as PNG (300 dpi) + PDF, with a CSV
of exactly the plotted numbers next to it. Figures whose inputs do not exist yet are skipped."""
import argparse
import json
import warnings

warnings.filterwarnings("ignore")

from cleanx import config

config.enter_app_dir()

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

INK, INK2, MUTED, GRID, AXIS = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
# Fixed entity -> colour mapping (validated categorical palette, fixed order, never cycled).
METHOD = {"centroid": ("CLEAN centroid", "#2a78d6", "o"),
          "exemplar": ("Exemplar (1-NN)", "#eb6834", "s"),
          "homology_tophit": ("Homology top hit", "#1baf7a", "^"),
          "kmean3": ("3-NN mean", "#eda100", "D"),
          "multiproto": ("Multi-prototype", "#e87ba4", "v")}
DATASET = {"new392": "New-392", "price149": "Price-149", "sp_temporal_test": "SP-temporal test"}
TIERS = ["common", "rare", "very_rare", "singleton"]
TIER_LABEL = {"common": "common\n(>10)", "rare": "rare\n(4-10)", "very_rare": "very rare\n(2-3)",
              "singleton": "singleton\n(1)", "unseen": "unseen\n(0)"}
ID_BINS = ["<30%", "30-50%", "50-70%", "70-90%", ">=90%"]


def style():
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"], "font.size": 8,
        "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "axes.titlecolor": INK, "axes.titlesize": 9,
        "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "axes.grid.axis": "y",
        "grid.color": GRID, "grid.linewidth": 0.6, "grid.linestyle": "-", "xtick.color": MUTED,
        "ytick.color": MUTED, "xtick.labelcolor": INK2, "ytick.labelcolor": INK2, "legend.frameon": False,
        "lines.linewidth": 1.5, "lines.markersize": 5, "figure.facecolor": "white",
        "axes.facecolor": "white", "savefig.facecolor": "white"})


def save(fig, name, table):
    config.CX_FIGS.mkdir(parents=True, exist_ok=True)
    fig.savefig(config.CX_FIGS / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(config.CX_FIGS / f"{name}.pdf", bbox_inches="tight")
    table.to_csv(config.CX_FIGS / f"{name}.csv", index=False)
    plt.close(fig)
    print(f"  wrote {name}.png/.pdf/.csv")


def fig_blind_spot():
    p = config.CX_RESULTS / "E1" / "E1_blind_spot.csv"
    if not p.is_file():
        return print("  skip fig1: run cx_01_audit.py")
    d = pd.read_csv(p)
    d = d[d.method == "clean_maxsep"].copy()
    if d.empty:
        return print("  skip fig1: no clean_maxsep rows in E1_blind_spot.csv")
    d["label"] = d.dataset.map(DATASET).fillna(d.dataset)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.0, 2.4), gridspec_kw={"wspace": 0.45})
    y = np.arange(len(d))
    h = 0.36
    a1.barh(y + h / 2 + 0.01, d["weighted_P(reported)"], height=h, color="#2a78d6", label="Reported (support-weighted)")
    a1.barh(y - h / 2 - 0.01, d["micro_P"], height=h, color="#eb6834", label="All predictions counted (micro)")
    for yy, v1, v2 in zip(y, d["weighted_P(reported)"], d["micro_P"]):
        a1.text(v1 + 0.01, yy + h / 2, f"{v1:.3f}", va="center", color=INK2, fontsize=7)
        a1.text(v2 + 0.01, yy - h / 2, f"{v2:.3f}", va="center", color=INK2, fontsize=7)
    a1.set_yticks(y, d.label)
    a1.set_xlim(0, 0.75)
    a1.set_xlabel("Precision at EC level 4")
    a1.grid(axis="x"); a1.grid(axis="y", visible=False)
    a1.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=1, fontsize=7)
    vis = d.FP - d.FP_outside_true_labelset
    a2.barh(y, vis, height=0.5, color="#2a78d6", label="Counted by reported precision")
    a2.barh(y, d.FP_outside_true_labelset, left=vis + 0.8, height=0.5, color="#eb6834",
            label="Invisible (EC absent from the test labels)")
    for yy, v, inv, share in zip(y, vis, d.FP_outside_true_labelset, d.share_invisible_to_weighted_P):
        a2.text(v + inv + 3, yy, f"{share:.0%} invisible", va="center", color=INK2, fontsize=7)
    a2.set_yticks(y, d.label)
    a2.set_xlabel("False-positive EC calls")
    a2.set_xlim(0, (d.FP.max()) * 1.45)
    a2.grid(axis="x"); a2.grid(axis="y", visible=False)
    a2.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=1, fontsize=7)
    save(fig, "fig1_blind_spot", d[["dataset", "weighted_P(reported)", "micro_P", "FP", "FP_outside_true_labelset",
                                    "share_invisible_to_weighted_P"]])


def _e2(backbone="split100"):
    p = config.CX_RESULTS / f"E2_{backbone}" / "E2_metrics_long.csv"
    return pd.read_csv(p) if p.is_file() else None


def fig_tiers(backbone="split100"):
    m = _e2(backbone)
    if m is None:
        return print("  skip fig2: run cx_05_geometry.py")
    t = m[(m.subset_type == "true_tier") & (m.depth == 4) & m.method.isin(["centroid", "exemplar", "homology_tophit"])]
    sets = [s for s in DATASET if s in set(t.dataset)]
    if not sets:
        return print("  skip fig2: no known dataset in E2 results")
    fig, axes = plt.subplots(1, len(sets), figsize=(2.5 * len(sets), 2.5), sharey=True, squeeze=False)
    rows = []
    for ax, s in zip(axes[0], sets):
        g = t[t.dataset == s]
        tiers = [x for x in TIERS + ["unseen"] if x in set(g.subset)]
        meths = [x for x in ["centroid", "exemplar", "homology_tophit"] if x in set(g.method)]
        w = 0.8 / len(meths)
        for k, me in enumerate(meths):
            v = [g[(g.subset == ti) & (g.method == me)].micro_F1.mean() for ti in tiers]
            ax.bar(np.arange(len(tiers)) + (k - (len(meths) - 1) / 2) * w, v, width=w - 0.03,
                   color=METHOD[me][1], label=METHOD[me][0])
            rows += [{"dataset": s, "tier": ti, "method": me, "micro_F1": vv} for ti, vv in zip(tiers, v)]
        ns = [int(g[(g.subset == ti)].n.iloc[0]) for ti in tiers]
        ax.set_xticks(np.arange(len(tiers)), [f"{TIER_LABEL[ti]}\nn={n}" for ti, n in zip(tiers, ns)], fontsize=6.5)
        ax.set_title(DATASET[s])
        ax.set_ylim(0, 1)
    axes[0][0].set_ylabel("Exact-EC micro-F1")
    axes[0][0].legend(loc="lower left", bbox_to_anchor=(0, 1.08), ncol=3, fontsize=7)
    save(fig, "fig2_tiers", pd.DataFrame(rows))


def fig_identity(backbone="split100"):
    m = _e2(backbone)
    if m is None or "identity" not in set(m.subset_type):
        return print("  skip fig3: needs cx_04_identity.py then cx_05_geometry.py")
    t = m[(m.subset_type == "identity") & (m.depth == 4)]
    sets = [s for s in DATASET if s in set(t.dataset)]
    if not sets:
        return print("  skip fig3: no known dataset with identity bins in E2 results")
    fig, axes = plt.subplots(1, len(sets), figsize=(2.5 * len(sets), 2.4), sharey=True, squeeze=False)
    rows = []
    for ax, s in zip(axes[0], sets):
        g = t[t.dataset == s]
        bins = [b for b in ID_BINS if b in set(g.subset)]
        for me in ["centroid", "exemplar", "homology_tophit"]:
            if me not in set(g.method):
                continue
            v = [g[(g.subset == b) & (g.method == me)].micro_F1.mean() for b in bins]
            ax.plot(range(len(bins)), v, color=METHOD[me][1], marker=METHOD[me][2], label=METHOD[me][0])
            rows += [{"dataset": s, "identity_bin": b, "method": me, "micro_F1": vv} for b, vv in zip(bins, v)]
        ns = [int(g[g.subset == b].n.iloc[0]) for b in bins]
        ax.set_xticks(range(len(bins)), [f"{b}\nn={n}" for b, n in zip(bins, ns)], fontsize=6.5)
        ax.set_title(DATASET[s])
        ax.set_ylim(0, 1)
        ax.set_xlabel("Max identity to training")
    axes[0][0].set_ylabel("Exact-EC micro-F1")
    axes[0][0].legend(loc="lower left", bbox_to_anchor=(0, 1.08), ncol=3, fontsize=7)
    save(fig, "fig3_identity", pd.DataFrame(rows))


def fig_ablation(backbone="split100"):
    m = _e2(backbone)
    if m is None:
        return print("  skip figS1: run cx_05_geometry.py")
    rows = []
    sets = [s for s in DATASET if s in set(m.dataset)]
    if not sets:
        return print("  skip figS1: no known dataset in E2 results")
    fig, axes = plt.subplots(1, len(sets), figsize=(2.5 * len(sets), 2.4), sharey=True, squeeze=False)
    for ax, s in zip(axes[0], sets):
        g = m[(m.dataset == s) & (m.subset_type == ("closed" if "closed" in set(m[m.dataset == s].subset_type) else "all"))]
        for me in ["centroid", "multiproto", "kmean3", "exemplar", "homology_tophit"]:
            if me not in set(g.method):
                continue
            v = [g[(g.depth == dd) & (g.method == me)].micro_F1.mean() for dd in (1, 2, 3, 4)]
            ax.plot([1, 2, 3, 4], v, color=METHOD[me][1], marker=METHOD[me][2], label=METHOD[me][0])
            rows += [{"dataset": s, "depth": dd, "method": me, "micro_F1": vv} for dd, vv in zip((1, 2, 3, 4), v)]
        ax.set_xticks([1, 2, 3, 4])
        ax.set_xlabel("EC digits required")
        ax.set_title(DATASET[s])
        ax.set_ylim(0, 1)
    axes[0][0].set_ylabel("Micro-F1 (seen-EC proteins)")
    axes[0][0].legend(loc="lower left", bbox_to_anchor=(0, 1.08), ncol=3, fontsize=7)
    save(fig, "figS1_depth_ablation", pd.DataFrame(rows))


def fig_nise(backbone="split100"):
    p = config.CX_RESULTS / f"E3_{backbone}" / "E3_stats.json"
    g = config.CX_RESULTS / f"E3_{backbone}" / "E3_ec_geometry.csv"
    if not (p.is_file() and g.is_file()):
        return print("  skip fig4: run cx_06_ec_biology.py")
    st = json.loads(p.read_text(encoding="utf-8"))
    link = st.get("exemplar_gain_by_true_EC_nise_status", {})
    if not link:
        return print("  skip fig4: E3 has no test-set link yet (run cx_05_geometry.py before cx_06_ec_biology.py)")
    rows = []
    for s, v in link.items():
        rows.append({"dataset": s, "group": "NISE-like EC", "mean_gain": v["gain_nise"], "n": v["n_nise"]})
        rows.append({"dataset": s, "group": "other EC", "mean_gain": v["gain_other"], "n": v["n_other"]})
    d = pd.DataFrame(rows).dropna()
    both = [s for s in DATASET if s in set(d.dataset) and d[(d.dataset == s) & (d.n > 0)].group.nunique() == 2]
    if not both:
        return print("  skip fig4: no dataset has proteins of both NISE-like and other ECs")
    d = d[d.dataset.isin(both)]
    fig, ax = plt.subplots(figsize=(1.6 + 1.3 * len(both), 2.4))
    x = np.arange(len(both))
    for k, (grp, col) in enumerate((("NISE-like EC", "#eb6834"), ("other EC", "#2a78d6"))):
        v = [d[(d.dataset == s) & (d.group == grp)].mean_gain.mean() for s in both]
        ax.bar(x + (k - 0.5) * 0.38, v, width=0.35, color=col, label=grp)
    ax.axhline(0, color=AXIS, linewidth=0.8)
    ticks = []
    for s in both:
        n1 = int(d[(d.dataset == s) & (d.group == "NISE-like EC")].n.sum())
        n2 = int(d[(d.dataset == s) & (d.group == "other EC")].n.sum())
        ticks.append(f"{DATASET[s]}\nn={n1} / {n2}")
    ax.set_xticks(x, ticks, fontsize=7)
    ax.set_ylabel("Per-protein F1 gain\n(exemplar - centroid)")
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.02), ncol=2, fontsize=7)
    save(fig, "fig4_nise_gain", d)


def fig_reliability():
    base = config.CX_RESULTS / "E4"
    need = [base / f"E4_features_{n}_{g}.csv" for n in ("cal", "sp_temporal_test") for g in ("centroid", "exemplar")]
    if not all(p.is_file() for p in need):
        return print("  skip fig5/6: run cx_07_reliability.py")
    from cleanx.calibration import MondrianIso, reliability, risk_coverage
    cal_e = pd.read_csv(need[1])
    te_c, te_e = pd.read_csv(need[2]), pd.read_csv(need[3])
    mond = MondrianIso(50).fit(cal_e.score, cal_e.correct4, cal_e.pred_tier)
    te_e["conf"] = mond.predict(te_e.score, te_e.pred_tier)
    tiers = [t for t in TIERS if (te_c.pred_tier == t).sum() >= 10]
    if not tiers:
        return print("  skip fig5-7: no predicted-EC tier has >= 10 test proteins")
    fig, axes = plt.subplots(1, len(tiers), figsize=(2.1 * len(tiers), 2.3), sharey=True, squeeze=False)
    rows = []
    for ax, t in zip(axes[0], tiers):
        for lab, f, conf, col, mk in (("CLEAN GMM confidence (released)", te_c, te_c.clean_gmm, "#2a78d6", "o"),
                                      ("Exemplar + tier-conditional isotonic", te_e, te_e.conf, "#eb6834", "s")):
            msk = (f.pred_tier == t).to_numpy()
            bins, ece = reliability(np.asarray(conf)[msk], f.correct4.to_numpy()[msk], n_bins=5, strategy="quantile")
            ax.plot([b["mean_conf"] for b in bins], [b["acc"] for b in bins], color=col, marker=mk,
                    label=lab)
            rows += [{"tier": t, "score": lab, "ECE": ece, **b} for b in bins]
        ax.plot([0, 1], [0, 1], color=MUTED, linewidth=0.8)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.set_title(f"predicted EC {t.replace('_', ' ')}\n(n={int((te_c.pred_tier == t).sum())} / "
                     f"{int((te_e.pred_tier == t).sum())})")
        ax.set_xlabel("Mean confidence")
        ax.grid(axis="x")
    axes[0][0].set_ylabel("Observed accuracy")
    axes[0][0].legend(loc="lower left", bbox_to_anchor=(0, 1.2), ncol=2, fontsize=7)
    save(fig, "fig5_reliability_by_tier", pd.DataFrame(rows))

    fig, ax = plt.subplots(figsize=(3.4, 2.6))
    rows = []
    for lab, f, conf, col in (("CLEAN centroid, GMM confidence", te_c, te_c.clean_gmm, "#2a78d6"),
                              ("Exemplar, tier-conditional isotonic", te_e, te_e.conf, "#eb6834")):
        cov, risk, aurc = risk_coverage(np.asarray(conf), f.correct4.to_numpy())
        ax.plot(cov, risk, color=col, label=f"{lab} (AURC {aurc:.3f})")
        rows += [{"score": lab, "coverage": c, "selective_risk": r} for c, r in zip(cov[::5], risk[::5])]
    ax.set_xlabel("Coverage (share of proteins annotated)")
    ax.set_ylabel("Selective risk (exact EC wrong)")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.grid(axis="x")
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), fontsize=7)
    save(fig, "fig6_risk_coverage", pd.DataFrame(rows))

    rc = base / "E4_risk_control.csv"
    if rc.is_file():
        r = pd.read_csv(rc)
        r = r[(r.dataset == "sp_temporal_test") & (r.group.str.startswith("pred_tier")) & (r.geometry == "exemplar")]
        alpha = 0.2 if 0.2 in set(r.alpha) else sorted(set(r.alpha))[0]
        r = r[r.alpha == alpha]
        if len(r):
            fig, ax = plt.subplots(figsize=(3.6, 2.4))
            tiers = [t for t in TIERS if f"pred_tier:{t}" in set(r.group)]
            for k, (rule, col) in enumerate((("global", "#2a78d6"), ("mondrian", "#eb6834"))):
                v = [r[(r.rule == rule) & (r.group == f"pred_tier:{t}")].selective_risk.mean() for t in tiers]
                cv = [r[(r.rule == rule) & (r.group == f"pred_tier:{t}")].coverage.mean() for t in tiers]
                xs = np.arange(len(tiers)) + (k - 0.5) * 0.38
                ax.bar(xs, np.nan_to_num(v), width=0.35, color=col,
                       label="one global threshold" if rule == "global" else "one threshold per tier")
                for xi, vi, ci_ in zip(xs, v, cv):
                    ax.text(xi, (0 if np.isnan(vi) else vi) + 0.01, "abstain" if np.isnan(vi) else f"cov {ci_:.0%}",
                            ha="center", va="bottom", fontsize=6, color=INK2, rotation=90)
            ax.axhline(alpha, color=MUTED, linewidth=0.8)
            ax.text(len(tiers) - 0.5, alpha + 0.01, f"target risk {alpha:g}", ha="right", fontsize=6.5, color=INK2)
            ax.set_xticks(range(len(tiers)), [TIER_LABEL[t] for t in tiers], fontsize=6.5)
            ax.set_ylabel("Realized selective risk (test)")
            ax.set_ylim(0, 1)
            ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2, fontsize=7)
            save(fig, "fig7_tier_risk_control", r)


def _pick(files, *tags):
    """First results folder among tags that holds every file in files (a partial or unrelated folder is skipped)."""
    for t in tags:
        if all((config.CX_RESULTS / t / f).is_file() for f in files):
            return t
    return None


def fig_ecbench_rescoring():
    """E5: reported (support-weighted) vs micro F1 for EC-Bench's 16 methods, and inflation vs ECs predicted."""
    tag = _pick(("E5_cards.csv",), "E5", "E5_pilot")
    if tag is None:
        return print("  skip fig8: run cx_09_evalcard.py (or --pilot)")
    c = pd.read_csv(config.CX_RESULTS / tag / "E5_cards.csv", dtype={"regime": str})
    ds = "test" if tag == "E5" else "price149"
    g = c[(c.source == "ecbench") & (c.dataset == ds) & (c.regime == "100") & (c.depth == 4) &
          (c.variant == "as_published")].copy()  # all methods, as in the pre-registered S4
    if g.empty:
        return print(f"  skip fig8: no EC-Bench rows for {ds} in {tag}")
    g = g.sort_values("reported_ECBench_F1")
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.0, 3.0), gridspec_kw={"wspace": 0.55, "width_ratios": [1.15, 1]})
    y = np.arange(len(g))
    for yy, r, m in zip(y, g.reported_ECBench_F1, g.u_micro_F1):
        a1.plot([r, m], [yy, yy], color=AXIS, linewidth=1.0, zorder=1)
    a1.scatter(g.reported_ECBench_F1, y, s=22, color="#2a78d6", marker="o", zorder=3,
               edgecolors="white", linewidths=0.6, label="Reported (support-weighted)")
    a1.scatter(g.u_micro_F1, y, s=22, color="#eb6834", marker="s", zorder=3,
               edgecolors="white", linewidths=0.6, label="Every FP counted (micro)")
    a1.set_yticks(y, g.method, fontsize=6.5)
    a1.set_xlabel("F1 at EC level 4")
    a1.set_xlim(0, max(0.05, float(max(g.reported_ECBench_F1.max(), g.u_micro_F1.max()))) * 1.12)
    a1.grid(axis="x"); a1.grid(axis="y", visible=False)
    a1.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=1, fontsize=6.5)
    infl = g.reported_ECBench_F1 - g.u_micro_F1
    a2.axhline(0, color=AXIS, linewidth=0.8)
    a2.scatter(g.mean_pred_labels, infl, s=24, color="#2a78d6", marker="o", edgecolors="white", linewidths=0.6)
    for _, r in g.assign(infl=infl).sort_values("infl", ascending=False).head(3).iterrows():
        a2.annotate(r.method, (r.mean_pred_labels, r.infl), xytext=(4, 2), textcoords="offset points",
                    fontsize=6, color=INK2)
    from scipy.stats import spearmanr
    rho = spearmanr(g.mean_pred_labels, infl)[0]
    a2.text(0.02, 0.97, f"Spearman rho = {rho:.2f} (n = {len(g)} methods)", transform=a2.transAxes,
            va="top", fontsize=6.5, color=INK2)
    a2.set_xlabel("Mean ECs predicted per protein")
    a2.set_ylabel("Reported F1 - micro F1")
    title = "EC-Bench test set" if ds == "test" else "EC-Bench Price-149 (pilot)"
    fig.suptitle(f"{title}, models trained on the 100% clustering", fontsize=8, color=INK, y=1.08)
    name = "fig8_ecbench_rescoring" if tag == "E5" else "figS2_ecbench_rescoring_pilot"
    save(fig, name, g.assign(inflation=infl)[["method", "reported_ECBench_F1", "u_micro_F1", "u_macro_union_F1",
                                              "mean_pred_labels", "inflation", "fp", "fp_invisible"]])


def fig_cert_budget():
    """E6: calibration proteins needed to certify selective risk, universally and per predicted tier."""
    tag = _pick(("E6_nstar_grid.csv", "E6_tier_budget.csv"), "E6", "E6_pilot")
    if tag is None:
        return print("  skip fig9: run cx_10_certbudget.py (or --pilot)")
    grid = pd.read_csv(config.CX_RESULTS / tag / "E6_nstar_grid.csv")
    bud = pd.read_csv(config.CX_RESULTS / tag / "E6_tier_budget.csv")
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.0, 2.7), gridspec_kw={"wspace": 0.4})
    styles = [(0.1, "#2a78d6", "o"), (0.2, "#eb6834", "s"), (0.3, "#1baf7a", "^")]
    rows = []
    for a, col, mk in styles:
        g = grid[(np.isclose(grid.alpha, a)) & np.isfinite(grid.n_stable)]
        if g.empty:
            continue
        a1.plot(g.true_selective_risk, g.n_stable, color=col, marker=mk, markevery=max(1, len(g) // 6),
                markersize=4, label=f"target risk {a:g}")
        last = g.iloc[-1]
        a1.annotate(f"{a:g}", (last.true_selective_risk, last.n_stable), xytext=(3, 0), textcoords="offset points",
                    fontsize=6.5, color=INK2, va="center")
        rows += [{"panel": "a", "alpha": a, "true_selective_risk": r.true_selective_risk, "n_needed": r.n_stable}
                 for r in g.itertuples()]
    a1.set_yscale("log")
    a1.set_xlabel("True selective risk of the accepted calls")
    a1.set_ylabel("Accepted calibration proteins\nneeded (power 0.8)")
    a1.grid(axis="x")
    if len(grid):
        a1.text(0.98, 0.03, f"per-test delta {grid.per_test_delta.iloc[0]:.4f}", transform=a1.transAxes,
                ha="right", va="bottom", fontsize=6, color=MUTED)
    a1.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=3, fontsize=6.5)
    b = bud[(bud.geometry == "centroid") & np.isclose(bud.alpha, 0.2)].set_index("tier").reindex(TIERS)
    x = np.arange(len(TIERS))
    need = b.n_needed.to_numpy(float)
    have = b.n_available.to_numpy(float)
    top = np.nanmax(np.concatenate([need[np.isfinite(need)], have[np.isfinite(have)], [10.0]])) * 4
    shown = np.where(np.isfinite(need), need, top)
    for xi, n_, h_ in zip(x, shown, have):
        a2.plot([xi, xi], [max(h_, 1), n_], color=AXIS, linewidth=1.0, zorder=1)
    a2.scatter(x, np.maximum(have, 1), s=26, color="#eb6834", marker="s", zorder=3, edgecolors="white",
               linewidths=0.6, label="available (calibration split)")
    fin = np.isfinite(need)
    a2.scatter(x[fin], shown[fin], s=26, color="#2a78d6", marker="o", zorder=3, edgecolors="white", linewidths=0.6,
               label="needed (lower bound)")
    if (~fin).any():  # unbounded: an open arrow at the top, never a dot that reads as a finite value
        a2.scatter(x[~fin], shown[~fin], s=34, facecolors="none", edgecolors="#2a78d6", marker="^", zorder=3,
                   linewidths=1.0, label="no certifiable operating point")
    for xi, n_, yrs in zip(x, need, b.years_of_accrual_for_deficit):
        txt = "not certifiable" if not np.isfinite(n_) else (f"~{yrs:.0f} y" if np.isfinite(yrs) and yrs >= 1 else "")
        if txt:
            a2.annotate(txt, (xi, n_ if np.isfinite(n_) else top), xytext=(0, 5), textcoords="offset points",
                        ha="center", fontsize=6, color=INK2)
    a2.set_yscale("log")
    a2.set_xticks(x, [TIER_LABEL[t] for t in TIERS], fontsize=6.5)
    a2.set_ylabel("Calibration proteins\n(target risk 0.2)")
    a2.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=1, fontsize=6.5)
    rows += [{"panel": "b", "tier": t, "n_needed": n_, "n_available": h_, "years": y_}
             for t, n_, h_, y_ in zip(TIERS, need, have, b.years_of_accrual_for_deficit)]
    name = "fig9_certification_budget" if tag == "E6" else "figS3_certification_budget_pilot"
    save(fig, name, pd.DataFrame(rows))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--backbone", default="split100")
    args = ap.parse_args()
    style()
    fig_blind_spot()
    fig_tiers(args.backbone)
    fig_identity(args.backbone)
    fig_ablation(args.backbone)
    fig_nise(args.backbone)
    fig_reliability()
    fig_ecbench_rescoring()
    fig_cert_budget()
    print(f"CX08 DONE -> {config.CX_FIGS}")


if __name__ == "__main__":
    main()
