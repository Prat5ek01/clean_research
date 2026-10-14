"""Render paper figures as PNG, PDF and source CSV.

Skip missing inputs or F1 estimates that disagree with the stored metrics."""
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
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

from cleanx.metrics import bootstrap, ci, evaluate

INK, INK2, MUTED, GRID, AXIS = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
# Fixed method colours; colour-vision contrast checked on 2026-09-28.
METHOD = {"centroid": ("CLEAN centroid", "#2a78d6", "o"),
          "exemplar": ("Exemplar (1-NN)", "#eb6834", "s"),
          "homology_tophit": ("Homology top hit", "#1baf7a", "^"),
          "kmean3": ("3-NN mean", "#eda100", "D"),
          "multiproto": ("Multi-prototype", "#e87ba4", "v"),
          "balanced3": ("Size-balanced 3-NN", "#4a3aa7", "P")}
FIG_METHODS = ["centroid", "balanced3", "exemplar", "homology_tophit"]
DATASET = {"new392": "New-392", "price149": "Price-149", "sp_temporal_test": "SP-temporal test"}
TIERS = ["common", "rare", "very_rare", "singleton"]
TIER_LABEL = {"common": "common\n(>10)", "rare": "rare\n(4-10)", "very_rare": "very rare\n(2-3)",
              "singleton": "singleton\n(1)", "unseen": "unseen\n(0)"}
ID_BINS = ["<30%", "30-50%", "50-70%", "70-90%", ">=90%"]
B_FIG, SEED = 2000, config.SEED  # the bootstrap protocol of cx_05
MIN_N = 10  # smaller groups are drawn without an interval and marked
MONDRIAN_MIN_N = 50  # cx_07 default (--mondrian-min-n)
TOOL = None  # tool behind the identity bins and the homology baseline, set in main()


def style():
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"], "font.size": 8,
        "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "axes.titlecolor": INK, "axes.titlesize": 9,
        "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "axes.grid.axis": "y",
        "grid.color": GRID, "grid.linewidth": 0.6, "grid.linestyle": "-", "xtick.color": MUTED,
        "ytick.color": MUTED, "xtick.labelcolor": INK2, "ytick.labelcolor": INK2, "legend.frameon": False,
        "lines.linewidth": 1.5, "lines.markersize": 5, "figure.facecolor": "white",
        "axes.facecolor": "white", "savefig.facecolor": "white"})


def note(fig, text):
    """A caption 6 pt under everything drawn in the axes (tick labels, axis labels, axes legends), wrapped to the
    figure width so that it never widens the saved figure beyond a journal page. The position is resolved when the
    figure is drawn, by the renderer that saves it, so it holds for any font (Arial on Windows, DejaVu elsewhere)."""
    import textwrap
    from matplotlib.text import Annotation
    from matplotlib.transforms import Bbox
    width = int(fig.get_figwidth() * 72 / 3.4)  # characters per line at 6 pt (about 3.4 pt per character)
    text = "\n".join(textwrap.fill(par, width) for par in text.split("\n"))
    axes = list(fig.axes)
    fig.add_artist(Annotation(text, xy=(0.5, 0), xycoords=lambda r: Bbox.union([a.get_tightbbox(r) for a in axes]),
                              xytext=(0, -6), textcoords="offset points", ha="center", va="top", fontsize=6,
                              color=INK2, annotation_clip=False))


def save(fig, name, table):
    config.CX_FIGS.mkdir(parents=True, exist_ok=True)
    fig.savefig(config.CX_FIGS / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(config.CX_FIGS / f"{name}.pdf", bbox_inches="tight")
    table.to_csv(config.CX_FIGS / f"{name}.csv", index=False)
    plt.close(fig)
    print(f"  wrote {name}.png/.pdf/.csv")


def identity_tool():
    """Tool that produced the identity bins and the homology baseline (cx_04 writes both in one run)."""
    tools = set()
    for s in DATASET:
        p = config.CX_IDENTITY / f"{s}.csv"
        if p.is_file():
            df = pd.read_csv(p, dtype=str, keep_default_na=False)
            tools |= set(df["tool"]) if "tool" in df else {"unrecorded"}
    if len(tools) > 1:
        print(f"  WARNING: the identity tables come from different tools {sorted(tools)}; re-run cx_04_identity.py")
    return tools.pop() if len(tools) == 1 else None


def mlabel(m):
    if m == "homology_tophit" and TOOL == "blast":
        return "BLASTp top hit"
    if m == "balanced1":
        return "Size-balanced 1-NN"
    return METHOD[m][0] if m in METHOD else m


def mcolor(m):
    return METHOD[m][1] if m in METHOD else INK2


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


def _deltas(backbone="split100"):
    p = config.CX_RESULTS / f"E2_{backbone}" / "E2_deltas.csv"
    return pd.read_csv(p) if p.is_file() else None


def _per_query(backbone="split100"):
    p = config.CX_RESULTS / f"E2_{backbone}" / "E2_per_query.csv"
    return pd.read_csv(p, dtype=str, keep_default_na=False) if p.is_file() else None


def _methods_in(q, wanted):
    """Methods scored on every protein of this dataset (a method missing for a dataset has empty f1 cells)."""
    return [m for m in wanted if f"f1_{m}" in q and (q[f"f1_{m}"] != "").all()]


def _scored(q, methods):
    """Truth and per-method predictions of one dataset, exactly as cx_05 scored them."""
    split = lambda col: [[e for e in s.split(";") if e] for s in col]
    ids = list(q.Entry)
    return dict(zip(ids, split(q.true_EC))), {m: dict(zip(ids, split(q[f"pred_{m}"]))) for m in methods}


def _f1_with_ci(groups, truth, preds):
    """Exact-EC micro-F1 and its paired-bootstrap 95% CI per group and method."""
    rows = []
    for g, sub in groups.items():
        bs = bootstrap(truth, preds, sub, depth=4, B=B_FIG, seed=SEED, keys=("micro_F1",))
        for m, P in preds.items():
            lo, hi = ci(bs[m]["micro_F1"])
            rows.append({"group": g, "n": len(sub), "method": m, "micro_F1": evaluate(truth, P, sub, 4)["micro_F1"],
                         "ci_lo": lo, "ci_hi": hi})
    return pd.DataFrame(rows)


def _mismatch(rows, m_long, dataset, subset_type):
    ref = m_long[(m_long.dataset == dataset) & (m_long.subset_type == subset_type) & (m_long.depth == 4)]
    ref = {(r.subset, r.method): r.micro_F1 for r in ref.itertuples()}
    return [(r.group, r.method) for r in rows.itertuples()
            if not abs(ref.get((r.group, r.method), np.nan) - r.micro_F1) <= 1e-9]


def _holm_p(dl, dataset, tier, method):
    """Holm-adjusted p of the pre-registered per-tier test (H3a) of method vs the centroid, or None."""
    if dl is None or "p_holm" not in dl:
        return None
    r = dl[(dl.dataset == dataset) & (dl.subset == f"tier:{tier}") & (dl.comparison == f"{method} - centroid") &
           (dl.metric == "micro_F1")]
    return float(r.p_holm.iloc[0]) if len(r) and np.isfinite(r.p_holm.iloc[0]) else None


def fig_tiers(backbone="split100"):
    q, m_long = _per_query(backbone), _e2(backbone)
    if q is None or m_long is None:
        return print("  skip fig2: run cx_05_geometry.py")
    sets = [s for s in DATASET if s in set(q.dataset)]
    if not sets:
        return print("  skip fig2: no known dataset in E2 results")
    dl = _deltas(backbone)
    per, unseen = {}, []
    for s in sets:
        g = q[q.dataset == s]
        meths = _methods_in(g, FIG_METHODS)
        truth, preds = _scored(g, meths)
        tiers = [t for t in TIERS + ["unseen"] if t in set(g.tier)]
        rows = _f1_with_ci({t: list(g.Entry[g.tier == t]) for t in tiers}, truth, preds)
        bad = _mismatch(rows, m_long, s, "true_tier")
        if bad:
            return print(f"  skip fig2: E2_per_query.csv does not reproduce E2_metrics_long.csv ({s}: {bad[:3]})")
        # unseen-EC proteins can only be right through a second, seen EC: reported in the note and the CSV
        if "unseen" in tiers:
            u = rows[rows.group == "unseen"]
            unseen.append(f"{DATASET[s]} n={int(u.n.iloc[0])}: {u.micro_F1.min():.2f}-{u.micro_F1.max():.2f}")
        per[s] = ([t for t in tiers if t != "unseen"], meths, rows.assign(dataset=s))
    fig, axes = plt.subplots(1, len(sets), figsize=(7.2, 2.6), sharey=True, squeeze=False,
                             gridspec_kw={"width_ratios": [len(per[s][0]) for s in sets], "wspace": 0.06})
    out, handles = [], {}
    for ax, s in zip(axes[0], sets):
        tiers, meths, rows = per[s]
        out += [{**rr, "tier": rr["group"], "p_holm_vs_centroid": None, "plotted": False}
                for rr in rows[rows.group == "unseen"].to_dict("records")]
        w = 0.84 / len(meths)
        for k, me in enumerate(meths):
            r = rows[rows.method == me].set_index("group").reindex(tiers)
            xs = np.arange(len(tiers)) + (k - (len(meths) - 1) / 2) * w
            handles[me] = ax.bar(xs, r.micro_F1, width=w - 0.02, color=METHOD[me][1], label=mlabel(me))
            ax.errorbar(xs, r.micro_F1, yerr=[r.micro_F1 - r.ci_lo, r.ci_hi - r.micro_F1], fmt="none", ecolor=INK,
                        elinewidth=0.6, capsize=0)
            ph = [_holm_p(dl, s, t, me) if me != "centroid" else None for t in tiers]
            for x_, top, p in zip(xs, r.ci_hi, ph):
                if p is not None and p < 0.05:
                    ax.text(x_, top + 0.01, "*", ha="center", va="bottom", fontsize=8, color=INK)
            out += [{**rr, "tier": t, "p_holm_vs_centroid": p, "plotted": True}
                    for t, rr, p in zip(tiers, r.reset_index().to_dict("records"), ph)]
        ns = [int(rows[rows.group == t].n.iloc[0]) for t in tiers]
        ax.set_xticks(np.arange(len(tiers)), [f"{TIER_LABEL[t]}\nn={n}" for t, n in zip(tiers, ns)], fontsize=6)
        ax.set_title(DATASET[s])
        ax.set_ylim(0, 1)
        ax.set_xlim(-0.55, len(tiers) - 0.45)
    axes[0][0].set_ylabel("Exact-EC micro-F1")
    fig.legend([handles[m] for m in FIG_METHODS if m in handles], [mlabel(m) for m in FIG_METHODS if m in handles],
               loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=4, fontsize=7)
    note(fig, "Tier: training proteins of the protein's rarest true EC. Whiskers: 95% bootstrap CI. * differs from "
         "CLEAN centroid (paired bootstrap, Holm-adjusted p < 0.05; the pre-registered per-tier tests)." +
         (f"\nNot shown, proteins with an EC absent from training (micro-F1 range over methods): {'; '.join(unseen)}."
          if unseen else ""))
    save(fig, "fig2_tiers", pd.DataFrame(out)[["dataset", "tier", "method", "n", "micro_F1", "ci_lo", "ci_hi",
                                               "p_holm_vs_centroid", "plotted"]])


def fig_identity(backbone="split100"):
    q, m_long = _per_query(backbone), _e2(backbone)
    if q is None or m_long is None or "identity" not in set(m_long.subset_type):
        return print("  skip fig3: needs cx_04_identity.py then cx_05_geometry.py")
    sets = [s for s in DATASET if s in set(m_long[m_long.subset_type == "identity"].dataset) and s in set(q.dataset)]
    if not sets:
        return print("  skip fig3: no known dataset with identity bins in E2 results")
    per = {}
    for s in sets:
        g = q[(q.dataset == s) & (q.tier != "unseen") & (q.identity_bin != "")]
        meths = _methods_in(q[q.dataset == s], FIG_METHODS)
        truth, preds = _scored(g, meths)
        bins = [b for b in ID_BINS if b in set(g.identity_bin)]
        rows = _f1_with_ci({b: list(g.Entry[g.identity_bin == b]) for b in bins}, truth, preds)
        bad = _mismatch(rows, m_long, s, "identity")
        if bad:
            return print(f"  skip fig3: E2_per_query.csv does not reproduce E2_metrics_long.csv ({s}: {bad[:3]})")
        per[s] = (bins, meths, rows.assign(dataset=s))
    fig, axes = plt.subplots(1, len(sets), figsize=(7.2, 2.5), sharey=True, squeeze=False,
                             gridspec_kw={"wspace": 0.08})
    out = []
    for ax, s in zip(axes[0], sets):
        bins, meths, rows = per[s]
        ns = np.array([int(rows[rows.group == b].n.iloc[0]) for b in bins])
        ok = ns >= MIN_N
        for k, me in enumerate(meths):
            r = rows[rows.method == me].set_index("group").reindex(bins)
            xs = np.arange(len(bins)) + (k - (len(meths) - 1) / 2) * 0.1
            col, mk = METHOD[me][1], METHOD[me][2]
            v = r.micro_F1.to_numpy()
            ax.plot(xs, np.where(ok, v, np.nan), color=col, marker=mk, markersize=4, label=mlabel(me))
            ax.errorbar(xs[ok], v[ok], yerr=[(r.micro_F1 - r.ci_lo)[ok], (r.ci_hi - r.micro_F1)[ok]], fmt="none",
                        ecolor=col, elinewidth=0.7, capsize=0)
            if (~ok).any():
                ax.plot(xs[~ok], v[~ok], linestyle="none", marker=mk, markersize=4, markerfacecolor="white",
                        markeredgecolor=col)
            out += [{**rr, "drawn_with_ci": bool(o)} for rr, o in zip(r.reset_index().to_dict("records"), ok)]
        ax.set_xticks(range(len(bins)), [f"{b}\nn={n}" for b, n in zip(bins, ns)], fontsize=6)
        ax.set_title(DATASET[s])
        ax.set_ylim(0, 1)
    suffix = {"blast": " (BLASTp)", "python": " (Smith-Waterman)", "m8": " (MMseqs2)"}.get(TOOL, "")
    axes[0][len(sets) // 2].set_xlabel(f"Maximum sequence identity to the training set{suffix}")
    axes[0][0].set_ylabel("Exact-EC micro-F1\n(seen-EC proteins)")
    axes[0][0].legend(loc="lower left", bbox_to_anchor=(0, 1.1), ncol=4, fontsize=7)
    note(fig, "Whiskers: 95% bootstrap CI." + (f" Open markers: bins with n < {MIN_N} (no interval, not connected)."
                                               if any((per[s][2].n < MIN_N).any() for s in sets) else ""))
    out = pd.DataFrame(out).rename(columns={"group": "identity_bin"})
    save(fig, "fig3_identity", out[["dataset", "identity_bin", "method", "n", "micro_F1", "ci_lo", "ci_hi",
                                    "drawn_with_ci"]])


def fig_primary_forest():
    """E2: the pre-registered primary endpoint with its sensitivity analyses, the descriptive geometries and the
    other test sets, as paired differences vs the CLEAN centroid with 95% bootstrap CIs."""
    d, d70 = _deltas("split100"), _deltas("split70")
    if d is None:
        return print("  skip fig10: run cx_05_geometry.py")
    T = "sp_temporal_test"
    # (label, kind, deltas table, dataset, subset, method); kind "head" is a block title
    spec = [("Primary endpoint (SP-temporal test)", "head", None, None, None, None),
            ("Size-balanced 3-NN", "primary", d, T, "closed", "balanced3"),
            ("Sensitivity analyses", "head", None, None, None, None),
            ("cluster bootstrap over ECs", "sensitivity", d, T, "closed|cluster_bootstrap_EC", "balanced3"),
            ("at most 5 proteins per EC", "sensitivity", d, T, "closed_capped5", "balanced3"),
            ("<50% identity to the CAL split", "sensitivity", d, T, "closed_cal_disjoint", "balanced3"),
            ("experimental evidence only (PE1)", "sensitivity", d, T, "closed_PE1", "balanced3"),
            ("split70 backbone", "sensitivity", d70, T, "closed", "balanced3"),
            ("Other methods (descriptive)", "head", None, None, None, None)] + \
           [(mlabel(m), "descriptive", d, T, "closed", m)
            for m in ("exemplar", "kmean3", "balanced1", "multiproto", "homology_tophit")] + \
           [("Other test sets (size-balanced 3-NN)", "head", None, None, None, None),
            ("New-392 (development set)", "other_set", d, "new392", "closed", "balanced3"),
            ("Price-149", "other_set", d, "price149", "closed", "balanced3")]
    rows, labels, kinds = [], [], []
    for lab, kind, df, ds, sub, m in spec:
        if kind == "head":
            labels.append(lab); kinds.append(kind)
            continue
        got = df[(df.dataset == ds) & (df.subset == sub) & (df.comparison == f"{m} - centroid")] \
            if df is not None else pd.DataFrame()
        if got.empty:
            continue
        for metric in ("macro_union_F1", "micro_F1"):
            r = got[got.metric == metric]
            if len(r):
                r = r.iloc[0]
                rows.append({"row": len(labels), "kind": kind, "label": lab.strip(), "dataset": ds, "subset": sub,
                             "method": m, "metric": metric, "n": int(r.n), "delta_obs": r.delta_obs,
                             "ci_lo": r.ci_lo, "ci_hi": r.ci_hi, "p_two_sided": r.p_two_sided})
        labels.append(f"{lab} (n={int(got.n.iloc[0])})"); kinds.append(kind)
    t = pd.DataFrame(rows)
    if t.empty or not (t.kind == "primary").any():
        return print("  skip fig10: no balanced3 - centroid row on SP-temporal test in E2_deltas.csv")
    heads = [k for k, kd in enumerate(kinds) if kd == "head"]
    ends = heads[1:] + [len(labels)]
    keep = sorted(set(t.row) | {h for h, e in zip(heads, ends) if any(h < r < e for r in t.row)})
    pos = {r: k for k, r in enumerate(keep)}
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 0.21 * len(keep) + 0.8), sharey=True)
    fig.subplots_adjust(left=0.36, right=0.87, wspace=0.85)  # room for the row labels and the number columns
    for ax, metric, title in zip(axes, ("macro_union_F1", "micro_F1"),
                                 ("Union-label macro-F1 (primary metric)", "Micro-F1 (key secondary)")):
        for r in t[t.metric == metric].itertuples():
            yy, col = pos[r.row], mcolor(r.method)
            ax.plot([r.ci_lo, r.ci_hi], [yy, yy], color=col, linewidth=1.2, solid_capstyle="butt")
            big = r.kind == "primary"
            ax.plot([r.delta_obs], [yy], marker="D" if big else "o", markersize=5 if big else 3.5, color=col,
                    markeredgecolor="white", markeredgewidth=0.5)
            ax.text(1.03, yy, f"{r.delta_obs:+.3f} [{r.ci_lo:.3f}, {r.ci_hi:.3f}]",
                    transform=ax.get_yaxis_transform(), va="center", ha="left", fontsize=5.5, color=INK2)
        ax.axvline(0, color=INK2, linewidth=0.8)
        ax.set_title(title, fontsize=8)
        ax.set_xlabel("Difference vs CLEAN centroid (EC4)")
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
        ax.grid(axis="x"); ax.grid(axis="y", visible=False)
    axes[0].set_yticks(range(len(keep)), [labels[r] for r in keep], fontsize=6.5)
    for k, r in enumerate(keep):
        if kinds[r] == "head":
            axes[0].get_yticklabels()[k].set_fontweight("bold")
            axes[0].get_yticklabels()[k].set_color(INK)
    axes[0].set_ylim(len(keep) - 0.5, -0.5)
    axes[0].tick_params(axis="y", length=0)
    note(fig, "Points: observed difference; bars: 95% paired bootstrap CI (B = 2000; clusters = ECs where stated). The "
         "primary row is the pre-registered confirmatory test; the first four sensitivity analyses were fixed before "
         "any test result existed; split70 (CLEAN's second released backbone) is a robustness check; all other rows "
         "are descriptive.")
    save(fig, "fig10_primary_endpoint", t.drop(columns="row"))


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
            ax.plot([1, 2, 3, 4], v, color=METHOD[me][1], marker=METHOD[me][2], label=mlabel(me))
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
        rows.append({"dataset": s, "group": "NISE-like EC", "mean_gain": v["gain_nise"], "n": v["n_nise"],
                     "mannwhitney_p": v.get("mannwhitney_p")})
        rows.append({"dataset": s, "group": "other EC", "mean_gain": v["gain_other"], "n": v["n_other"],
                     "mannwhitney_p": v.get("mannwhitney_p")})
    d = pd.DataFrame(rows).dropna(subset=["mean_gain"])
    both = [s for s in DATASET if s in set(d.dataset) and d[(d.dataset == s) & (d.n > 0)].group.nunique() == 2]
    if not both:
        return print("  skip fig4: no dataset has proteins of both NISE-like and other ECs")
    d = d[d.dataset.isin(both)]
    fig, ax = plt.subplots(figsize=(1.6 + 1.3 * len(both), 2.4))
    x = np.arange(len(both))
    # Grey denotes EC groups; method colours are reserved for other figures.
    for k, (grp, col) in enumerate((("NISE-like EC", INK2), ("other EC", MUTED))):
        v = [d[(d.dataset == s) & (d.group == grp)].mean_gain.mean() for s in both]
        ax.bar(x + (k - 0.5) * 0.38, v, width=0.35, color=col, label=grp)
    ax.axhline(0, color=AXIS, linewidth=0.8)
    ticks = []
    for s in both:
        n1 = int(d[(d.dataset == s) & (d.group == "NISE-like EC")].n.sum())
        n2 = int(d[(d.dataset == s) & (d.group == "other EC")].n.sum())
        pv = d[d.dataset == s].mannwhitney_p.iloc[0]
        ticks.append(f"{DATASET[s]}\nn={n1} / {n2}" + (f"\nMann-Whitney p={pv:.2f}" if pv == pv and pv is not None
                                                       else ""))
    ax.set_xticks(x, ticks, fontsize=7)
    ax.set_ylabel("Per-protein F1 gain\n(exemplar - centroid)")
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.02), ncol=2, fontsize=7)
    save(fig, "fig4_nise_gain", d)


def _acc_ci(bins):
    """Exact 95% CI of each reliability bin's accuracy."""
    from cleanx.calibration import risk_ci
    out = []
    for b in bins:
        k = int(round(b["acc"] * b["n"]))
        lo, hi = risk_ci(b["n"] - k, b["n"])
        out.append((1 - hi, 1 - lo))
    return out


def fig_reliability():
    base = config.CX_RESULTS / "E4"
    feat = {(n, g): base / f"E4_features_{n}_{g}.csv" for n in ("cal", "sp_temporal_test")
            for g in ("centroid", "exemplar", "balanced3")}
    if not (feat["cal", "centroid"].is_file() and feat["sp_temporal_test", "centroid"].is_file()):
        return print("  skip fig5/6: run cx_07_reliability.py")
    from cleanx.calibration import MondrianIso, reliability, risk_coverage
    conf, te = {}, {}
    for g in ("centroid", "exemplar", "balanced3"):
        if feat["cal", g].is_file() and feat["sp_temporal_test", g].is_file():
            c, t = pd.read_csv(feat["cal", g]), pd.read_csv(feat["sp_temporal_test", g])
            te[g] = t
            # the tier-conditional (Mondrian) isotonic map of cx_07, fitted on the calibration split only
            conf[g, "tier"] = MondrianIso(MONDRIAN_MIN_N).fit(c.score, c.correct4, c.pred_tier).predict(t.score,
                                                                                                      t.pred_tier)
    te_c = te["centroid"]
    conf["centroid", "gmm"] = te_c.clean_gmm.to_numpy()
    y = te_c.correct4.to_numpy().astype(float)
    tiers = [t for t in TIERS if (te_c.pred_tier == t).sum() >= MIN_N]
    if not tiers:
        return print("  skip fig5-7: no predicted-EC tier has >= 10 test proteins")
    blue = METHOD["centroid"][1]
    fig, axes = plt.subplots(1, len(tiers), figsize=(1.8 * len(tiers), 2.3), sharey=True, squeeze=False,
                             gridspec_kw={"wspace": 0.12})
    rows = []
    kinds = (("gmm", "CLEAN GMM confidence (as released)", "--", "white"),
             ("tier", "Tier-conditional isotonic, fitted on CAL", "-", blue))
    for ax, t in zip(axes[0], tiers):
        msk = (te_c.pred_tier == t).to_numpy()
        eces = {}
        for key, lab, ls, fc in kinds:
            cf = np.asarray(conf["centroid", key])[msk]
            bins, _ = reliability(cf, y[msk], n_bins=5, strategy="quantile")
            _, eces[key] = reliability(cf, y[msk])  # the reported ECE: 10 equal-width bins, as cx_07
            cis = _acc_ci(bins)
            xs, ys = [b["mean_conf"] for b in bins], [b["acc"] for b in bins]
            ax.plot(xs, ys, color=blue, linestyle=ls, marker="o" if key == "gmm" else "s", markersize=4,
                    markerfacecolor=fc, markeredgecolor=blue, label=lab)
            err = [[a - lo for a, (lo, hi) in zip(ys, cis)], [hi - a for a, (lo, hi) in zip(ys, cis)]]
            ax.errorbar(xs, ys, yerr=err, fmt="none", ecolor=blue, elinewidth=0.6, capsize=0, alpha=0.8)
            rows += [{"tier": t, "confidence": lab, "n_tier": int(msk.sum()), "ECE_reported": eces[key], **b,
                      "acc_ci_lo": lo, "acc_ci_hi": hi} for b, (lo, hi) in zip(bins, cis)]
        ax.text(0.03, 0.97, f"ECE released {eces['gmm']:.2f}\nECE tier-cond. {eces['tier']:.2f}",
                transform=ax.transAxes, va="top", ha="left", fontsize=6, color=INK2)
        ax.plot([0, 1], [0, 1], color=MUTED, linewidth=0.8)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.set_xticks([0, 0.5, 1], ["0", "0.5", "1"])
        ax.set_title(f"predicted EC {t.replace('_', ' ')}\n(n={int(msk.sum())})", fontsize=8)
        ax.set_xlabel("Mean confidence")
        ax.grid(axis="x")
    axes[0][0].set_ylabel("Observed accuracy (exact EC)")
    axes[0][0].legend(loc="lower left", bbox_to_anchor=(0, 1.25), ncol=2, fontsize=7)
    note(fig, "CLEAN centroid top-1 predictions on SP-temporal test (the same proteins for both confidences); "
         "5 equal-count bins; whiskers: exact 95% CI.")
    save(fig, "fig5_reliability_by_tier", pd.DataFrame(rows))

    fig, ax = plt.subplots(figsize=(3.6, 2.7))
    rows = []
    curves = [("centroid", "gmm", "CLEAN centroid, GMM confidence as released", "--"),
              ("centroid", "tier", "CLEAN centroid, tier-conditional isotonic", "-"),
              ("exemplar", "tier", "Exemplar (1-NN), tier-conditional isotonic", "-"),
              ("balanced3", "tier", "Size-balanced 3-NN, tier-conditional isotonic", "-")]
    for g, key, lab, ls in curves:
        if (g, key) not in conf:
            continue
        cov, risk, aurc = risk_coverage(np.asarray(conf[g, key]), te[g].correct4.to_numpy())
        ax.plot(cov, risk, color=METHOD[g][1], linestyle=ls, linewidth=1.3, label=f"{lab} (AURC {aurc:.3f})")
        rows += [{"geometry": g, "confidence": key, "AURC": aurc, "coverage": c, "selective_risk": r}
                 for c, r in zip(cov[::5], risk[::5])]
    ax.set_xlabel("Coverage (share of test proteins annotated)")
    ax.set_ylabel("Selective risk (exact EC wrong)")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.grid(axis="x")
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), fontsize=6.5)
    save(fig, "fig6_risk_coverage", pd.DataFrame(rows))

    rc = base / "E4_risk_control.csv"
    if rc.is_file():
        fig_tier_risk_control(pd.read_csv(rc))


def fig_tier_risk_control(r):
    """Realized per-tier selective risk on SP-temporal test at alpha = 0.2, one panel per geometry."""
    from cleanx.calibration import risk_ci
    r = r[(r.dataset == "sp_temporal_test") & ((r.group == "all") | r.group.str.startswith("pred_tier:"))].copy()
    if r.empty:
        return print("  skip fig7: no SP-temporal test rows in E4_risk_control.csv")
    alpha = 0.2 if np.isclose(r.alpha, 0.2).any() else sorted(set(r.alpha))[0]
    r = r[np.isclose(r.alpha, alpha)]
    if "errors" not in r:  # results written before the error count was reported
        r["errors"] = [int(round(s * a)) if a else 0 for s, a in zip(r.selective_risk.fillna(0), r.accepted)]
    lo_hi = [risk_ci(int(e), int(a)) for e, a in zip(r.errors, r.accepted)]
    r["risk_ci_lo"], r["risk_ci_hi"] = [x[0] for x in lo_hi], [x[1] for x in lo_hi]
    geoms = [g for g in ("centroid", "exemplar", "balanced3") if g in set(r.geometry)]
    groups = ["all"] + [f"pred_tier:{t}" for t in TIERS if f"pred_tier:{t}" in set(r.group)]
    fig, axes = plt.subplots(1, len(geoms), figsize=(2.4 * len(geoms), 2.6), sharey=True, squeeze=False,
                             gridspec_kw={"wspace": 0.08})
    short = {"common": "common", "rare": "rare", "very_rare": "very\nrare", "singleton": "singleton"}
    for ax, g in zip(axes[0], geoms):
        col = METHOD[g][1]
        for k, (rule, mk, fc) in enumerate((("global", "o", "white"), ("mondrian", "s", col))):
            s = r[(r.geometry == g) & (r.rule == rule)].set_index("group").reindex(groups)
            xs = np.arange(len(groups)) + (k - 0.5) * 0.36
            for x_, row in zip(xs, s.itertuples()):
                if not np.isfinite(row.accepted):  # group absent for this rule
                    continue
                ax.text(x_, 0.985, f"{row.coverage:.0%}", ha="center", va="top", fontsize=5, color=INK2, rotation=90)
                if row.accepted == 0:  # the rule abstains on the whole group: no risk to show
                    continue
                ax.plot([x_, x_], [row.risk_ci_lo, row.risk_ci_hi], color=col, linewidth=0.9)
                ax.plot([x_], [row.selective_risk], marker=mk, markersize=4.5, markerfacecolor=fc,
                        markeredgecolor=col, linestyle="none")
        ax.axhline(alpha, color=MUTED, linewidth=0.8, linestyle="--")
        ax.axvline(0.5, color=GRID, linewidth=0.8)
        ax.set_xticks(range(len(groups)), ["all"] + [short[x.split(":")[1]] for x in groups[1:]], fontsize=6)
        ax.set_xlim(-0.55, len(groups) - 0.45)
        ax.set_ylim(0, 1.0)
        ax.set_title(mlabel(g), fontsize=8)
    axes[0][len(geoms) // 2].set_xlabel("Tier of the predicted EC")
    axes[0][0].set_ylabel("Realized selective risk on test\n(95% exact CI)")
    handles = [Line2D([], [], marker="o", color=INK2, markerfacecolor="white", linestyle="none", markersize=4.5),
               Line2D([], [], marker="s", color=INK2, markerfacecolor=INK2, linestyle="none", markersize=4.5),
               Line2D([], [], color=MUTED, linewidth=0.8, linestyle="--")]
    fig.legend(handles, ["one global threshold", "one threshold per predicted-EC tier (Mondrian)",
                         f"target risk {alpha:g}"], loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3, fontsize=7)
    note(fig, "Numbers at the top: coverage (share of the group's test proteins accepted; 0% = the rule abstains on "
         f"the whole group). SGR thresholds certified on SP-temporal CAL at risk {alpha:g}, delta "
         f"{r.delta.iloc[0]:g} (split over the {len(groups) - 1} tiers for the per-tier rule).")
    cols = ["geometry", "rule", "alpha", "group", "n", "accepted", "coverage", "errors", "selective_risk", "risk_ci_lo",
            "risk_ci_hi", "theta"]
    save(fig, "fig7_tier_risk_control", r[[c for c in cols if c in r]])


def _pick(files, *tags):
    """First results folder among tags that holds every file in files (a partial or unrelated folder is skipped)."""
    for t in tags:
        if all((config.CX_RESULTS / t / f).is_file() for f in files):
            return t
    return None


def _reversals(tag, ds):
    """Rank reversals of the pre-registered H5 family (EC4 micro-F1, 100% regime, as published)."""
    p = config.CX_RESULTS / tag / "E5_pairs.csv"
    if not p.is_file():
        return None
    x = pd.read_csv(p, dtype={"regime": str})
    return x[(x.dataset == ds) & (x.regime == "100") & (x.variant == "as_published") & (x.depth == 4) &
             (x.metric == "micro_F1")]


def fig_ecbench_rescoring(tag=None):
    """E5: reported (support-weighted) vs micro F1 for EC-Bench's 16 methods, and inflation vs ECs predicted.
    tag: results folder; by default the confirmatory E5, else the disclosed pilot."""
    tag = tag or _pick(("E5_cards.csv",), "E5", "E5_pilot")
    if tag is None:
        return print("  skip fig8: run cx_09_evalcard.py (or --pilot)")
    c = pd.read_csv(config.CX_RESULTS / tag / "E5_cards.csv", dtype={"regime": str})
    ds = "test" if tag == "E5" else "price149"
    g = c[(c.source == "ecbench") & (c.dataset == ds) & (c.regime == "100") & (c.depth == 4) &
          (c.variant == "as_published")].copy()  # all methods, as in the pre-registered S4
    if g.empty:
        return print(f"  skip fig8: no EC-Bench rows for {ds} in {tag}")
    g = g.sort_values("reported_ECBench_F1")
    rev = _reversals(tag, ds)
    sig = rev[rev.significant_reversal.astype(bool)] if rev is not None else pd.DataFrame()
    marked = set(sig.A) | set(sig.B) if len(sig) else set()
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.5, 3.0), gridspec_kw={"wspace": 0.55, "width_ratios": [1.15, 1]})
    y = np.arange(len(g))
    for yy, r, m in zip(y, g.reported_ECBench_F1, g.u_micro_F1):
        a1.plot([r, m], [yy, yy], color=AXIS, linewidth=1.0, zorder=1)
    a1.scatter(g.reported_ECBench_F1, y, s=22, color="#2a78d6", marker="o", zorder=3,
               edgecolors="white", linewidths=0.6, label="Reported (support-weighted)")
    a1.scatter(g.u_micro_F1, y, s=22, color="#eb6834", marker="s", zorder=3,
               edgecolors="white", linewidths=0.6, label="Every FP counted (micro)")
    a1.set_yticks(y, [f"{m} †" if m in marked else m for m in g.method], fontsize=6.5)
    a1.set_xlabel("F1 at EC level 4")
    a1.set_xlim(0, max(0.05, float(max(g.reported_ECBench_F1.max(), g.u_micro_F1.max()))) * 1.12)
    a1.grid(axis="x"); a1.grid(axis="y", visible=False)
    a1.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=1, fontsize=6.5)
    infl = g.reported_ECBench_F1 - g.u_micro_F1
    a2.axhline(0, color=AXIS, linewidth=0.8)
    a2.scatter(g.mean_pred_labels, infl, s=24, color="#2a78d6", marker="o", edgecolors="white", linewidths=0.6)
    gi = g.assign(infl=infl, abst=g.abstained / g.n)
    for _, r in gi.sort_values("infl", ascending=False).head(3).iterrows():
        a2.annotate(r.method, (r.mean_pred_labels, r.infl), xytext=(4, 2), textcoords="offset points",
                    fontsize=6, color=INK2)
    ab = gi[gi.abst > 0.5]
    if len(ab):
        rr = ab.loc[ab.mean_pred_labels.idxmax()]
        a2.annotate(f"{', '.join(ab.sort_values('mean_pred_labels').method)}\n(abstain on "
                    f"{ab.abst.min():.0%}-{ab.abst.max():.0%} of proteins)",
                    (rr.mean_pred_labels, rr.infl), xytext=(8, -4), textcoords="offset points", fontsize=6,
                    color=INK2, va="center")
    from scipy.stats import spearmanr
    rho = spearmanr(g.mean_pred_labels, infl)[0]
    a2.text(0.03, 0.75, f"Spearman rho = {rho:.2f}\n(n = {len(g)} methods)", transform=a2.transAxes,
            va="center", fontsize=6.5, color=INK2)
    a2.set_xlabel("Mean ECs predicted per protein")
    a2.set_ylabel("Reported F1 - micro F1")
    title = "EC-Bench test set" if ds == "test" else "EC-Bench Price-149 (pilot)"
    fig.suptitle(f"{title}, models trained on the 100% clustering", fontsize=8, color=INK, y=1.08)
    if rev is not None and len(rev):
        txt = f"{int(rev.reversed.sum())} of {len(rev)} method pairs change order between the two F1s; "
        shown = [f"{r.A} vs {r.B} (q = {r.q_bh:.3f})" for r in sig.sort_values("q_bh").itertuples()]
        more = f"; +{len(shown) - 3} more" if len(shown) > 3 else ""
        txt += ("† significant (BH q < 0.05, sign-consistent): " + "; ".join(shown[:3]) + more) if shown \
            else "none significantly (BH q < 0.05)"
        note(fig, txt)
    name = "fig8_ecbench_rescoring" if tag == "E5" else "figS2_ecbench_rescoring_pilot"
    save(fig, name, gi.assign(inflation=infl, significant_reversal=g.method.isin(marked))[
        ["method", "reported_ECBench_F1", "u_micro_F1", "u_macro_union_F1", "mean_pred_labels", "abstained", "n",
         "inflation", "fp", "fp_invisible", "significant_reversal"]])


def fig_cert_budget(tag=None):
    """E6: calibration proteins needed to certify selective risk, universally and per predicted tier.
    tag: results folder; by default the confirmatory E6, else the disclosed pilot."""
    tag = tag or _pick(("E6_nstar_grid.csv", "E6_tier_budget.csv"), "E6", "E6_pilot")
    if tag is None:
        return print("  skip fig9: run cx_10_certbudget.py (or --pilot)")
    grid = pd.read_csv(config.CX_RESULTS / tag / "E6_nstar_grid.csv")
    bud = pd.read_csv(config.CX_RESULTS / tag / "E6_tier_budget.csv")
    acc_p = config.CX_RESULTS / tag / "E6_accrual.csv"
    acc = pd.read_csv(acc_p) if acc_p.is_file() else pd.DataFrame(columns=["group", "per_year"])
    rate = dict(zip(acc.group, acc.per_year))
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.0, 2.8), gridspec_kw={"wspace": 0.42})
    styles = [(0.1, INK, "o"), (0.2, INK2, "s"), (0.3, MUTED, "^")]  # one grey ramp: darker = stricter target
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
    have = b.n_available.to_numpy(float)
    need = b.n_needed.to_numpy(float)
    lo = have + b.deficit_ci_lo.to_numpy(float)
    hi = have + b.deficit_ci_hi.to_numpy(float)
    fin = np.concatenate([v[np.isfinite(v)] for v in (need, have, lo, hi)] + [np.array([10.0])])
    top = float(np.nanmax(fin)) * 4
    blue = METHOD["centroid"][1]
    yrs_lo = []
    for xi, h_, n_, l_, u_, t in zip(x, have, need, lo, hi, TIERS):
        per = rate.get(f"pred_tier:{t}", np.nan)
        yl = (l_ - h_) / per if np.isfinite(l_) and l_ > h_ and per > 0 else np.nan
        yrs_lo.append(yl)
        if np.isfinite(l_):
            a2.plot([xi, xi], [max(min(h_, l_), 1), max(l_, 1)], color=AXIS, linewidth=1.0, zorder=1)
            a2.plot([xi, xi], [max(l_, 1), u_ if np.isfinite(u_) else top], color=blue, linewidth=1.4, zorder=2)
            if not np.isfinite(u_):
                a2.plot([xi], [top], marker="^", markersize=4, color=blue, zorder=3)
            if np.isfinite(yl) and yl >= 1:  # floor keeps the "at least" statement true
                a2.annotate(f"≥{np.floor(yl):.0f} y", (xi, l_), xytext=(5, 0), textcoords="offset points",
                            fontsize=6, color=INK2, va="center", ha="left")
        if not np.isfinite(n_) and np.isfinite(h_):
            a2.annotate("not certifiable\n(point estimate)", (xi, top), xytext=(0, 6), textcoords="offset points",
                        ha="center", fontsize=5.5, color=INK2)
    a2.scatter(x, np.maximum(have, 1), s=26, color=INK2, marker="s", zorder=3, edgecolors="white", linewidths=0.6,
               label="available (calibration split)")
    ok = np.isfinite(need)
    a2.scatter(x[ok], need[ok], s=26, color=blue, marker="o", zorder=4, edgecolors="white", linewidths=0.6,
               label="needed (point estimate)")
    a2.plot([], [], color=blue, linewidth=1.4, label="needed, 95% bootstrap CI\n(▲ = unbounded)")
    a2.set_yscale("log")
    a2.set_xticks(x, [TIER_LABEL[t] for t in TIERS], fontsize=6.5)
    a2.set_ylabel("Calibration proteins, CLEAN centroid\n(target risk 0.2)")
    a2.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2, fontsize=6.5)
    rows += [{"panel": "b", "tier": t, "n_available": h_, "n_needed": n_, "needed_ci_lo": l_, "needed_ci_hi": u_,
              "years_point": y_, "years_at_ci_lo": yl}
             for t, h_, n_, l_, u_, y_, yl in zip(TIERS, have, need, lo, hi, b.years_of_accrual_for_deficit, yrs_lo)]
    if rate:
        note(fig, "≥ y: years of Swiss-Prot accrual that the deficit at its lower 95% bound represents (the "
             "tier's share of the new enzymes per year).")
    name = "fig9_certification_budget" if tag == "E6" else "figS3_certification_budget_pilot"
    save(fig, name, pd.DataFrame(rows))


def main():
    global TOOL
    ap = argparse.ArgumentParser()
    ap.add_argument("--backbone", default="split100")
    args = ap.parse_args()
    style()
    TOOL = identity_tool()
    fig_blind_spot()
    fig_tiers(args.backbone)
    fig_identity(args.backbone)
    fig_primary_forest()
    fig_ablation(args.backbone)
    fig_nise(args.backbone)
    fig_reliability()
    fig_ecbench_rescoring()
    fig_cert_budget()
    # the disclosed pilots (Supplement) are redrawn with the same code whenever the confirmatory results exist
    if _pick(("E5_cards.csv",), "E5") and _pick(("E5_cards.csv",), "E5_pilot"):
        fig_ecbench_rescoring("E5_pilot")
    if _pick(("E6_nstar_grid.csv", "E6_tier_budget.csv"), "E6") and \
            _pick(("E6_nstar_grid.csv", "E6_tier_budget.csv"), "E6_pilot"):
        fig_cert_budget("E6_pilot")
    print(f"CX08 DONE -> {config.CX_FIGS}")


if __name__ == "__main__":
    main()
