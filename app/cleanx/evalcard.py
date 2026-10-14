"""Report EC metrics, omitted false positives and label-space invariance checks.

See docs/code_reference.md#evalcard."""
import numpy as np

from . import labels as L
from .metrics import SCHEMES, evaluate

PHANTOM = "99.99.99.{}"  # never a valid EC: guaranteed outside every label set


def label_space(truth, ids, depth=4):
    return {L.truncate(e, depth) for i in ids for e in truth[i]}


def project(preds, space, depth=4):
    """Drop predicted labels outside `space` (CLEAN's protocol: labels unknown to the fitted binarizer vanish)."""
    return {q: [e for e in (p or []) if L.truncate(e, depth) in space] for q, p in preds.items()}


def card(truth, preds, ids, depth=4):
    """All schemes on the union label space ('u_') and after projection onto the test label set ('t_').

    reported_CLEAN_* = support-weighted on the projected predictions (CLEAN, CLEAN-Contact, ProtDETR,
    EnzymeHunter); reported_ECBench_* = support-weighted on the union label space. By Prop. 1 both equal
    u_weighted_*; the card asserts it."""
    space = label_space(truth, ids, depth)
    u = evaluate(truth, preds, ids, depth)
    t = evaluate(truth, project(preds, space, depth), ids, depth)
    for k in ("P", "R", "F1"):
        for s in ("weighted", "macro_true"):
            if abs(u[f"{s}_{k}"] - t[f"{s}_{k}"]) > 1e-12:
                raise AssertionError(f"Proposition 1 violated for {s}_{k}: {u[f'{s}_{k}']} vs {t[f'{s}_{k}']}")
    n_pred = u["n_pred_labels"]
    out = {"n": u["n"], "depth": depth, "n_true_classes": len(space), "n_pred_labels": n_pred,
           "mean_pred_labels": n_pred / max(u["n"], 1), "abstained": u["abstained"], "tp": u["tp"], "fp": u["fp"],
           "fn": u["fn"], "fp_invisible": u["fp_outside_true_labelset"],
           "invisible_share": u["fp_outside_true_labelset"] / u["fp"] if u["fp"] else float("nan"),
           "subset_acc": u["subset_acc"]}
    for s in SCHEMES:
        for k in ("P", "R", "F1"):
            out[f"u_{s}_{k}"] = u[f"{s}_{k}"]
            out[f"t_{s}_{k}"] = t[f"{s}_{k}"]
    for k in ("P", "R", "F1"):
        out[f"reported_CLEAN_{k}"] = t[f"weighted_{k}"]
        out[f"reported_ECBench_{k}"] = u[f"weighted_{k}"]
    out["inflation_F1"] = out["reported_ECBench_F1"] - out["u_micro_F1"]
    return out


def augment_out_of_space(preds, ids, n_extra=3, seed=0):
    """Add 0..n_extra phantom labels (outside every label set) to each protein's predictions.

    Returns (augmented predictions, id -> number of phantom labels added)."""
    rng = np.random.default_rng(seed)
    out = {q: list(p or []) for q, p in preds.items()}
    added = {}
    for i in ids:
        out.setdefault(i, [])
        k = int(rng.integers(0, n_extra + 1))
        out[i] += [PHANTOM.format(int(x)) for x in rng.choice(1000, size=k, replace=False)]
        added[i] = k
    return out, added


def invariance_check(truth, preds, ids, depth=4, n_extra=3, seed=0):
    """Numerical check of Props. 1-2 on real predictions. Returns a dict; 'ok' is False on any violation."""
    a = card(truth, preds, ids, depth)
    aug, per = augment_out_of_space(preds, ids, n_extra, seed)
    b = card(truth, aug, ids, depth)
    added = b["n_pred_labels"] - a["n_pred_labels"]
    hit_augmented = any(per[i] > 0 and {L.truncate(e, depth) for e in truth[i]}
                        & {L.truncate(e, depth) for e in (preds.get(i) or [])} for i in ids)
    inv = ["reported_CLEAN_P", "reported_CLEAN_R", "reported_CLEAN_F1", "reported_ECBench_P",
           "reported_ECBench_R", "reported_ECBench_F1", "u_macro_true_F1"] + [f"t_{s}_F1" for s in SCHEMES]
    max_change = max(abs(a[k] - b[k]) for k in inv)
    res = {"phantom_labels_added": int(added), "max_change_invariant_scores": max_change,
           "micro_P_before": a["u_micro_P"], "micro_P_after": b["u_micro_P"],
           "samples_P_before": a["u_samples_P"], "samples_P_after": b["u_samples_P"],
           "macro_union_F1_before": a["u_macro_union_F1"], "macro_union_F1_after": b["u_macro_union_F1"]}
    ok = max_change < 1e-12 and abs(a["u_micro_R"] - a["reported_ECBench_R"]) < 1e-12  # Prop. 1, Lemma 3
    ok &= b["u_macro_union_F1"] <= a["u_macro_union_F1"] + 1e-15
    if added > 0 and a["tp"] > 0:
        ok &= b["u_micro_P"] < a["u_micro_P"]
    if hit_augmented:
        ok &= b["u_samples_P"] < a["u_samples_P"]
    res["ok"] = bool(ok)
    return res


def gaming_curve(truth, ranked, ids, ks=(1, 2, 3, 5, 10), depth=4):
    """ranked: id -> ECs best-first. One card per k (the first k ECs predicted for every protein)."""
    rows = []
    for k in ks:
        rows.append({"k": k, **card(truth, {q: list(ranked.get(q, []))[:k] for q in ids}, ids, depth)})
    return rows
