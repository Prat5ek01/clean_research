"""EC evaluation card: every averaging scheme under both label-space conventions, the false positives each
one cannot see, and a numerical check of the invariance results below.

Notation. Proteins i = 1..n with true EC sets T_i and predicted sets P_i (EC strings truncated to the evaluated
depth). L_T = union of the T_i (the evaluated label set); L_U = union of the T_i and P_i. For a class c:
TP_c = #{i: c in T_i and P_i}, FP_c = #{i: c in P_i, not T_i}, FN_c = #{i: c in T_i, not P_i}, support
s_c = TP_c + FN_c, S = sum_c s_c. An out-of-label-set augmentation adds sets A_i with A_i disjoint from L_T.

Proposition 1 (invariance). Under any out-of-label-set augmentation, (TP_c, FP_c, FN_c) is unchanged for every
c in L_T. Hence every score that depends only on these counts for c in L_T is unchanged, in particular
  (i)  every averaging (micro, macro, weighted, samples) computed after projecting predictions onto L_T
       (the protocol of CLEAN's get_eval_metrics, which fits its MultiLabelBinarizer on the test labels), and
  (ii) the support-weighted and the true-label-macro average computed on L_U (e.g. EC-Bench's scoring),
       because every class outside L_T has s_c = 0 and therefore weight s_c / S = 0, and is excluded from the
       true-label macro mean.
Proof: A_i is disjoint from L_T, so for c in L_T membership of c in P_i is unchanged; for (ii) the classes that
change all have zero support.

Proposition 2 (sensitivity). If some A_i is non-empty and sum_c TP_c > 0, micro precision on L_U strictly
decreases; samples-averaged precision strictly decreases if some protein with a correct prediction receives an
added label; the union-label macro F1 cannot increase, and strictly decreases when it was positive and an added
class was not yet in L_U (the class enters the mean with F1 = 0; classes outside L_T already had F1 = 0).

Lemma 3. On any label space containing L_T, support-weighted recall equals micro recall
(sum_c (s_c/S)(TP_c/s_c) = sum_c TP_c / S). The whole difference between the reported and the micro scores is
therefore on the precision side.

Proposition 4 (over-prediction). Adding a correct label c to P_i raises TP_c and lowers FN_c, so it weakly
raises P_c, R_c and F1_c; adding a wrong label inside L_T weakly lowers F1_c; adding a wrong label outside L_T
changes nothing (Prop. 1). Top-k expansion therefore raises the support-weighted F1 whenever its gains in
recall outweigh the costs of the wrong labels that happen to fall inside L_T, and the smaller L_T is relative to
the predictor's vocabulary, the larger the share of wrong labels that are free."""
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
