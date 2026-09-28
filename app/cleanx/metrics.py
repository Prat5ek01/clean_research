"""Multi-label EC metrics with an explicit label space.

'weighted' reproduces the number CLEAN reports (sklearn support-weighted average over classes).
Classes that are predicted but have no true instance get weight 0 there, so false positives on EC
numbers absent from the evaluated set's labels never reduce it; 'macro_true' (mean over classes
with a true instance) shares that blind spot. Both can be raised by over-predicting. 'micro',
'macro_union' (mean over true + predicted classes) and 'samples' count every prediction."""
import numpy as np
from scipy import sparse

from .labels import truncate

SCHEMES = ("micro", "weighted", "macro_true", "macro_union", "samples")


def _sets(truth, preds, ids, depth):
    T = [{truncate(e, depth) for e in truth[i]} for i in ids]
    P = [{truncate(e, depth) for e in (preds.get(i) or [])} for i in ids]
    return T, P


def _matrices(T, P):
    classes = sorted(set().union(*T, *P)) if T else []
    cidx = {c: k for k, c in enumerate(classes)}
    rows = {"tp": ([], []), "fp": ([], []), "fn": ([], [])}
    for i, (t, p) in enumerate(zip(T, P)):
        for c in t & p:
            rows["tp"][0].append(i); rows["tp"][1].append(cidx[c])
        for c in p - t:
            rows["fp"][0].append(i); rows["fp"][1].append(cidx[c])
        for c in t - p:
            rows["fn"][0].append(i); rows["fn"][1].append(cidx[c])
    shape = (len(T), max(len(classes), 1))
    mats = {k: sparse.csr_matrix((np.ones(len(r)), (r, c)), shape=shape) for k, (r, c) in rows.items()}
    return classes, mats


def _div(a, b):
    return np.divide(a, b, out=np.zeros_like(a, dtype=float), where=b > 0)


def _class_scores(tp, fp, fn):
    """tp/fp/fn: (..., C) per-class counts -> dict of scheme -> (P, R, F1) arrays over leading dims."""
    tp, fp, fn = (np.asarray(x, dtype=float) for x in (tp, fp, fn))
    Pc, Rc = _div(tp, tp + fp), _div(tp, tp + fn)
    Fc = _div(2 * Pc * Rc, Pc + Rc)
    sup = tp + fn
    seen = (tp + fp + fn) > 0
    true_cls = sup > 0
    out = {}
    TP, FP, FN = tp.sum(-1), fp.sum(-1), fn.sum(-1)
    mP, mR = _div(TP, TP + FP), _div(TP, TP + FN)
    out["micro"] = (mP, mR, _div(2 * mP * mR, mP + mR))
    w = sup.sum(-1)
    out["weighted"] = tuple(_div((sup * x).sum(-1), w) for x in (Pc, Rc, Fc))
    nt = true_cls.sum(-1)
    out["macro_true"] = tuple(_div((x * true_cls).sum(-1), nt) for x in (Pc, Rc, Fc))
    ns = seen.sum(-1)
    out["macro_union"] = tuple(_div((x * seen).sum(-1), ns) for x in (Pc, Rc, Fc))
    return out


def evaluate(truth, preds, ids, depth=4):
    """Point estimates for one method on one id subset. preds: id -> list of EC strings."""
    T, P = _sets(truth, preds, ids, depth)
    classes, M = _matrices(T, P)
    tp, fp, fn = (np.asarray(M[k].sum(0)).ravel() for k in ("tp", "fp", "fn"))
    sc = _class_scores(tp, fp, fn)
    ti = np.array([len(t & p) for t, p in zip(T, P)], float)
    np_ = np.array([len(p) for p in P], float)
    nt = np.array([len(t) for t in T], float)
    sP, sR = _div(ti, np_), _div(ti, nt)
    sc["samples"] = (sP.mean() if len(ids) else 0.0, sR.mean() if len(ids) else 0.0,
                     _div(2 * ti, np_ + nt).mean() if len(ids) else 0.0)
    true_labels = set().union(*T) if T else set()
    fp_oov = sum(1 for t, p in zip(T, P) for c in p - t if c not in true_labels)
    res = {"n": len(ids), "depth": depth, "tp": int(tp.sum()), "fp": int(fp.sum()), "fn": int(fn.sum()),
           "fp_outside_true_labelset": int(fp_oov),
           "n_pred_labels": int(np_.sum()), "n_true_labels": int(nt.sum()),
           "abstained": int((np_ == 0).sum()),
           "subset_acc": float(np.mean([t == p for t, p in zip(T, P)])) if len(ids) else 0.0}
    for s in SCHEMES:
        p_, r_, f_ = sc[s]
        res[f"{s}_P"], res[f"{s}_R"], res[f"{s}_F1"] = float(p_), float(r_), float(f_)
    return res


def bootstrap(truth, preds_by_method, ids, depth=4, B=2000, seed=0,
              keys=("micro_F1", "macro_union_F1", "macro_true_F1", "weighted_F1", "samples_F1"), clusters=None):
    """Paired nonparametric bootstrap over proteins (same resamples for every method).

    clusters: optional id -> cluster label. Clusters are resampled with replacement and every protein takes its
    cluster's multiplicity (cluster bootstrap for correlated proteins, e.g. sharing an EC). With every protein
    in its own cluster, in id order, the resamples are identical to the protein bootstrap.

    Returns {method: {key: (B,) array}}."""
    rng = np.random.default_rng(seed)
    N = len(ids)
    if clusters is None:
        W = rng.multinomial(N, np.full(N, 1.0 / N), size=B).astype(float)
    else:
        labs = list(dict.fromkeys(clusters[i] for i in ids))
        pos = {c: k for k, c in enumerate(labs)}
        Wc = rng.multinomial(len(labs), np.full(len(labs), 1.0 / len(labs)), size=B).astype(float)
        W = Wc[:, [pos[clusters[i]] for i in ids]]
    out = {}
    for m, preds in preds_by_method.items():
        T, P = _sets(truth, preds, ids, depth)
        _, M = _matrices(T, P)
        tp, fp, fn = (np.asarray(W @ M[k]) for k in ("tp", "fp", "fn"))
        sc = _class_scores(tp, fp, fn)
        ti = np.array([len(t & p) for t, p in zip(T, P)], float)
        np_ = np.array([len(p) for p in P], float)
        nt = np.array([len(t) for t in T], float)
        # weighted mean over the resample: the weights sum to N for the protein bootstrap but not for clusters
        samples_f1 = (W @ _div(2 * ti, np_ + nt)) / W.sum(1)
        res = {}
        for key in keys:
            scheme, metric = key.rsplit("_", 1)
            if scheme == "samples":
                res[key] = samples_f1
            else:
                res[key] = sc[scheme][{"P": 0, "R": 1, "F1": 2}[metric]]
        out[m] = res
    return out


def ci(x, level=0.95):
    a = (1 - level) / 2
    return float(np.quantile(x, a)), float(np.quantile(x, 1 - a))


def paired_delta(bs, m_new, m_ref, key):
    """(bootstrap mean of the difference, 95% percentile CI, two-sided bootstrap p). Report the observed
    full-sample difference (observed_delta) as the effect; for macro-averaged keys the bootstrap mean is biased
    because resampling changes which classes are present."""
    d = bs[m_new][key] - bs[m_ref][key]
    lo, hi = ci(d)
    p = 2 * min(np.mean(d <= 0), np.mean(d >= 0))
    return float(d.mean()), lo, hi, float(min(p, 1.0))


def observed_delta(truth, preds_new, preds_ref, ids, key, depth=4):
    """Full-sample point estimate of metric(new) - metric(ref) on the same ids."""
    return float(evaluate(truth, preds_new, ids, depth)[key] - evaluate(truth, preds_ref, ids, depth)[key])


def holm(pvals):
    """Holm step-down adjusted p-values (family-wise error control), in input order."""
    p = np.asarray(pvals, float)
    m = len(p)
    if m == 0:
        return p
    o = np.argsort(p, kind="stable")
    adj = np.minimum(1.0, np.maximum.accumulate((m - np.arange(m)) * p[o]))
    out = np.empty(m)
    out[o] = adj
    return out


def bh(pvals):
    """Benjamini-Hochberg adjusted p-values (false-discovery-rate control), in input order."""
    p = np.asarray(pvals, float)
    m = len(p)
    if m == 0:
        return p
    o = np.argsort(p, kind="stable")
    adj = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]
    out = np.empty(m)
    out[o] = np.minimum(adj, 1.0)
    return out


def group_ids(ids, groups):
    """groups: id -> label. Returns label -> [ids] in first-seen order."""
    out = {}
    for i in ids:
        out.setdefault(groups[i], []).append(i)
    return out
