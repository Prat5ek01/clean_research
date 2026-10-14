"""Paired method contrasts and independently resampled stratum differences.

Identity standardisation uses pooled bin weights and resamples each stratum/bin."""
import numpy as np

from . import metrics as M


def two_sided_p(d):
    """Two-sided bootstrap p-value of a difference: twice the smaller tail mass at zero, capped at 1."""
    d = np.asarray(d, float)
    return float(min(1.0, 2 * min(np.mean(d <= 0), np.mean(d >= 0))))


def delta(truth, pa, pb, ids, key, B, seed, depth=4, clusters=None):
    """Observed metric(A) - metric(B) on ids and its (B,) paired-bootstrap distribution."""
    bs = M.bootstrap(truth, {"a": pa, "b": pb}, ids, depth=depth, B=B, seed=seed, keys=(key,), clusters=clusters)
    obs = M.evaluate(truth, pa, ids, depth)[key] - M.evaluate(truth, pb, ids, depth)[key]
    return float(obs), bs["a"][key] - bs["b"][key]


def summarize(obs, boot):
    lo, hi = M.ci(boot)
    return {"estimate": float(obs), "ci_lo": lo, "ci_hi": hi, "p_two_sided": two_sided_p(boot),
            "se_boot": float(np.std(boot, ddof=1)), "mde80": mde80(np.std(boot, ddof=1))}


def mde80(se):
    """Minimum detectable effect at 80% power, two-sided alpha 0.05: (z_0.975 + z_0.80) * SE = 2.80 * SE."""
    return float((1.959963984540054 + 0.8416212335729143) * se)


def did(truth, pa, pb, ids0, ids1, key, B, seed, depth=4, clusters=None):
    """DiD = delta(ids1) - delta(ids0), strata resampled independently with seeds seed (S0) and seed + 1 (S1).

    Returns {'delta0', 'delta1', 'did'}, each a summarize() dict, plus the bootstrap arrays."""
    o0, b0 = delta(truth, pa, pb, ids0, key, B, seed, depth, clusters)
    o1, b1 = delta(truth, pa, pb, ids1, key, B, seed + 1, depth, clusters)
    return {"delta0": summarize(o0, b0), "delta1": summarize(o1, b1), "did": summarize(o1 - o0, b1 - b0),
            "boot0": b0, "boot1": b1}


def merge_bins(order, counts, min_n):
    """Merge ordered bins until every stratum has >= min_n proteins in every merged bin.

    order: bin labels from low to high identity; counts: {stratum: {bin: n}}. Rule (fixed): scan from the highest
    bin downwards and merge a short bin into the next lower one; then, while the lowest merged bin is short, merge it
    into the next higher one. Returns [[bins], ...] from low to high, or None if even one bin is short."""
    groups = [[b] for b in order]

    def short(g):
        return any(sum(c.get(b, 0) for b in g) < min_n for c in counts.values())

    i = len(groups) - 1
    while i >= 1:
        if short(groups[i]):
            groups[i - 1] = groups[i - 1] + groups[i]
            del groups[i]
        i -= 1
    while len(groups) > 1 and short(groups[0]):
        groups[1] = groups[0] + groups[1]
        del groups[0]
    if short(groups[0]):
        return None
    return groups


def standardized(truth, pa, pb, strata, bin_of, order, key, B, seed, depth=4, min_n=10):
    """Identity-standardised deltas of every stratum and their bootstrap distributions.

    strata: {name: [ids]} (disjoint); bin_of: id -> bin label; order: bins low -> high. Weights w are the pooled
    share of each merged bin over all strata. Cell (stratum k, merged bin j) is resampled with seed + 1000*k + j.
    Returns None if the bins cannot be merged to min_n per cell, else {'groups', 'weights', 'strata': {name:
    (observed, (B,) array)}}."""
    counts = {s: {} for s in strata}
    for s, ids in strata.items():
        for i in ids:
            counts[s][bin_of[i]] = counts[s].get(bin_of[i], 0) + 1
    groups = merge_bins(order, counts, min_n)
    if groups is None:
        return None
    total = sum(len(ids) for ids in strata.values())
    w = np.array([sum(counts[s].get(b, 0) for s in strata for b in g) for g in groups], float) / total
    out = {}
    for k, (s, ids) in enumerate(strata.items()):
        obs, boot = 0.0, np.zeros(B)
        for j, g in enumerate(groups):
            cell = [i for i in ids if bin_of[i] in g]
            o, b = delta(truth, pa, pb, cell, key, B, seed + 1000 * k + j, depth)
            obs += w[j] * o
            boot += w[j] * b
        out[s] = (float(obs), boot)
    return {"groups": ["+".join(g) for g in groups], "weights": w.tolist(), "strata": out}
