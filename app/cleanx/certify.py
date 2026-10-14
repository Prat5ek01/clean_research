"""Compute fixed-operating-point binomial power and calibration budgets.

Budgets use the stable power crossing within a finite search window.

See docs/code_reference.md#certify."""
import math

import numpy as np
from scipy.stats import binom

from .calibration import sgr_threshold


def k_max(n, alpha, delta):
    """Largest certifiable error count for n accepted proteins (vectorised over n); -1 = none."""
    n = np.atleast_1d(np.asarray(n, dtype=np.int64))
    k = binom.ppf(delta, n, alpha).astype(np.int64) - 1
    k = np.clip(k, -1, n - 1)
    # guard floating-point edge cases of ppf: enforce cdf(k) < delta <= cdf(k + 1)
    up = (k + 1 <= n - 1) & (binom.cdf(k + 1, n, alpha) < delta)
    k[up] += 1
    down = (k >= 0) & (binom.cdf(k, n, alpha) >= delta)
    k[down] -= 1
    return k


def power(n, r, alpha, delta):
    """Probability that n accepted proteins with true selective risk r certify alpha (vectorised over n)."""
    n = np.atleast_1d(np.asarray(n, dtype=np.int64))
    k = k_max(n, alpha, delta)
    out = np.where(k >= 0, binom.cdf(np.maximum(k, 0), n, r), 0.0)
    return out


_NSTAR = {}


def n_star(r, alpha, delta, target=0.8, n_max=200000):
    """(n_first, n_stable) accepted proteins needed; (inf, inf) if the target is not reached by n_max.

    The search window grows x4 until the power over its last fifth stays >= target + 0.02 (past the saw-tooth),
    so the cost is proportional to the answer; results are cached."""
    key = (round(float(r), 12), float(alpha), float(delta), float(target), int(n_max))
    if key in _NSTAR:
        return _NSTAR[key]
    res = (math.inf, math.inf)
    if r < alpha:
        n0 = max(1, int(math.floor(math.log(delta) / math.log(1 - alpha))) + 1)  # below n0 nothing certifies
        N = max(1000, 4 * n0)
        while True:
            N = min(N, n_max)
            n = np.arange(n0, N + 1)
            pw = power(n, r, alpha, delta)
            ok = pw >= target
            settled = len(pw) > 5 and (pw[int(0.8 * len(pw)):] >= target + 0.02).all()
            if ok.any() and (settled or N >= n_max):
                first = int(n[np.argmax(ok)])
                bad = np.nonzero(~ok)[0]
                if not len(bad):
                    stable = first
                elif bad[-1] + 1 < len(n):
                    stable = int(n[bad[-1] + 1])
                else:
                    stable = math.inf
                res = (first, stable)
                break
            if N >= n_max:
                break
            N *= 4
    _NSTAR[key] = res
    return res


def per_test_delta(delta, groups, m):
    return delta / (groups * max(1, math.ceil(math.log2(max(m, 2)))))


def group_budget(r, coverage, alpha, delta, groups=4, target=0.8, n_max=200000, max_iter=60):
    """Calibration proteins m a group needs so that its accepted share (coverage c, true selective risk r)
    certifies alpha with probability >= target under Mondrian SGR. Returns (m, n_accepted, per_test_delta).

    Iterates m <- ceil(n_stable(delta'(m)) / c) from m = 2. delta'(m) never increases with m, so the sequence is
    non-decreasing and stops at the smallest self-consistent budget (m == ceil(n_accepted / c))."""
    if r >= alpha or coverage <= 0:
        return math.inf, math.inf, float("nan")
    m = 2
    for _ in range(max_iter):
        dp = per_test_delta(delta, groups, m)
        _, n_acc = n_star(r, alpha, dp, target, n_max)
        if not math.isfinite(n_acc):
            return math.inf, math.inf, dp
        m_new = max(m, int(math.ceil(n_acc / coverage)))
        if m_new == m:
            return m, n_acc, dp
        m = m_new
    raise RuntimeError(f"group_budget did not converge (r={r}, c={coverage}, alpha={alpha})")


def percentile_ci(x, level=0.95):
    """Percentile interval of bootstrap budgets that may be +inf (resamples that cannot certify). numpy's linear
    interpolation between two infinite values returns nan; that quantile is +inf."""
    x = np.asarray(x, float)
    if not len(x):
        return float("nan"), float("nan")
    a = round((1 - level) / 2, 12)  # exactly 0.025 for 95%: (1 - 0.95) / 2 is 0.025000000000000022 in floating point
    out = []
    for q in (a, round(1 - a, 12)):
        v = float(np.quantile(x, q))
        out.append(math.inf if math.isnan(v) else v)
    return out[0], out[1]


def accepted_counts(score, correct, coverages):
    """(k, errors) among the top ceil(c * n) proteins by score (larger = more confident), per coverage c."""
    s = np.asarray(score, float)
    err = 1.0 - np.asarray(correct, float)
    o = np.argsort(-s, kind="stable")
    cum = np.concatenate([[0.0], np.cumsum(err[o])])
    ks = np.array([int(math.ceil(c * len(s))) for c in coverages])
    return ks, cum[ks]


def best_budget(score, correct, alpha, delta, groups=4, target=0.8, coverages=None, n_max=200000,
                min_accepted=10):
    """Smallest group budget over the coverage grid (the most favourable operating point).

    The selective risk at coverage c is the Jeffreys posterior mean (errors + 0.5) / (k + 1) over the k accepted
    calibration proteins, and coverages with k < min_accepted are skipped (a raw 0/3 is no estimate). Because the
    coverage is chosen with hindsight (an oracle threshold), the budget is a lower bound on what SGR's own
    threshold search needs; cx_10 reports SGR's simulated certification rate at the budget.
    status: 'estimated', 'risk_above_alpha' (no coverage reaches a risk below alpha), or
    'too_few_calibration_proteins' (no coverage has min_accepted proteins)."""
    cov = np.round(np.arange(0.1, 0.9001, 0.1), 2) if coverages is None else np.asarray(coverages, float)
    ks, errs = accepted_counts(score, correct, cov)
    best = (math.inf, float("nan"), float("nan"), math.inf)
    usable = [(c, k, e) for c, k, e in zip(cov, ks, errs) if k >= min_accepted]
    for c, k, e in usable:
        r = (e + 0.5) / (k + 1)
        if r >= alpha:
            continue
        m, n_acc, _ = group_budget(float(r), float(c), alpha, delta, groups, target, n_max)
        if m < best[0]:
            best = (m, float(c), float(r), n_acc)
    status = ("too_few_calibration_proteins" if not usable else
              ("estimated" if math.isfinite(best[0]) else "risk_above_alpha"))
    return {"n_needed": best[0], "best_coverage": best[1], "risk_at_best": best[2], "n_accepted_needed": best[3],
            "status": status}


def simulate_sgr(m, coverage, r, alpha, delta, groups=4, reps=500, r_rest=0.9, seed=0):
    """Empirical certification rate of Mondrian SGR on synthetic groups of m proteins whose top `coverage`
    fraction (by score) has error rate r and the rest r_rest."""
    rng = np.random.default_rng(seed)
    hits = 0
    for _ in range(reps):
        s = rng.uniform(size=m)
        top = s >= np.quantile(s, 1 - coverage)
        err = np.where(top, rng.uniform(size=m) < r, rng.uniform(size=m) < r_rest)
        theta, _, _ = sgr_threshold(s, ~err, alpha, delta / groups)
        hits += math.isfinite(theta)
    return hits / reps
