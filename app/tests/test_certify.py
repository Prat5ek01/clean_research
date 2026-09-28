"""Known-answer tests for cleanx.certify (exact binomial certification budgets)."""
import math

import numpy as np

from cleanx import certify as C
from cleanx.calibration import binom_ucb


def test_kmax_matches_clopper_pearson_bound():
    for n in (5, 10, 17, 50, 137, 400):
        for a in (0.1, 0.2, 0.3):
            for d in (0.1, 0.025, 0.003):
                km = int(C.k_max(n, a, d)[0])
                for k in range(0, min(n - 1, km + 6) + 1):
                    assert (binom_ucb(k, n, d) < a) == (k <= km), (n, a, d, k, km)


def test_nstar_at_zero_risk_has_closed_form():
    for a, d in ((0.2, 0.025), (0.1, 0.1), (0.3, 0.00357), (0.05, 0.01)):
        expect = math.floor(math.log(d) / math.log(1 - a)) + 1
        first, stable = C.n_star(0.0, a, d)
        assert first == expect and stable == expect, (a, d, first, stable, expect)


def test_nstar_infinite_at_or_above_alpha():
    assert C.n_star(0.2, 0.2, 0.025) == (math.inf, math.inf)
    assert C.n_star(0.35, 0.3, 0.1) == (math.inf, math.inf)


def test_nstar_monotone_in_risk():
    prev = 0
    for r in np.arange(0.0, 0.19, 0.02):
        s = C.n_star(float(r), 0.2, 0.01)[1]
        assert s >= prev, (r, s, prev)
        prev = s


def test_power_matches_monte_carlo():
    rng = np.random.default_rng(0)
    n, r, a, d = 300, 0.12, 0.2, 0.025
    pw = float(C.power(n, r, a, d)[0])
    km = int(C.k_max(n, a, d)[0])
    mc = float((rng.binomial(n, r, size=40000) <= km).mean())
    assert abs(pw - mc) < 0.01, (pw, mc)


def test_group_budget_is_a_fixed_point():
    m, n_acc, dp = C.group_budget(0.05, 0.5, 0.2, 0.1, groups=4)
    assert m == math.ceil(n_acc / 0.5)
    assert abs(dp - 0.1 / (4 * math.ceil(math.log2(m)))) < 1e-15
    assert C.n_star(0.05, 0.2, dp)[1] == n_acc


def test_best_budget_status_and_shrinkage():
    few = C.best_budget(np.linspace(1, 0, 8), np.ones(8, bool), 0.2, 0.1)
    assert few["status"] == "too_few_calibration_proteins" and not math.isfinite(few["n_needed"])
    bad = C.best_budget(np.linspace(1, 0, 200), np.zeros(200, bool), 0.2, 0.1)
    assert bad["status"] == "risk_above_alpha"
    good = C.best_budget(np.linspace(1, 0, 200), np.arange(200) < 150, 0.2, 0.1)
    assert good["status"] == "estimated" and good["risk_at_best"] > 0  # shrinkage: 0 errors never gives r = 0
    k = math.ceil(good["best_coverage"] * 200)
    assert abs(good["risk_at_best"] - 0.5 / (k + 1)) < 1e-12


def test_sgr_simulation_extremes():
    assert C.simulate_sgr(400, 0.5, 0.0, 0.2, 0.1, groups=1, reps=60, seed=1) >= 0.95
    assert C.simulate_sgr(400, 0.5, 0.35, 0.2, 0.1, groups=1, reps=60, seed=1) <= 0.1
