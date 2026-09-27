import math
import pickle

import numpy as np
from scipy.stats import beta

from . import config


def load_gmms(path=config.GMM_PKL):
    with open(path, "rb") as fh:
        return pickle.load(fh)


def gmm_confidence(dist, gmms):
    """CLEAN's released confidence (CLEAN.evaluate.infer_confidence_gmm) for each distance."""
    from CLEAN.evaluate import infer_confidence_gmm
    return np.array([infer_confidence_gmm(float(d), gmms) for d in np.ravel(dist)])


def gmm_confidence_fast(dist, gmms):
    """Same computation as gmm_confidence, one predict_proba call per mixture component model."""
    x = np.asarray(dist, float).reshape(-1, 1)
    out = []
    for g in gmms:
        a, b = g.means_
        out.append(g.predict_proba(x)[:, 0 if a[0] < b[0] else 1])
    return np.mean(out, axis=0)


def reliability(conf, correct, n_bins=10, strategy="uniform"):
    conf = np.asarray(conf, float)
    y = np.asarray(correct, float)
    if strategy == "uniform":
        edges = np.linspace(0, 1, n_bins + 1)
    else:
        edges = np.unique(np.quantile(conf, np.linspace(0, 1, n_bins + 1)))
    idx = np.clip(np.searchsorted(edges, conf, side="right") - 1, 0, len(edges) - 2)
    rows, ece = [], 0.0
    for b in range(len(edges) - 1):
        m = idx == b
        if not m.any():
            continue
        rows.append({"lo": float(edges[b]), "hi": float(edges[b + 1]), "n": int(m.sum()),
                     "mean_conf": float(conf[m].mean()), "acc": float(y[m].mean())})
        ece += m.mean() * abs(conf[m].mean() - y[m].mean())
    return rows, float(ece)


def brier(conf, correct):
    return float(np.mean((np.asarray(conf, float) - np.asarray(correct, float)) ** 2))


def risk_coverage(conf, correct):
    """Selective risk at every coverage level (most confident first) and the area under it (AURC)."""
    conf = np.asarray(conf, float)
    err = 1.0 - np.asarray(correct, float)
    o = np.argsort(-conf, kind="stable")
    cum = np.cumsum(err[o])
    n = np.arange(1, len(conf) + 1)
    risk = cum / n
    return n / len(conf), risk, float(risk.mean()) if len(conf) else float("nan")


class IsoCalibrator:
    """Isotonic map from a score (larger = more confident) to P(correct)."""

    def fit(self, score, correct):
        from sklearn.isotonic import IsotonicRegression
        self.iso = IsotonicRegression(increasing=True, out_of_bounds="clip", y_min=0.0, y_max=1.0)
        self.iso.fit(np.asarray(score, float), np.asarray(correct, float))
        return self

    def predict(self, score):
        return self.iso.predict(np.asarray(score, float))


class MondrianIso:
    """Separate isotonic calibration per group (e.g. tier of the predicted EC); small groups
    fall back to the pooled calibrator."""

    def __init__(self, min_n=50):
        self.min_n = min_n

    def fit(self, score, correct, groups):
        score, correct, groups = np.asarray(score, float), np.asarray(correct, float), np.asarray(groups)
        self.pooled = IsoCalibrator().fit(score, correct)
        self.per = {}
        for g in np.unique(groups):
            m = groups == g
            if m.sum() >= self.min_n:
                self.per[g] = IsoCalibrator().fit(score[m], correct[m])
        return self

    def predict(self, score, groups):
        score, groups = np.asarray(score, float), np.asarray(groups)
        out = self.pooled.predict(score)
        for g, cal in self.per.items():
            m = groups == g
            if m.any():
                out[m] = cal.predict(score[m])
        return out


def binom_ucb(k, n, delta):
    """One-sided (1-delta) Clopper-Pearson upper bound on an error rate after k errors in n."""
    if n == 0:
        return 1.0
    if k >= n:
        return 1.0
    return float(beta.ppf(1 - delta, k + 1, n - k))


def sgr_threshold(score, correct, alpha, delta):
    """Selection with Guaranteed Risk (Geifman & El-Yaniv, NeurIPS 2017, Algorithm 1).

    Accept predictions with score >= theta. Returns (theta, bound, coverage_on_cal); theta=+inf
    means no threshold could be certified (abstain on everything). With probability >= 1-delta
    the selective risk of the returned rule is <= bound (< alpha when certified)."""
    s = np.sort(np.asarray(score, float))
    y = np.asarray(correct, bool)
    sc = np.asarray(score, float)
    m = len(s)
    if m == 0:
        return math.inf, 1.0, 0.0
    k = max(1, math.ceil(math.log2(m)))
    zmin, zmax = 0, m - 1
    best = (math.inf, 1.0, 0.0)
    for _ in range(k):
        z = math.ceil((zmin + zmax) / 2)
        theta = s[z]
        acc = sc >= theta
        n = int(acc.sum())
        errs = int((~y[acc]).sum())
        b = binom_ucb(errs, n, delta / k)
        if b < alpha:
            if theta < best[0]:
                best = (float(theta), b, n / m)
            zmax = z
        else:
            zmin = z
        if zmax - zmin <= 0:
            break
    return best


def apply_threshold(score, theta):
    return np.asarray(score, float) >= theta


def selective_summary(accept, correct):
    accept = np.asarray(accept, bool)
    correct = np.asarray(correct, bool)
    n = len(accept)
    na = int(accept.sum())
    return {"n": n, "accepted": na, "coverage": na / n if n else float("nan"),
            "selective_risk": float((~correct[accept]).mean()) if na else float("nan")}
