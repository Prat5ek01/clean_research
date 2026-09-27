import hashlib
import json

import numpy as np

from . import config


def _sqdist(Q, X, X_sq):
    d2 = np.einsum("ij,ij->i", Q, Q)[:, None] + X_sq[None, :] - 2.0 * (Q @ X.T)
    np.maximum(d2, 0.0, out=d2)
    return d2


def centroid_distances(Q, index):
    return np.sqrt(_sqdist(Q, index.centroids, index.C_sq))


def exemplar_distances(Q, index, chunk=32):
    """Per-EC distance to the nearest training member (1-NN restricted to each EC)."""
    out = np.empty((len(Q), len(index.ecs)))
    for s in range(0, len(Q), chunk):
        d2 = _sqdist(Q[s:s + chunk], index.E, index.E_sq)
        out[s:s + chunk] = np.sqrt(np.minimum.reduceat(d2, index.starts, axis=1))
    return out


def kmean_distances(Q, index, k=config.KMEAN_K, chunk=32):
    """Per-EC mean distance to the k nearest members (all members when the EC has <= k)."""
    out = np.empty((len(Q), len(index.ecs)))
    small = index.sizes <= k
    big = np.nonzero(~small)[0]
    for s in range(0, len(Q), chunk):
        d = np.sqrt(_sqdist(Q[s:s + chunk], index.E, index.E_sq))
        sums = np.add.reduceat(d, index.starts, axis=1)
        blk = out[s:s + chunk]
        blk[:, small] = sums[:, small] / index.sizes[small][None, :]
        for j in big:
            a = index.starts[j]
            n = index.sizes[j]
            blk[:, j] = np.partition(d[:, a:a + n], k - 1, axis=1)[:, :k].mean(1)
    return out


def _mp_cache_path(index, params, seed):
    tag = hashlib.md5(json.dumps({"p": params, "s": seed, "b": index.backbone,
                                  "n": int(index.E.shape[0])}, sort_keys=True).encode()).hexdigest()[:10]
    return config.CX_CACHE / f"multiproto_{index.backbone}_{tag}.npz"


def build_multiproto(index, params=None, seed=config.SEED, use_cache=True, verbose=True):
    """k-means sub-prototypes per EC, kept only when the split is well separated (silhouette)."""
    from sklearn.cluster import KMeans
    from sklearn.metrics import silhouette_score
    params = dict(config.MULTIPROTO if params is None else params)
    path = _mp_cache_path(index, params, seed)
    if use_cache and path.is_file():
        z = np.load(path)
        return {k: z[k] for k in z.files}
    protos, owner = [], []
    nk = np.ones(len(index.ecs), dtype=np.int64)
    sil = np.full(len(index.ecs), np.nan)
    for j, (a, n) in enumerate(zip(index.starts, index.sizes)):
        centers = index.centroids[j:j + 1]
        if n >= params["min_n"]:
            X = index.E[a:a + n]
            best = None
            for k in range(2, params["kmax"] + 1):
                if n < k * params["min_size"]:
                    break
                km = KMeans(n_clusters=k, n_init=10, random_state=seed).fit(X)
                counts = np.bincount(km.labels_, minlength=k)
                if counts.min() < max(params["min_size"], int(np.ceil(params["min_frac"] * n))):
                    continue
                try:
                    s = silhouette_score(X, km.labels_, sample_size=min(int(n), 1000), random_state=seed)
                except ValueError:
                    continue
                if best is None or s > best[0]:
                    best = (s, k, km.cluster_centers_)
            if best is not None:
                sil[j] = best[0]
                if best[0] >= params["sil_min"]:
                    centers = best[2]
                    nk[j] = best[1]
        protos.append(centers)
        owner.extend([j] * len(centers))
        if verbose and (j + 1) % 1000 == 0:
            print(f"  multiproto: {j + 1}/{len(index.ecs)} ECs", flush=True)
    P = np.concatenate(protos)
    owner = np.asarray(owner)
    starts = np.searchsorted(owner, np.arange(len(index.ecs)))
    out = {"P": P, "owner": owner, "starts": starts, "nk": nk, "sil": sil}
    config.CX_CACHE.mkdir(parents=True, exist_ok=True)
    np.savez(path, **out)
    return out


def multiproto_distances(Q, mp, chunk=256):
    P = mp["P"]
    P_sq = np.einsum("ij,ij->i", P, P)
    out = np.empty((len(Q), len(mp["starts"])))
    for s in range(0, len(Q), chunk):
        d2 = _sqdist(Q[s:s + chunk], P, P_sq)
        out[s:s + chunk] = np.sqrt(np.minimum.reduceat(d2, mp["starts"], axis=1))
    return out


def balanced_distances(Q, index, k=1, lam=0.0, base=None):
    """Class-balanced local-density score: d_k(q, c) * n_c**lam.

    With class-conditional kNN density estimates p_c(q) ~ k / (n_c * r_k^D), a decision with equal
    class priors ranks classes by n_c * r_k^D, i.e. by r_k * n_c^(1/D). lam = 1/D_eff is tuned on the
    calibration split; lam = 0 is plain exemplar (k=1) / local-mean (k=3) inference, whose minimum
    over many members systematically favours populous EC classes."""
    if base is None:
        base = exemplar_distances(Q, index) if k == 1 else kmean_distances(Q, index, k=k)
    return base * (index.sizes.astype(float) ** lam)[None, :]


def scores_for(method, Q, index, mp=None, lam=None):
    if method == "centroid":
        return centroid_distances(Q, index)
    if method == "exemplar":
        return exemplar_distances(Q, index)
    if method == "kmean3":
        return kmean_distances(Q, index, k=3)
    if method == "multiproto":
        return multiproto_distances(Q, mp if mp is not None else build_multiproto(index))
    if method in ("balanced1", "balanced3"):
        if lam is None:
            raise ValueError(f"{method} needs a tuned lambda")
        return balanced_distances(Q, index, k=int(method[-1]), lam=lam)
    raise ValueError(f"unknown method {method}")


def maxsep_predict(S, ecs, ids, topn=10):
    """CLEAN max-separation on any distance-like score matrix (smaller = better)."""
    from CLEAN.evaluate import maximum_separation
    order = np.argsort(S, axis=1, kind="stable")[:, :topn]
    preds = {}
    for i, q in enumerate(ids):
        o = order[i]
        try:
            k = int(maximum_separation(S[i, o], True, False))
        except IndexError:
            k = 0
        preds[q] = [(ecs[j], float(S[i, j])) for j in o[:k + 1]]
    return preds


def top_candidates(S, ecs, k=10):
    order = np.argsort(S, axis=1, kind="stable")[:, :k]
    return order, np.take_along_axis(S, order, axis=1)
