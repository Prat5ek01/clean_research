"""Maximum sequence identity of query proteins to the CLEAN training set, and a homology-transfer
baseline (EC of the best local-alignment hit).

pid_short = identical aligned residues / length of the shorter sequence (conservative against short
local hits); pid_aln = identical / alignment columns (MMseqs2 'fident')."""
import time

import numpy as np
from scipy import sparse

_AA = "ACDEFGHIKLMNPQRSTVWY"
_LUT = np.full(256, 20, dtype=np.int64)
for _i, _a in enumerate(_AA):
    _LUT[ord(_a)] = _i
_NK = 21 ** 3


def _kmers(seq):
    a = _LUT[np.frombuffer(seq.encode("ascii", "replace"), dtype=np.uint8)]
    if len(a) < 3:
        return np.empty(0, dtype=np.int64)
    return np.unique(a[:-2] * 441 + a[1:-1] * 21 + a[2:])


def _kmer_matrix(seqs):
    parts = [_kmers(s) for s in seqs]
    lens = np.array([len(p) for p in parts], dtype=np.int64)
    indptr = np.concatenate([[0], np.cumsum(lens)])
    idx = np.concatenate(parts) if parts else np.empty(0, np.int64)
    M = sparse.csr_matrix((np.ones(len(idx), np.float32), idx, indptr), shape=(len(seqs), _NK))
    return M, lens


def kmer_candidates(q_seqs, t_seqs, topn=30, t_chunk=40000, q_chunk=256, verbose=True):
    """Top-n training candidates per query by cosine similarity of 3-mer sets."""
    Mq, nq = _kmer_matrix(q_seqs)
    best_i = np.full((len(q_seqs), topn), -1, dtype=np.int64)
    best_s = np.full((len(q_seqs), topn), -1.0)
    t0 = time.time()
    for ts in range(0, len(t_seqs), t_chunk):
        Mt, nt = _kmer_matrix(t_seqs[ts:ts + t_chunk])
        MtT = Mt.T.tocsr()
        for qs in range(0, len(q_seqs), q_chunk):
            S = (Mq[qs:qs + q_chunk] @ MtT).toarray()
            S /= np.sqrt(np.maximum(nq[qs:qs + q_chunk, None], 1) * np.maximum(nt[None, :], 1))
            k = min(topn, S.shape[1])
            cand = np.argpartition(-S, k - 1, axis=1)[:, :k]
            cs = np.take_along_axis(S, cand, 1)
            allc = np.concatenate([best_i[qs:qs + q_chunk], cand + ts], axis=1)
            alls = np.concatenate([best_s[qs:qs + q_chunk], cs], axis=1)
            keep = np.argsort(-alls, axis=1, kind="stable")[:, :topn]
            best_i[qs:qs + q_chunk] = np.take_along_axis(allc, keep, 1)
            best_s[qs:qs + q_chunk] = np.take_along_axis(alls, keep, 1)
        if verbose:
            print(f"  k-mer prefilter: {min(ts + t_chunk, len(t_seqs))}/{len(t_seqs)} training seqs "
                  f"({time.time() - t0:.0f}s)", flush=True)
    return best_i, best_s


def _aligner():
    from Bio import Align
    from Bio.Align import substitution_matrices
    al = Align.PairwiseAligner()
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score = -11
    al.extend_gap_score = -1
    al.mode = "local"
    return al


def _clean(seq):
    return "".join(c if c in _AA else "X" for c in seq.upper())


def align_stats(al, q, t):
    """(score, identical, aln_columns) of the best local alignment."""
    aln = al.align(q, t)[0]
    segs_q, segs_t = aln.aligned
    ident = 0
    cols = 0
    prev_qe = prev_te = None
    for (qs, qe), (ts, te) in zip(segs_q, segs_t):
        a = np.frombuffer(q[qs:qe].encode(), dtype=np.uint8)
        b = np.frombuffer(t[ts:te].encode(), dtype=np.uint8)
        ident += int((a == b).sum())
        cols += qe - qs
        if prev_qe is not None:
            cols += (qs - prev_qe) + (ts - prev_te)
        prev_qe, prev_te = qe, te
    return float(aln.score), ident, cols


def python_search(q_ids, q_seqs, t_ids, t_seqs, topn=30, verbose=True):
    """Returns list of dicts per query: max_pid_short, max_pid_aln, top_hit, top_score."""
    q_seqs = [_clean(s) for s in q_seqs]
    t_clean = [_clean(s) for s in t_seqs]
    cand, _ = kmer_candidates(q_seqs, t_clean, topn=topn, verbose=verbose)
    al = _aligner()
    rows = []
    t0 = time.time()
    for i, (qid, q) in enumerate(zip(q_ids, q_seqs)):
        best = {"Entry": qid, "max_pid_short": 0.0, "max_pid_aln": 0.0, "top_hit": "", "top_score": 0.0,
                "top_hit_pid_short": 0.0}
        for j in cand[i]:
            if j < 0:
                continue
            t = t_clean[j]
            score, ident, cols = align_stats(al, q, t)
            pid_s = ident / max(1, min(len(q), len(t)))
            pid_a = ident / max(1, cols)
            best["max_pid_short"] = max(best["max_pid_short"], pid_s)
            best["max_pid_aln"] = max(best["max_pid_aln"], pid_a)
            if score > best["top_score"]:
                best.update(top_hit=t_ids[j], top_score=score, top_hit_pid_short=pid_s)
        rows.append(best)
        if verbose and (i + 1) % 200 == 0:
            print(f"  aligned {i + 1}/{len(q_ids)} queries ({time.time() - t0:.0f}s)", flush=True)
    return rows


def parse_m8(path, q_ids):
    """MMseqs2 easy-search output with
    --format-output query,target,fident,alnlen,qlen,tlen,bits,evalue"""
    import pandas as pd
    cols = ["query", "target", "fident", "alnlen", "qlen", "tlen", "bits", "evalue"]
    df = pd.read_csv(path, sep="\t", header=None, names=cols)
    df["ident"] = (df.fident * df.alnlen).round()
    df["pid_short"] = df.ident / np.minimum(df.qlen, df.tlen).clip(lower=1)
    rows = []
    g = df.groupby("query")
    for q in q_ids:
        if q not in g.groups:
            rows.append({"Entry": q, "max_pid_short": 0.0, "max_pid_aln": 0.0, "top_hit": "",
                         "top_score": 0.0, "top_hit_pid_short": 0.0})
            continue
        d = g.get_group(q)
        top = d.loc[d.bits.idxmax()]
        rows.append({"Entry": q, "max_pid_short": float(d.pid_short.max()), "max_pid_aln": float(d.fident.max()),
                     "top_hit": str(top.target), "top_score": float(top.bits),
                     "top_hit_pid_short": float(top.pid_short)})
    return rows


def identity_bin(pid):
    if pid < 0.30:
        return "<30%"
    if pid < 0.50:
        return "30-50%"
    if pid < 0.70:
        return "50-70%"
    if pid < 0.90:
        return "70-90%"
    return ">=90%"


IDENTITY_BINS = ["<30%", "30-50%", "50-70%", "70-90%", ">=90%"]
