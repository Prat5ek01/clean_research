"""Run reproducible BLASTp searches and reduce hits for baseline evaluation.

Rank hits by bit score, E-value and accession. Identity uses the shorter sequence.

See docs/code_reference.md#blast."""
import gzip
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd

OUTFMT = "6 qseqid sseqid pident length qlen slen bitscore evalue nident"
COLS = ["qseqid", "sseqid", "pident", "length", "qlen", "slen", "bitscore", "evalue", "nident"]
NUM = ["pident", "length", "qlen", "slen", "bitscore", "evalue", "nident"]
EVALUE = 10.0
MAX_TARGET_SEQS = 500
KEEP = 10


def seq_digest(ids, seqs):
    """SHA-256 over the sorted (id, sequence) pairs: identifies a sequence set independently of file format,
    row order and line endings."""
    h = hashlib.sha256()
    for i in sorted(ids):
        h.update(f"{i}\t{seqs[i]}\n".encode("utf-8"))
    return h.hexdigest()


def parse_lines(lines):
    """BLAST tabular lines (OUTFMT) -> DataFrame: original text in COLS plus numeric x_<col>, pid_short, pid_aln."""
    recs = [ln.rstrip("\r\n").split("\t") for ln in lines if ln.strip()]
    bad = [r for r in recs if len(r) != len(COLS)]
    if bad:
        raise ValueError(f"BLAST output has {len(bad)} lines without {len(COLS)} columns (expected -outfmt "
                         f"'{OUTFMT}'), e.g. {bad[0]}")
    df = pd.DataFrame(recs, columns=COLS)
    for c in NUM:
        df[f"x_{c}"] = df[c].astype(float) if len(df) else pd.Series(dtype=float)
    df["pid_short"] = df.x_nident / np.minimum(df.x_qlen, df.x_slen).clip(lower=1)
    df["pid_aln"] = df.x_nident / df.x_length.clip(lower=1)
    return df


def read_hits(path):
    raw = Path(path).read_bytes()
    text = gzip.decompress(raw).decode("utf-8") if str(path).endswith(".gz") else raw.decode("utf-8")
    return parse_lines(text.splitlines())


def rank_hits(df):
    """Hits sorted best-first within each query: bit score desc, E-value asc, subject id asc."""
    return df.sort_values(["qseqid", "x_bitscore", "x_evalue", "sseqid"], ascending=[True, False, True, True])


def reduce_hits(df, keep=KEEP):
    """Per query: the `keep` best hits by bit score plus the maximum-identity hit (if not already kept)."""
    r = rank_hits(df)
    top = r.groupby("qseqid", sort=False).head(keep)
    best_pid = r.sort_values(["qseqid", "pid_short", "x_bitscore", "sseqid"],
                             ascending=[True, False, False, True]).groupby("qseqid", sort=False).head(1)
    out = pd.concat([top, best_pid])
    out = out[~out.index.duplicated(keep="first")]
    return rank_hits(out)


def summarize(df, q_ids):
    """Per query id (order kept): max_pid_short, max_pid_aln, top_hit, top_score (bits), top_evalue,
    top_hit_pid_short. Queries without any hit get zeros and no top hit (the baseline abstains)."""
    r = rank_hits(df)
    first = r.groupby("qseqid", sort=False).head(1).set_index("qseqid")
    mx = df.groupby("qseqid").agg(max_pid_short=("pid_short", "max"), max_pid_aln=("pid_aln", "max"))
    rows = []
    for q in q_ids:
        if q not in first.index:
            rows.append({"Entry": q, "max_pid_short": 0.0, "max_pid_aln": 0.0, "top_hit": "", "top_score": 0.0,
                         "top_evalue": np.nan, "top_hit_pid_short": 0.0})
            continue
        t = first.loc[q]
        rows.append({"Entry": q, "max_pid_short": float(mx.loc[q, "max_pid_short"]),
                     "max_pid_aln": float(mx.loc[q, "max_pid_aln"]), "top_hit": str(t.sseqid),
                     "top_score": float(t.x_bitscore), "top_evalue": float(t.x_evalue),
                     "top_hit_pid_short": float(t.pid_short)})
    return rows


def write_hits(df, path):
    """Kept lines, original text, gzip with mtime 0 so identical content gives identical bytes."""
    text = "".join("\t".join(row) + "\n" for row in df[COLS].itertuples(index=False, name=None))
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_bytes(gzip.compress(text.encode("utf-8"), mtime=0))
