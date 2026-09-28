"""Known-answer tests for cleanx.blast (BLASTp baseline): hit ranking and ties, maximum identity over the shorter
sequence, reduction, deterministic gzip, and sequence-set digests."""
import tempfile
from pathlib import Path

from cleanx import blast as BL

LINES = [
    "q1\ts1\t50.000\t100\t200\t150\t120\t1e-20\t50",  # pid_short 50/150
    "q1\ts2\t90.000\t40\t200\t40\t80.0\t1e-10\t36",   # short local hit: pid_short 36/40 = 0.9 (the maximum)
    "q1\ts3\t60.000\t100\t200\t300\t120\t1e-25\t60",  # ties s1 on bits, lower E-value -> top hit
    "q2\ts4\t30.000\t80\t100\t90\t40.0\t1e-3\t24",
]


def test_summary_known_answer():
    rows = {r["Entry"]: r for r in BL.summarize(BL.parse_lines(LINES), ["q1", "q2", "q3"])}
    q1, q2, q3 = rows["q1"], rows["q2"], rows["q3"]
    assert (q1["top_hit"], q1["top_score"], q1["top_evalue"]) == ("s3", 120.0, 1e-25)
    assert abs(q1["top_hit_pid_short"] - 60 / 200) < 1e-12
    assert abs(q1["max_pid_short"] - 0.9) < 1e-12 and abs(q1["max_pid_aln"] - 0.9) < 1e-12
    assert q2["top_hit"] == "s4" and abs(q2["max_pid_short"] - 24 / 90) < 1e-12
    assert q3["top_hit"] == "" and q3["max_pid_short"] == 0.0  # no hit: the baseline abstains


def test_reduce_keeps_top_hits_and_the_maximum_identity_hit():
    red = BL.reduce_hits(BL.parse_lines(LINES), keep=1)
    assert list(zip(red.qseqid, red.sseqid)) == [("q1", "s3"), ("q1", "s2"), ("q2", "s4")]
    full = {r["Entry"]: r for r in BL.summarize(BL.parse_lines(LINES), ["q1", "q2"])}
    part = {r["Entry"]: r for r in BL.summarize(red, ["q1", "q2"])}
    assert full == part  # reduction loses nothing cx_04 uses


def test_write_read_round_trip_is_byte_identical():
    df = BL.parse_lines(LINES)
    with tempfile.TemporaryDirectory() as td:
        a, b = Path(td) / "a.tsv.gz", Path(td) / "b.tsv.gz"
        BL.write_hits(df, a)
        BL.write_hits(BL.read_hits(a), b)
        assert a.read_bytes() == b.read_bytes()  # original text kept, gzip without timestamp
        assert list(BL.read_hits(a).sseqid) == ["s1", "s2", "s3", "s4"]


def test_sequence_digest_ignores_order_and_detects_changes():
    s = {"a": "MKV", "b": "MLL"}
    assert BL.seq_digest(["a", "b"], s) == BL.seq_digest(["b", "a"], s)
    assert BL.seq_digest(["a", "b"], s) != BL.seq_digest(["a", "b"], {"a": "MKV", "b": "MLI"})
