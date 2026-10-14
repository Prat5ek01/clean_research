"""Run BLASTp and verify sequence digests for the homology baseline.

See docs/code_reference.md#cx-11-blast."""
import argparse
import json
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone

from cleanx import config

config.enter_app_dir()

import pandas as pd

from cleanx import blast as BL
from cleanx import labels as L

BLAST_DIR = config.CX_IDENTITY / "blast"
DB_DIR = config.CX_CACHE / "blastdb"
MANIFEST = BLAST_DIR / "manifest.json"
SETS = ["new392", "price149", "sp_temporal_cal", "sp_temporal_test"]
PAIR = "sp_temporal_test:sp_temporal_cal"
fails = []


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f" -- {detail}" if detail else ""), flush=True)
    if not ok:
        fails.append(name)
    return ok


def training():
    df = pd.read_csv(config.BACKBONES["split100"]["csv"], sep="\t", dtype=str, keep_default_na=False)
    return list(df.Entry), dict(zip(df.Entry, df.Sequence))


def load_set(name):
    ids, _, seqs = L.load_named_dataset(name)
    return ids, seqs


def file_digest(path):
    import hashlib
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(cmd):
    r = subprocess.run([str(c) for c in cmd], capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"command failed ({r.returncode}): {' '.join(map(str, cmd))}\n{r.stderr[-2000:]}")
    return r


def blast_version():
    exe = shutil.which("blastp")
    if exe is None or shutil.which("makeblastdb") is None:
        return None
    return run(["blastp", "-version"]).stdout.splitlines()[0].strip()


def make_db(name, ids, seqs):
    DB_DIR.mkdir(parents=True, exist_ok=True)
    fasta = DB_DIR / f"{name}.fasta"
    L.write_fasta(fasta, ids, seqs)
    run(["makeblastdb", "-in", fasta, "-dbtype", "prot", "-out", DB_DIR / name])
    return DB_DIR / name


def search(q_ids, q_seqs, db, target_ids, tag, threads):
    """blastp -> reduced hit table (DataFrame) and the exact command."""
    DB_DIR.mkdir(parents=True, exist_ok=True)
    qfa, raw = DB_DIR / f"q_{tag}.fasta", DB_DIR / f"q_{tag}.tsv"
    L.write_fasta(qfa, q_ids, q_seqs)
    cmd = ["blastp", "-query", qfa, "-db", db, "-outfmt", BL.OUTFMT, "-evalue", BL.EVALUE,
           "-max_target_seqs", BL.MAX_TARGET_SEQS, "-num_threads", threads, "-out", raw]
    t0 = time.time()
    run(cmd)
    df = BL.parse_lines(raw.read_text(encoding="utf-8").splitlines())
    raw.unlink()
    qfa.unlink()
    unknown = set(df.sseqid) - set(target_ids)
    if unknown:
        raise SystemExit(f"{tag}: {len(unknown)} subject ids are not target accessions (e.g. {sorted(unknown)[:3]});"
                         f" BLAST rewrote the ids")
    stray = set(df.qseqid) - set(q_ids)
    if stray:
        raise SystemExit(f"{tag}: {len(stray)} query ids not in the query set (e.g. {sorted(stray)[:3]})")
    print(f"  {tag}: {len(q_ids)} queries, {df.qseqid.nunique()} with hits, {len(df)} HSPs "
          f"({time.time() - t0:.0f}s)", flush=True)
    rel = [str(c).replace(str(config.APP_DIR) + "/", "").replace(str(config.APP_DIR) + "\\", "") for c in cmd]
    return BL.reduce_hits(df), rel


def verify(manifest):
    """Digests of the local sequences and of the stored tables must equal the manifest."""
    t_ids, t_seqs = training()
    check("training sequences (split100) are the ones searched", BL.seq_digest(t_ids, t_seqs) ==
          manifest["database"]["digest"], f"{len(t_ids)} sequences")
    for name, info in manifest["sets"].items():
        ids, seqs = load_set(name)
        check(f"{name}: local sequences are the ones searched", BL.seq_digest(ids, seqs) == info["digest"],
              f"{len(ids)} proteins")
        path = BLAST_DIR / info["file"]
        check(f"{name}: stored hit table intact", path.is_file() and file_digest(path) == info["file_sha256"],
              info["file"])
    pair = manifest.get("pair")
    if pair:
        for side in ("query", "target"):
            ids, seqs = load_set(pair[side])
            check(f"pair {pair['query']} vs {pair['target']}: {side} sequences are the ones searched",
                  BL.seq_digest(ids, seqs) == pair[f"{side}_digest"])
        path = BLAST_DIR / pair["file"]
        check("pair: stored hit table intact", path.is_file() and file_digest(path) == pair["file_sha256"],
              pair["file"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sets", nargs="+", default=SETS)
    ap.add_argument("--pair", default=PAIR)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--verify", action="store_true", help="check stored tables against local data; no BLAST")
    ap.add_argument("--smoke", action="store_true", help="verify + re-search 5 New-392 queries if blastp exists")
    args = ap.parse_args()
    config.ensure_dirs()
    t0 = time.time()

    if args.verify or args.smoke:
        if not MANIFEST.is_file():
            check("BLAST manifest present", False, f"{MANIFEST} missing; run cx_11_blast.py where BLAST+ is installed")
        else:
            manifest = L.load_json(MANIFEST)
            print(f"BLAST tables from {manifest['blast_version']} ({manifest['created_utc']})")
            verify(manifest)
            if args.smoke:
                ver = blast_version()
                if ver is None:
                    print("  SKIP re-search: blastp/makeblastdb not installed here (verification above still applies)")
                else:
                    t_ids, t_seqs = training()
                    db = make_db("split100", t_ids, t_seqs)
                    ids, seqs = load_set("new392")
                    ids = ids[:5]
                    red, _ = search(ids, seqs, db, t_ids, "new392_smoke", args.threads)
                    BL.write_hits(red, BLAST_DIR / "smoke" / "new392_smoke.tsv.gz")
                    new = pd.DataFrame(BL.summarize(red, ids)).set_index("Entry")
                    old = pd.DataFrame(BL.summarize(BL.read_hits(BLAST_DIR / manifest["sets"]["new392"]["file"]),
                                                    ids)).set_index("Entry")
                    same = (new.top_hit.equals(old.top_hit) and (new.top_score - old.top_score).abs().max() < 1e-9
                            and (new.max_pid_short - old.max_pid_short).abs().max() < 1e-12)
                    check("re-search of 5 New-392 queries reproduces the stored top hits and identities", same,
                          ver)
        tag = "CX11 SMOKE" if args.smoke else "CX11 VERIFY"
        print(f"\n{tag} {'FAILED' if fails else 'DONE'} in {time.time() - t0:.0f}s")
        sys.exit(1 if fails else 0)

    ver = blast_version()
    if ver is None:
        raise SystemExit("blastp/makeblastdb not found. Install NCBI BLAST+ (e.g. apt-get install ncbi-blast+) or use "
                         "the committed tables with cx_04_identity.py --tool blast")
    print(f"{ver}; {args.threads} threads")
    t_ids, t_seqs = training()
    db = make_db("split100", t_ids, t_seqs)
    manifest = {"blast_version": ver, "created_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
                "parameters": {"outfmt": BL.OUTFMT, "evalue": BL.EVALUE, "max_target_seqs": BL.MAX_TARGET_SEQS,
                               "kept_per_query": f"{BL.KEEP} best by bit score + maximum-identity hit",
                               "scoring": "blastp defaults (BLOSUM62, gap open 11, extend 1, comp. stats 2)",
                               "makeblastdb": "-dbtype prot (no -parse_seqids, so subject ids stay verbatim)"},
                "database": {"name": "split100", "n": len(t_ids), "digest": BL.seq_digest(t_ids, t_seqs)},
                "sets": {}, "commands": {}}
    for name in args.sets:
        ids, seqs = load_set(name)
        red, cmd = search(ids, seqs, db, t_ids, name, args.threads)
        path = BLAST_DIR / f"{name}.tsv.gz"
        BL.write_hits(red, path)
        manifest["sets"][name] = {"n": len(ids), "digest": BL.seq_digest(ids, seqs), "file": path.name,
                                  "file_sha256": file_digest(path), "n_with_hits": int(red.qseqid.nunique()),
                                  "rows_kept": int(len(red))}
        manifest["commands"][name] = " ".join(cmd)
    if args.pair:
        qn, tn = args.pair.split(":")
        q_ids, q_seqs = load_set(qn)
        tt_ids, tt_seqs = load_set(tn)
        if set(q_ids) & set(tt_ids):
            raise SystemExit(f"{qn} and {tn} share accessions; they must be disjoint")
        pdb = make_db(tn, tt_ids, tt_seqs)
        red, cmd = search(q_ids, q_seqs, pdb, tt_ids, f"{qn}_vs_{tn}", args.threads)
        path = BLAST_DIR / f"{qn}_vs_{tn}.tsv.gz"
        BL.write_hits(red, path)
        manifest["pair"] = {"query": qn, "target": tn, "query_digest": BL.seq_digest(q_ids, q_seqs),
                            "target_digest": BL.seq_digest(tt_ids, tt_seqs), "file": path.name,
                            "file_sha256": file_digest(path), "n_with_hits": int(red.qseqid.nunique())}
        manifest["commands"]["pair"] = " ".join(cmd)
    L.save_json(manifest, MANIFEST)
    print(f"\nCX11 DONE in {time.time() - t0:.0f}s -> {BLAST_DIR}")


if __name__ == "__main__":
    main()
