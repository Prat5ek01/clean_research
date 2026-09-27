"""ESM-1b mean embeddings for the built sets (resumable), plus a mandatory equivalence smoke test
against the embeddings CLEAN already used for New-392."""
import argparse
import shutil
import sys
import time
import warnings

warnings.filterwarnings("ignore")

from cleanx import config

config.enter_app_dir()

from cleanx import labels as L
from cleanx.esm_embed import compare_embeddings, embed_records, load_model

VERIFY_DIR = config.CX_CACHE / "esm_verify"


def verify(k, threads):
    ids, _, seqs = L.read_dataset(config.LEGACY_SETS["new392"])
    ids = sorted(ids, key=lambda i: len(seqs[i]))
    pick = [ids[0], ids[len(ids) // 2], ids[-1]][:k] if k <= 3 else ids[:: max(1, len(ids) // k)][:k]
    if VERIFY_DIR.exists():
        shutil.rmtree(VERIFY_DIR)
    model, alphabet = load_model(threads)
    t0 = time.time()
    embed_records([(i, seqs[i]) for i in pick], VERIFY_DIR, model, alphabet, overwrite=True)
    sec = (time.time() - t0) / len(pick)
    mad, mcos = compare_embeddings(pick, VERIFY_DIR, config.ESM_DIR)
    ok = mad < 1e-3 and mcos > 0.99999
    lens = [len(seqs[i]) for i in pick]
    print(f"[{'PASS' if ok else 'FAIL'}] recomputed ESM-1b embeddings match CLEAN's stored ones for "
          f"{len(pick)} New-392 proteins (lengths {lens}): max |diff| {mad:.2e}, min cosine {mcos:.7f}")
    print(f"       CPU speed: {sec:.1f} s/protein at median length ~{sorted(lens)[len(lens) // 2]} aa")
    return ok, model, alphabet


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sets", nargs="+", default=["sp_temporal_cal", "sp_temporal_test", "sp_nonenzyme"])
    ap.add_argument("--verify", type=int, default=3, help="New-392 proteins re-embedded as the equivalence test")
    ap.add_argument("--verify-only", action="store_true")
    ap.add_argument("--limit", type=int, default=0, help="embed only the first N proteins per set (smoke)")
    ap.add_argument("--threads", type=int, default=0)
    ap.add_argument("--toks-per-batch", type=int, default=4096)
    args = ap.parse_args()
    config.ensure_dirs()
    if args.verify > 0:
        ok, model, alphabet = verify(args.verify, args.threads)
        if not ok:
            print("Embedding equivalence FAILED - do not embed new sets until this passes.")
            sys.exit(1)
    else:
        print("  WARNING: equivalence test skipped (--verify 0)")
        model, alphabet = load_model(args.threads)
    if args.verify_only:
        return
    for name in args.sets:
        fa = config.CX_FASTA / f"{name}.fasta"
        if not fa.is_file():
            print(f"  skip {name}: {fa} not found (run cx_02_build_temporal.py)")
            continue
        recs = L.read_fasta(fa)
        if args.limit:
            recs = recs[:args.limit]
        print(f"== {name}: {len(recs)} proteins")
        embed_records(recs, config.ESM_DIR, model, alphabet, toks_per_batch=args.toks_per_batch)
    missing = []
    for name in args.sets:
        fa = config.CX_FASTA / f"{name}.fasta"
        if fa.is_file():
            recs = L.read_fasta(fa)[:args.limit or None]
            missing += [i for i, _ in recs if not (config.ESM_DIR / f"{i}.pt").is_file()]
    print(f"\nCX03 {'DONE' if not missing else 'INCOMPLETE'}: {len(missing)} embeddings missing")
    sys.exit(1 if missing else 0)


if __name__ == "__main__":
    main()
