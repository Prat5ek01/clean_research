"""ESM-1b mean embeddings computed exactly like esm/scripts/extract.py --include mean
(layer 33, residues 1..min(L,1022), toks_per_batch batching), saved in CLEAN's .pt format."""
import os
import time
from pathlib import Path

import torch

MODEL_NAME = "esm1b_t33_650M_UR50S"


def load_model(threads=None):
    import esm
    if threads:
        torch.set_num_threads(int(threads))
    print(f"  loading {MODEL_NAME} (needs ~8-10 GB free RAM while loading) ...", flush=True)
    t0 = time.time()
    model, alphabet = esm.pretrained.load_model_and_alphabet(MODEL_NAME)
    model.eval()
    print(f"  model loaded in {time.time() - t0:.0f}s; torch threads = {torch.get_num_threads()}", flush=True)
    return model, alphabet


def embed_records(records, out_dir, model=None, alphabet=None, toks_per_batch=4096, truncation=1022,
                  threads=None, overwrite=False, verbose=True):
    """records: [(id, sequence)]. Writes <out_dir>/<id>.pt; skips existing files unless overwrite."""
    from esm import FastaBatchedDataset
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    todo = [(i, s) for i, s in records if overwrite or not (out_dir / f"{i}.pt").is_file()]
    if verbose:
        print(f"  {len(records) - len(todo)} already embedded, {len(todo)} to do -> {out_dir}", flush=True)
    if not todo:
        return 0
    if model is None:
        model, alphabet = load_model(threads)
    layer = model.num_layers
    ds = FastaBatchedDataset([i for i, _ in todo], [s for _, s in todo])
    batches = ds.get_batch_indices(toks_per_batch, extra_toks_per_seq=1)
    loader = torch.utils.data.DataLoader(ds, collate_fn=alphabet.get_batch_converter(truncation),
                                         batch_sampler=batches)
    done, t0 = 0, time.time()
    with torch.no_grad():
        for b, (labels, strs, toks) in enumerate(loader):
            out = model(toks, repr_layers=[layer], return_contacts=False)
            reps = out["representations"][layer]
            for k, label in enumerate(labels):
                L = min(truncation, len(strs[k]))
                res = {"label": label, "mean_representations": {layer: reps[k, 1:L + 1].mean(0).clone()}}
                tmp = out_dir / f".{label}.pt.tmp"
                torch.save(res, tmp)
                os.replace(tmp, out_dir / f"{label}.pt")
            done += len(labels)
            if verbose and (b % 10 == 0 or done == len(todo)):
                el = time.time() - t0
                eta = el / done * (len(todo) - done)
                print(f"  batch {b + 1}/{len(batches)}: {done}/{len(todo)} seqs, {el / 60:.1f} min elapsed, "
                      f"ETA {eta / 60:.1f} min", flush=True)
    return done


def compare_embeddings(ids, dir_a, dir_b):
    """Max absolute difference and min cosine between two sets of saved mean embeddings."""
    def vec(path):
        obj = torch.load(str(path), map_location="cpu")
        return (obj["mean_representations"][33] if isinstance(obj, dict) else obj).double().reshape(-1)

    worst_abs, worst_cos = 0.0, 1.0
    for i in ids:
        a = vec(Path(dir_a) / f"{i}.pt")
        b = vec(Path(dir_b) / f"{i}.pt")
        worst_abs = max(worst_abs, float((a - b).abs().max()))
        worst_cos = min(worst_cos, float(torch.nn.functional.cosine_similarity(a, b, dim=0)))
    return worst_abs, worst_cos
