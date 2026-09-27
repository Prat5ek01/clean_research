from pathlib import Path

import numpy as np
import torch

from . import config


def load_esm_means(ids, esm_dir=config.ESM_DIR):
    esm_dir = Path(esm_dir)
    missing = [i for i in ids if not (esm_dir / f"{i}.pt").is_file()]
    if missing:
        raise FileNotFoundError(
            f"{len(missing)} ESM-1b embeddings missing in {esm_dir} (e.g. {missing[:5]}); "
            f"compute them with cx_03_embed.py")
    X = np.empty((len(ids), 1280), dtype=np.float32)
    for k, i in enumerate(ids):
        obj = torch.load(str(esm_dir / f"{i}.pt"), map_location="cpu")
        v = obj["mean_representations"][33] if isinstance(obj, dict) else obj
        X[k] = v.detach().float().reshape(-1).numpy()
    return X


class TrainIndex:
    """Frozen CLEAN training embeddings, grouped into contiguous per-EC row blocks.

    Row order of the released embedding file follows CLEAN.utils.get_ec_id_dict (EC keys in
    first-appearance order of the CSV); only block boundaries are recoverable, not which
    training accession sits on which row inside a block."""

    def __init__(self, backbone="split100"):
        from CLEAN.utils import get_ec_id_dict
        spec = config.BACKBONES[backbone]
        self.backbone = backbone
        self.id_ec, ec_id = get_ec_id_dict(str(spec["csv"]))
        self.ecs = list(ec_id.keys())
        self.col = {e: i for i, e in enumerate(self.ecs)}
        self.sizes = np.array([len(ec_id[e]) for e in self.ecs], dtype=np.int64)
        self.starts = np.zeros_like(self.sizes)
        self.starts[1:] = np.cumsum(self.sizes)[:-1]
        emb = torch.load(str(spec["emb"]), map_location="cpu")
        if not torch.is_tensor(emb):
            raise TypeError(f"{spec['emb']} is not a tensor")
        E = emb.detach().to(torch.float64).numpy()
        del emb
        if E.shape[0] != int(self.sizes.sum()):
            raise ValueError(f"{spec['emb']}: {E.shape[0]} rows but {int(self.sizes.sum())} (EC, id) pairs "
                             f"in {spec['csv']}; embeddings and CSV do not belong together")
        self.E = np.ascontiguousarray(E)
        self.E_sq = np.einsum("ij,ij->i", self.E, self.E)
        self.row_ec = np.repeat(np.arange(len(self.ecs)), self.sizes)
        self.centroids = np.add.reduceat(self.E, self.starts, axis=0) / self.sizes[:, None]
        self.C_sq = np.einsum("ij,ij->i", self.centroids, self.centroids)
        self.counts = dict(zip(self.ecs, self.sizes.tolist()))
        self.tiers = np.array([config.tier_of_count(int(n)) for n in self.sizes])
        self.pth = spec["pth"]
        self.out_dim = spec["out_dim"]
        self._model = None

    def model(self):
        if self._model is None:
            from CLEAN.model import LayerNormNet
            m = LayerNormNet(512, self.out_dim, "cpu", torch.float32)
            m.load_state_dict(torch.load(str(self.pth), map_location="cpu"))
            m.eval()
            self._model = m
        return self._model

    def embed(self, ids, esm_dir=config.ESM_DIR, batch=1024):
        X = load_esm_means(ids, esm_dir)
        out = []
        with torch.no_grad():
            for s in range(0, len(X), batch):
                out.append(self.model()(torch.from_numpy(X[s:s + batch])).numpy())
        return np.concatenate(out).astype(np.float64) if out else np.zeros((0, self.out_dim))

    def train_sequences(self):
        """Entry -> sequence for the backbone's training CSV (needed for identity search only)."""
        import pandas as pd
        df = pd.read_csv(config.BACKBONES[self.backbone]["csv"], sep="\t", dtype=str, keep_default_na=False)
        return dict(zip(df["Entry"], df["Sequence"]))
