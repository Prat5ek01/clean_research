from pathlib import Path

import pandas as pd

from . import config
from . import labels as L
from .identity import identity_bin
from .nomenclature import forward_counts, forward_map, parse_enzyme_dat

ENZYME_DAT = config.CX_DOWNLOADS / "enzyme.dat"


def load_fmap(required=False):
    if not ENZYME_DAT.is_file():
        if required:
            raise FileNotFoundError(f"{ENZYME_DAT} missing; run cx_02_build_temporal.py first")
        return None
    return forward_map(parse_enzyme_dat(ENZYME_DAT))


class EvalSet:
    """A labelled evaluation set. Temporal sets carry current (post-2022) EC labels, so predictions
    from the 2022-era CLEAN vocabulary are mapped through the ENZYME transfer table before scoring."""

    def __init__(self, name, fmap=None):
        self.name = name
        self.ids, self.truth, self.seqs = L.load_named_dataset(name)
        self.temporal = config.is_temporal(name)
        self.fmap = fmap if self.temporal else None
        if self.temporal and self.fmap is None:
            self.fmap = load_fmap(required=True)

    def score_labels(self, preds):
        lab = L.labels_only(preds)
        return L.harmonize_all(lab, self.truth, self.fmap) if self.fmap else lab

    def tiers(self, index):
        counts = forward_counts(index.counts, self.fmap) if self.fmap else index.counts
        return {i: L.query_tier(self.truth[i], counts) for i in self.ids}

    def identity(self):
        p = config.CX_IDENTITY / f"{self.name}.csv"
        if not p.is_file():
            return None
        df = pd.read_csv(p, dtype={"Entry": str})
        m = dict(zip(df.Entry, df.max_pid_short))
        if any(i not in m for i in self.ids):
            return None
        return {i: identity_bin(float(m[i])) for i in self.ids}


def write_rows(rows, path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    df.to_csv(path, index=False)
    return df


def fmt(x, nd=3):
    return "nan" if x != x else f"{x:.{nd}f}"
