import os
from pathlib import Path

APP_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = APP_DIR / "data"
PRETRAINED_DIR = DATA_DIR / "pretrained"
ESM_DIR = DATA_DIR / "esm_data"
CX_DATA = DATA_DIR / "cleanx"
CX_DOWNLOADS = CX_DATA / "downloads"
CX_SETS = CX_DATA / "sets"
CX_FASTA = CX_DATA / "fasta"
CX_IDENTITY = CX_DATA / "identity"
CX_CACHE = CX_DATA / "cache"
CX_RESULTS = APP_DIR / "results" / "cleanx"
CX_PREDS = CX_RESULTS / "preds"
CX_FIGS = CX_RESULTS / "figures"

GMM_PKL = PRETRAINED_DIR / "gmm_ensumble.pkl"

BACKBONES = {
    "split100": {"csv": DATA_DIR / "split100.csv", "emb": PRETRAINED_DIR / "100.pt",
                 "pth": PRETRAINED_DIR / "split100.pth", "out_dim": 128},
    "split70": {"csv": DATA_DIR / "split70.csv", "emb": PRETRAINED_DIR / "70.pt",
                "pth": PRETRAINED_DIR / "split70.pth", "out_dim": 128},
}

LEGACY_SETS = {
    "new392": DATA_DIR / "new.csv",
    "price149": DATA_DIR / "datasets" / "price.csv",
}

# Tier thresholds reproduce data/ec_train_frequency.csv (verified by cx_00_check.py).
TIER_ORDER = ["unseen", "singleton", "very_rare", "rare", "common"]
TIER_RANK = {t: i for i, t in enumerate(TIER_ORDER)}


def tier_of_count(n):
    if n <= 0:
        return "unseen"
    if n == 1:
        return "singleton"
    if n <= 3:
        return "very_rare"
    if n <= 10:
        return "rare"
    return "common"


# split100 = UniProt release 2022_01; its newest regular entries were created 2022-02-23.
SNAPSHOT_DATE = "2022-03-01"
MAX_LEN = 1022
MIN_LEN = 30
SEED = 20260927

MULTIPROTO = {"min_n": 20, "kmax": 4, "sil_min": 0.25, "min_frac": 0.10, "min_size": 3}
KMEAN_K = 3


def dataset_csv(name):
    if name in LEGACY_SETS:
        return LEGACY_SETS[name]
    return CX_SETS / f"{name}.csv"


def is_temporal(name):
    return name not in LEGACY_SETS


def ensure_dirs():
    for d in (CX_DATA, CX_DOWNLOADS, CX_SETS, CX_FASTA, CX_IDENTITY, CX_CACHE,
              CX_RESULTS, CX_PREDS, CX_FIGS):
        d.mkdir(parents=True, exist_ok=True)


def enter_app_dir():
    """CLEAN's own helpers use paths relative to app/, so every script runs from there."""
    os.chdir(APP_DIR)
