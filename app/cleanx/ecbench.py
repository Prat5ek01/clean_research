"""EC-Bench (Davoudi et al., Bioinformatics Advances 2026, doi:10.1093/bioadv/vbag004; MIT licence) publishes
the per-protein predictions of 16 EC predictors on its test set (468 proteins) and on Price-149, for models
trained on its 100%- and 30%-identity clusterings. E5 re-scores them. Files are fetched from a pinned commit and
recorded with their SHA-256; labels are parsed exactly as in EC-Bench's code/results.ipynb (tokens split on ',',
stripped, '-' = no prediction, EC strings truncated to the evaluated depth)."""
import csv
import hashlib
from pathlib import Path

from . import config
from . import downloads as dl

COMMIT = "5e98a2c632c6644f06975229ededbd09f4866809"
RAW = "https://raw.githubusercontent.com/dsaeedeh/EC-Bench/{commit}/data/{name}"
FILES = {("test", "100"): "results_test_100_go_3.csv", ("test", "30"): "results_test_30_go_3.csv",
         ("price149", "100"): "results_price_149_100_go_3.csv", ("price149", "30"): "results_price_149_30_go_3.csv"}
DEST = config.CX_DOWNLOADS / "ec-bench"
# SHA-256 of the four files at COMMIT (first downloaded 2026-09-28): the confirmatory run must read these bytes.
SHA256 = {"results_test_100_go_3.csv": "fb13fce97e81d79948b844cccdba9ee5e9cf47dfc5c9f0205a79ee1752234f87",
          "results_test_30_go_3.csv": "2cf5964a9b678d29d17ffb533494daeec5c495bf650de2c53743ce809ee8d72a",
          "results_price_149_100_go_3.csv": "099a0d83f84978555ddf27d1b429e502e853d1cc416e0f00695a4e54ce319a82",
          "results_price_149_30_go_3.csv": "157c79b73c3cf7abc862ac18346c5e5f4890d79e91dd012bf283e90a93b4f259"}

# The 16 methods EC-Bench plots, in its order (code/results.ipynb, f1_scores()).
METHODS = ["catfam", "deepec", "deepectransformer-regular", "deepectransformer-learnt", "ecpred", "ecrecer", "priam",
           "blastp", "clean-supcon", "proteinbert-regular", "proteinbert-learnt", "enzbert-regular", "enzbert-learnt",
           "majority-voting-regular", "majority-voting-learnt", "stacked"]

# Reproduction gate: support-weighted F1 at EC levels 1 and 4 on results_test_100_go_3.csv, as printed in the
# stored output of code/results.ipynb (cell calling plot(file, 'Test', '30', 'weighted')) at COMMIT.
PRINTED_TEST100 = {
    1: [0.1639364228413389, 0.34876687401220874, 0.6112500136413195, 0.5948146183556976, 0.7830016110057062,
        0.8435933182498203, 0.20830336717317194, 0.8419759037391658, 0.8378788625442523, 0.7400707870838628,
        0.7469782684929972, 0.8703182296293466, 0.8734297810338719, 0.8760049147636545, 0.8922915025071707,
        0.9059797466498624],
    4: [0.07421960270683257, 0.12687415767582522, 0.09515205379645458, 0.11521429890781562, 0.27277972168596387,
        0.3592416909611381, 0.09901552278762495, 0.502319740055168, 0.25084899092216373, 0.3940553494693946,
        0.41420540553083074, 0.46832080998537395, 0.47452269775652095, 0.507676739711872, 0.5147264391300335,
        0.6389903621983571]}


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def fetch(dest=DEST, force=False):
    """Download the four prediction files (cached). Returns {(dataset, regime): path}."""
    dest = Path(dest)
    out = {}
    for key, name in FILES.items():
        path = dest / name
        if force or not path.is_file():
            url = RAW.format(commit=COMMIT, name=name)
            print(f"  downloading EC-Bench {name} @ {COMMIT[:7]}", flush=True)
            headers = dl._fetch(url, path)
            dl._manifest(path, url, headers, {"commit": COMMIT, "sha256": _sha256(path), "licence": "MIT"})
        out[key] = path
    return verify(out)


def verify(paths):
    """Fail loudly unless every file has its pinned SHA-256. Returns paths unchanged."""
    for path in paths.values():
        path = Path(path)
        digest = _sha256(path)
        if digest != SHA256[path.name]:
            raise ValueError(f"{path}: SHA-256 {digest[:12]}... is not the pinned {SHA256[path.name][:12]}... "
                             f"(EC-Bench commit {COMMIT[:7]}); delete the file and re-run, or use --force-download")
    return paths


def parse_cell(cell):
    """EC-Bench cell -> list of EC tokens ('-' and blanks removed, order kept, duplicates dropped)."""
    out = []
    for tok in str(cell).split(","):
        t = tok.strip()
        if t and t != "-" and t.lower() != "nan" and t not in out:
            out.append(t)
    return out


def read_results(path, methods=METHODS):
    """Returns (ids, truth, {method: {id: [EC]}}). Fails loudly on missing columns or duplicate ids."""
    with open(path, encoding="utf-8-sig", newline="") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise ValueError(f"{path}: empty")
    missing = [m for m in list(methods) + ["id", "ec_number"] if m not in rows[0]]
    if missing:
        raise ValueError(f"{path}: missing columns {missing}")
    ids = [r["id"].strip() for r in rows]
    if len(set(ids)) != len(ids):
        raise ValueError(f"{path}: duplicate ids")
    truth = {i: parse_cell(r["ec_number"]) for i, r in zip(ids, rows)}
    preds = {m: {i: parse_cell(r[m]) for i, r in zip(ids, rows)} for m in methods}
    return ids, truth, preds


def write_synthetic(path, n=60, seed=0, methods=METHODS):
    """Synthetic file in EC-Bench's format (for smoke/integration tests only)."""
    import numpy as np
    rng = np.random.default_rng(seed)
    vocab = [f"{a}.{b}.{c}.{d}" for a in range(1, 4) for b in range(1, 3) for c in range(1, 3) for d in range(1, 6)]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["id"] + list(methods) + ["ec_number"])
        for i in range(n):
            t = list(rng.choice(vocab[:25], size=int(rng.integers(1, 3)), replace=False))
            row = []
            for k, m in enumerate(methods):
                if rng.random() < 0.1:
                    row.append("-")
                    continue
                size = 1 + int(rng.integers(0, 1 + k % 4))
                p = list(rng.choice(vocab, size=size, replace=False))
                if rng.random() < 0.3 + 0.04 * k:
                    p[0] = t[0]
                row.append(",".join(p))
            w.writerow([f"S{i:04d}"] + row + [",".join(t)])
    return path
