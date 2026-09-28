# CLEAN Project — Deep-Dive Overview

> **Status (2026-09-28):** this overview predates the `cleanx` audit layer (`app/cleanx/`, `app/cx_*.py`).
> For the current research plan, results and pre-registrations see [`PUBLICATION_PLAN.md`](PUBLICATION_PLAN.md);
> for the commands, see [`app/cleanx/RUNBOOK.md`](app/cleanx/RUNBOOK.md).

> This is a companion document to the upstream [`README.md`](README.md). The upstream file is the
> original project's install/usage guide. This file explains **what the project actually does**,
> **how the pieces fit together**, **what has been added on top of the vanilla codebase in this
> local copy**, and **what the current findings are**. Written 2026-08-10 from a full read of the
> repo (source, notebooks, data, git history, audit outputs).

---

## 1. What this project is

**CLEAN** (**C**ontrastive **L**earning **E**nabled **A**nnotation of **E**nzymes... officially
"Contrastive Learning enabled Enzyme ANnotation") is the official code release for:

> Yu, Cui, Li, Luo, Jiang, Zhao. *"Enzyme function prediction using contrastive learning."*
> **Science** 379(6639), 1358–1363 (2023). DOI: [10.1126/science.adf2465](https://www.science.org/doi/10.1126/science.adf2465)

The problem it solves: given a protein's amino-acid sequence, predict its **Enzyme Commission
(EC) number** — the standard 4-part hierarchical code (e.g. `3.5.2.6`) that describes exactly
what chemical reaction an enzyme catalyzes. This matters because most proteins discovered by
sequencing have no experimentally verified function, and EC number prediction is the main way
computational biology assigns one at scale (drug discovery, metabolic engineering, synthetic
biology, etc.).

CLEAN's contribution over prior tools (e.g. BLAST-based annotation transfer) is that it uses
**metric/contrastive learning** to build an embedding space where sequences with the same EC
number cluster tightly together — including for **rare and previously "orphan" EC classes**
that have very few or even a single labeled training example, which is where sequence-similarity
methods historically fail.

This local checkout is a clone of the upstream [`tttianhao/CLEAN`](https://github.com/tttianhao/CLEAN)
repository, **plus a custom evaluation/audit layer** (see §5) that was added on top to
independently reproduce and stress-test the paper's reported numbers — this looks like the actual
active work in this "chem project" folder.

---

## 2. How the method works (pipeline)

```
  Protein sequence (FASTA)
          │
          ▼
  ┌───────────────────┐
  │  ESM-1b            │   Meta AI's pretrained protein language model.
  │  (frozen, 650M)    │   Turns a sequence into a 1280-dim embedding
  └───────────────────┘   (mean-pooled per-residue representations).
          │
          ▼
  ┌───────────────────┐
  │  CLEAN projection  │   A small 3-layer LayerNorm MLP (1280 → 512 → 512 → 128)
  │  head (trained)    │   trained with a contrastive loss so that sequences
  └───────────────────┘   sharing an EC number land close together and
          │                sequences with different EC numbers are pushed apart.
          ▼
  128-dim (or 256-dim for SupCon-Hard) "CLEAN embedding"
          │
          ▼
  ┌───────────────────┐
  │  EC cluster        │   For every EC number seen in training, average the
  │  centers           │   embeddings of all its training sequences into one
  └───────────────────┘   "cluster center" vector.
          │
          ▼
  ┌───────────────────┐
  │  Nearest-center    │   At inference time, embed the query sequence and
  │  EC calling         │   measure Euclidean distance to every EC cluster
  │  (max-sep / p-val)  │   center. Decide which EC number(s), if any, to call.
  └───────────────────┘
          │
          ▼
  Predicted EC number(s) + distance/confidence score
```

### 2.1 The two training losses (`app/src/CLEAN/losses.py`, `train-*.py`)

- **Triplet margin loss** (`train-triplet.py`): classic anchor/positive/negative triplet loss —
  pull anchor and same-EC positive together, push a hard-mined different-EC negative away by a
  margin. Output embedding dim = 128.
- **SupCon-Hard loss** (`train-supconH.py`, `losses.py`): a supervised-contrastive variant that
  uses *multiple* positives and hard-mined negatives per anchor in one batch (paper's better-
  performing option, especially on small training splits, at the cost of longer training). Output
  embedding dim = 256, and requires flipping `LayerNormNet(512, 128, ...)` → `LayerNormNet(512,
  256, ...)` in `infer.py` before use.

Before training, singleton EC numbers (only one training sequence) are handled by a mutation
step (`mutate_single_seq_ECs`) that generates synthetic positive examples so the contrastive loss
always has *something* to pull together.

### 2.2 The two inference algorithms (`app/src/CLEAN/infer.py`)

- **`infer_maxsep`** — deterministic, no hyperparameters. Takes the query's 10 nearest EC cluster
  centers by distance, looks for the largest gap ("maximum separation") in that sorted distance
  list, and calls every EC number above the gap. This is the primary/recommended algorithm and is
  what `CLEAN_infer_fasta.py` uses by default.
- **`infer_pvalue`** — samples `nk_random` thousand random training sequences per EC number to
  build a null/background distance distribution, then calls an EC number significant if the
  query's distance beats a `p_value` threshold. Tunable, but slower and slightly weaker in the
  paper's own numbers (confirmed again in this repo's audit, §5).

### 2.3 Confidence calibration (`app/gmm.py`, `app/src/CLEAN/uncertainty.py`)

Raw embedding distance isn't a probability. `gmm.py` fits a 2-component Gaussian Mixture Model
per EC number over (positive-class distance, hard-negative distance) samples, producing a
distance → confidence-score (0–1) mapping. The ensemble of these GMMs is shipped pretrained as
`app/data/pretrained/gmm_ensumble.pkl` and is what `CLEAN_infer_fasta.py` / `CLEAN_inference.py`
use to turn distances into the "confidence" numbers in their output CSVs (see the audit's caveat
about this in §5.3 — this repo's own analysis suggests the calibration is only loosely
discriminative).

---

## 3. Repository layout

```
CLEAN/
├── README.md                     Upstream install/usage guide (official)
├── PROJECT_OVERVIEW.md            ← this file
├── CLEAN.png                      Paper figure / method diagram
├── Dockerfile / entrypoint.sh     CPU-only container build; downloads pretrained
│                                  weights via gdown, runs CLEAN_infer_fasta.py on
│                                  data/inputs/init.fasta by default
├── NON-EXCLUSIVE RESEARCH USE
│   LICENSE FOR CLEAN SOFTWARE.pdf Research-only license (NOT permissive open source —
│                                  commercial/production use is restricted, see §6)
│
└── app/                           Everything actually runs from here
    ├── src/CLEAN/                 The installable `CLEAN` python package
    │   ├── model.py               VanillaNet / LayerNormNet / BatchNormNet / InstanceNorm
    │   │                          projection-head architectures (LayerNormNet is what's used)
    │   ├── losses.py              SupCon-Hard loss implementation
    │   ├── dataloader.py          Triplet/SupCon batch sampling incl. hard-negative mining
    │   ├── distance_map.py        Pairwise distance + EC-cluster-center computation
    │   ├── infer.py               infer_maxsep / infer_pvalue (see §2.2)
    │   ├── evaluate.py            Metric computation (precision/recall/F1/AUC), CSV writers
    │   ├── uncertainty.py         get_dist() — pos/neg distance sampling used to fit the GMM
    │   └── utils.py               ESM embedding retrieval, FASTA/CSV helpers, seeding, etc.
    ├── build.py                   `python build.py install` → installs src/CLEAN as a package
    ├── build/, __pycache__/, src/CLEAN.egg-info/   Build artifacts from that install step
    │
    ├── train-triplet.py           Train with triplet margin loss
    ├── train-supconH.py           Train with SupCon-Hard loss
    ├── gmm.py                     Fit the GMM confidence ensemble (§2.3)
    │
    ├── CLEAN_infer_fasta.py       Main "just give me predictions" entry point:
    │                              FASTA in data/inputs/ → results/inputs/<name>_maxsep.csv,
    │                              using the pretrained split100 model + GMM confidence
    ├── CLEAN_inference.py         Larger-scale variant: streams ESM-1b embedding on the fly
    │                              in batches (no need to pre-save all esm_data/*.pt files),
    │                              writes chunked results/inputs/<name>_<start>_<end>.csv
    ├── CLEAN_inference.sh         Example CLI invocation of CLEAN_inference.py
    ├── DEMO.ipynb                 Notebook walkthrough: embed → infer_pvalue + infer_maxsep
    │                              on the "new" (392-sequence) benchmark set
    ├── CLEAN_infer_on_the_fly.ipynb  Notebook version of CLEAN_inference.py's on-the-fly path
    │
    ├── score_audit.py             ★ CUSTOM — not part of upstream CLEAN (see §5)
    ├── AUDIT_NOTES.md             ★ CUSTOM — audit findings summary (see §5)
    │
    ├── data/
    │   ├── split10/30/50/70/100.csv   Training splits at different sequence-identity
    │   │                               clustering thresholds (Entry / EC number / Sequence, TSV)
    │   ├── new.csv, new.fasta          The "New-392" held-out benchmark (392 enzymes,
    │   │                               discovered after training data cutoff — used to
    │   │                               reproduce the paper's headline numbers)
    │   ├── ec_train_frequency.csv      ★ CUSTOM — per-EC training-set count + rarity tier
    │   │                               (5243 rows: common / rare / very_rare / singleton)
    │   ├── inputs/init.fasta, price.fasta   Small example FASTA inputs (init.fasta is what
    │   │                               the Docker entrypoint runs by default)
    │   ├── datasets/halogenase.csv, price.csv   Case-study datasets from the paper
    │   │                               (novel/hard enzyme families used to show CLEAN's edge
    │   │                               over BLAST-style annotation transfer)
    │   ├── esm_data/                   541 precomputed ESM-1b embeddings (*.pt, one per
    │   │                               UniProt accession) — avoids recomputing embeddings
    │   │                               for sequences already looked up
    │   └── pretrained/                 Downloaded model weights (already present locally):
    │       ├── split70.pth, split100.pth    CLEAN projection-head weights (trained on the
    │       │                                70%- and 100%-identity clustering splits)
    │       ├── 70.pt, 100.pt                Precomputed EC cluster-center embeddings for
    │       │                                each split (skips recomputing all 5000+ centers)
    │       └── gmm_ensumble.pkl             The GMM confidence ensemble (§2.3)
    │
    ├── results/                   Prediction output CSVs (Entry, predicted EC(s) + score)
    │   ├── new_maxsep.csv, new_pvalue.csv     Predictions for the New-392 benchmark
    │   └── inputs/price_maxsep.csv
    │
    └── audit_out/                 ★ CUSTOM — score_audit.py's output tables (see §5)
        ├── maxsep_by_depth.csv
        ├── maxsep_by_tier.csv
        ├── maxsep_tier_x_depth.csv
        └── maxsep_score_separation.csv
```

The `app/esm/` directory (referenced by the README) is where the upstream Meta `facebookresearch/esm`
repo gets `git clone`'d during install — it's the ESM-1b model code, not part of CLEAN itself.

---

## 4. What's already set up in this local copy

Unlike a fresh clone, this checkout has already done the expensive setup steps:

- **Pretrained weights are downloaded** (`app/data/pretrained/`: `split70.pth`, `split100.pth`,
  `70.pt`, `100.pt`, `gmm_ensumble.pkl`, plus the original `pretrained.zip`), so `infer_maxsep(...,
  pretrained=True)` works without hitting the Google Drive links in the README.
- **541 ESM-1b embeddings are cached** in `app/data/esm_data/*.pt`, covering (at least) the
  New-392 benchmark set, so `retrive_esm1b_embedding(...)` doesn't need to rerun for those IDs.
- **Predictions have already been generated** for the New-392 set under both algorithms
  (`results/new_maxsep.csv`, `results/new_pvalue.csv`) and for the `price` case-study set
  (`results/inputs/price_maxsep.csv`).
- A Python virtualenv exists at `app/.venv`.

`git status` shows the tracked-file diffs are limited to regenerated build metadata
(`src/CLEAN.egg-info/*`) and a re-run `price_maxsep.csv`; the meaningful new work
(`AUDIT_NOTES.md`, `score_audit.py`, `data/ec_train_frequency.csv`, `audit_out/`) is **untracked**
— i.e. not yet committed to git.

---

## 5. The audit layer (the actual work-in-progress here)

This is the part that isn't in the upstream repo. It reads as an independent verification +
stress-test of CLEAN's claims, going beyond the paper's own aggregate precision/recall/F1/AUC to
ask two sharper questions: *"does it actually get the exact reaction right, not just the
family?"* and *"is rare/novel-enzyme performance as good as the paper's overall number implies?"*

### 5.1 Reproduction check

`AUDIT_NOTES.md` confirms the published New-392 numbers reproduce **exactly** from the
pretrained split100 weights in this repo:

| Algorithm | Precision | Recall | F1 | AUC |
|---|---|---|---|---|
| max-separation | 0.596 | 0.479 | 0.497 | 0.739 |
| p-value (1e-5, 20k) | 0.558 | 0.477 | 0.482 | 0.737 |

392 test sequences, 177 distinct true EC numbers, 5242 EC cluster centers from the 241,025-pair
split100 training set. This is a sanity check that the local environment/weights match the paper
before trusting any further analysis.

### 5.2 `score_audit.py` — the finer-grained scoring tool

Given a ground-truth CSV, a CLEAN predictions CSV, and the EC-frequency index, it computes
**micro-averaged multi-label precision/recall/F1** sliced three ways:

1. **By EC digit depth (1–4)** — an EC number is a 4-part hierarchy (`class.subclass.sub-
   subclass.serial`). Truncating predictions/truth to `depth` digits before scoring answers "does
   the model get the enzyme's broad chemistry right even when it misses the exact reaction?"
2. **By rarity tier** — each test sample is tagged by its *rarest* true EC number's training-set
   frequency (`common` / `rare` / `very_rare` / `singleton` / `unseen`, from
   `ec_train_frequency.csv`), answering "does accuracy hold up on the EC numbers CLEAN was
   specifically designed to help with (few/no training examples)?"
3. **Both crossed together** — the "core audit table," rarity tier × depth.

It also reports **score separation**: are distance/confidence scores actually lower (more
confident) on correct calls than on wrong ones? Run it with:

```bash
cd app
python score_audit.py                                    # scores results/new_maxsep.csv
python score_audit.py --preds results/new_pvalue.csv --tag pvalue
```

Results land in `audit_out/{tag}_by_depth.csv`, `{tag}_by_tier.csv`, `{tag}_tier_x_depth.csv`,
`{tag}_score_separation.csv`, and are also printed as formatted tables to stdout.

### 5.3 Current findings (from `audit_out/`, max-separation, New-392)

**Depth breakdown** — accuracy degrades sharply as more digits (more chemical specificity) are
required:

| Depth | Meaning | Precision | Recall | F1 |
|---|---|---|---|---|
| 1 | broad class | 0.942 | 0.937 | 0.939 |
| 2 | sub-class | 0.901 | 0.897 | 0.899 |
| 3 | sub-sub-class | 0.846 | 0.865 | 0.855 |
| 4 | exact reaction | **0.518** | **0.479** | **0.498** |

CLEAN is very good at getting the general enzyme family right (F1 ≈ 0.94 at depth 1) but the
paper's headline F1 ≈ 0.50 is really an *exact-reaction* number — accuracy roughly halves going
from "right family" to "right specific reaction."

**Rarity breakdown** (depth 4, exact EC):

| Tier | n samples | Precision | Recall | F1 |
|---|---|---|---|---|
| common | 217 | 0.586 | 0.659 | 0.621 |
| rare | 79 | 0.451 | 0.383 | 0.414 |
| very_rare | 58 | 0.565 | 0.379 | 0.453 |
| singleton | 38 | 0.182 | 0.131 | **0.152** |

Performance on **singleton** EC numbers (only one example in the entire 241k-pair training set)
collapses to F1 ≈ 0.15 — far below the aggregate 0.497. This is a meaningful nuance: CLEAN's
whole pitch is handling low-data/orphan enzyme classes better than similarity search, but on the
hardest tier (truly one-shot) it's still struggling in absolute terms, even if it may still beat
BLAST-style baselines there (that comparison isn't computed in this repo).

**Score separation** — is the model's distance score a trustworthy confidence signal?

| Group | n | Mean distance |
|---|---|---|
| Correct predictions | 241 | 4.33 |
| Wrong predictions | 224 | 5.99 |

There *is* a directional signal (correct calls sit closer on average), but the gap is modest
relative to the spread, meaning distance alone is a fairly noisy confidence proxy at the
individual-prediction level — consistent with why the project also carries a GMM calibration
step (§2.3), though the audit doesn't yet directly test the calibrated GMM confidence values by
this same correct-vs-wrong split (`score_audit.py`'s docstring flags that the score's meaning
differs between raw max-sep distance and GMM-calibrated confidence — something to be careful
about if extending this analysis).

### 5.4 Open ends visible from the current state

- `audit_out/` currently only has results for `maxsep`; `results/new_pvalue.csv` exists but
  hasn't been run through `score_audit.py --tag pvalue` yet (no `pvalue_*.csv` files in
  `audit_out/`) — the depth/tier/score-separation comparison between the two inference algorithms
  is incomplete.
  - **Correction (checked 2026-08-10 in this pass):** `pvalue_*.csv` files do **not** exist yet in
    `audit_out/` — confirmed by directory listing, so this gap is real, not stale.
- The custom audit files (`score_audit.py`, `AUDIT_NOTES.md`, `data/ec_train_frequency.csv`,
  `audit_out/`) are untracked in git — worth committing if this analysis is meant to be kept /
  shared.
- `gmm.py` writes to `./gmm_test/`, a directory that doesn't currently exist in the repo listing —
  running it as-is will need that directory created first (or the script's output path checked).

---

## 6. Running things

### 6.1 One-time setup

```bash
cd app
pip install -r requirements.txt
git clone https://github.com/facebookresearch/esm.git
python build.py install
```

(In this checkout, `data/pretrained/` and much of `data/esm_data/` are already populated — see §4
— so the large downloads in the upstream README may not be needed again.)

### 6.2 Predict EC numbers for a new FASTA file

Drop a file in `app/data/inputs/`, then:

```bash
cd app
python CLEAN_infer_fasta.py --fasta_data <name-without-extension>
```

Output: `results/inputs/<name>_maxsep.csv` (with GMM-calibrated confidence scores, since
`CLEAN_infer_fasta.py` passes `gmm='./data/pretrained/gmm_ensumble.pkl'` under the hood).

For large batches without pre-saving every ESM embedding to disk, use `CLEAN_inference.py` /
`CLEAN_inference.sh` instead (streams ESM-1b in batches, GPU-aware via `--gpu_id`).

### 6.3 Re-run the audit

```bash
cd app
python score_audit.py                                   # maxsep
python score_audit.py --preds results/new_pvalue.csv --tag pvalue   # fill the pvalue gap noted in §5.4
```

### 6.4 Train a new model

```bash
cd app
python ./train-triplet.py --training_data split10 --model_name split10_triplet --epoch 2500
# or
python ./train-supconH.py --training_data split10 --model_name split10_supconH --epoch 1500 --n_pos 9 --n_neg 30 -T 0.1
```

See the upstream [`README.md`](README.md) §3 for the required one-time preprocessing
(`mutate_single_seq_ECs`, `compute_esm_distance`) before the first training run on a given split.

### 6.5 Docker (CPU)

```bash
docker pull moleculemaker/clean-image-amd64
```

Builds from `Dockerfile`, which bakes in weight download + `python build.py install`, and
`entrypoint.sh` runs `CLEAN_infer_fasta.py --fasta_data init` by default. Needs >12 GB container
memory per the upstream README.

---

## 7. Licensing & citation

The repo ships a **`NON-EXCLUSIVE RESEARCH USE LICENSE FOR CLEAN SOFTWARE.pdf`** — this is a
research-use license, not a standard permissive open-source license (MIT/Apache/BSD). Read it
before any commercial or production use.

If referencing this work:

```bibtex
@article{doi:10.1126/science.adf2465,
  author = {Tianhao Yu and Haiyang Cui and Jianan Canal Li and Yunan Luo and Guangde Jiang and Huimin Zhao},
  title = {Enzyme function prediction using contrastive learning},
  journal = {Science},
  volume = {379},
  number = {6639},
  pages = {1358-1363},
  year = {2023},
  doi = {10.1126/science.adf2465}
}
```

A hosted web server (default parameters only) is available at the
[MMLI AlphaSynthesis tools site](https://clean.frontend.mmli1.ncsa.illinois.edu/configuration).
