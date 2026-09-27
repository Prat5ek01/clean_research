# cleanx RUNBOOK — exact commands, smoke tests, pass criteria

Companion to [`../../PUBLICATION_PLAN.md`](../../PUBLICATION_PLAN.md) (the *why*). This file is the *how*.
Every step has a **smoke test with pass/fail criteria**; run the full step only after the smoke passes.

## 0. Where and how to run

Open **PowerShell or cmd** and go to the app folder:

```bat
cd "C:\Users\pdlpr\OneDrive\Desktop\chem project\CLEAN\app"
```

Run every script through the launcher **`.\cx.bat <script> [args]`**. Why: Windows Application
Control on this PC blocks `.venv\Scripts\python.exe`; `cx.bat` runs the same Python 3.10 base
interpreter the venv was built from, with the venv's packages on `PYTHONPATH`. It works from both
PowerShell and cmd.

Verified environment (2026-09-27): Python 3.10.20, numpy 1.22.3, pandas 1.4.2, scipy 1.7.3,
scikit-learn 1.2.0, torch 1.11.0+cpu, fair-esm 2.0.0, biopython 1.83 (added for the identity step,
installed with `--no-deps` so nothing else changed).

Rules that avoid wasted runs:
- One heavy job at a time. Another Claude session (QCE project) also runs Python jobs on this PC;
  check Task Manager before long runs. Running jobs in parallel has also triggered Application
  Control DLL blocks on scipy before.
- The ESM-1b step needs **~8–10 GB of free RAM while the model loads**. Close browsers/IDEs first.
- Do not re-run `cx_02_build_temporal.py --force` unless you intend to: it re-downloads UniProt,
  and a newer UniProt release changes the benchmark. The current build is UniProt **2026_03**,
  ENZYME **02-Sep-2026**, Pfam **38.2**.

## 1. Status of each step (as of 2026-09-27)

| Step | Script | Status | Output |
|---|---|---|---|
| E0 checks | `cx_00_check.py` | ✅ done, 9/9 PASS | console |
| E1 evaluation audit (legacy) | `cx_01_audit.py` | ✅ done | `results/cleanx/E1/` |
| Temporal benchmark build | `cx_02_build_temporal.py --with-pfam` | ✅ done | `data/cleanx/sets/`, `data/cleanx/downloads/` |
| Identity + homology baseline | `cx_04_identity.py` | ✅ done (Python fallback) | `data/cleanx/identity/`, `results/cleanx/preds/*_homology_tophit.csv` |
| ESM-1b embeddings, temporal sets | `cx_03_embed.py` | ⏳ **your step** (~2–4 h CPU) | `data/esm_data/*.pt` |
| E2 inference geometry | `cx_05_geometry.py` | legacy ✅; temporal ⏳ needs embeddings | `results/cleanx/E2_split100/` |
| E3 NISE biology | `cx_06_ec_biology.py` | run after E2 | `results/cleanx/E3_split100/` |
| E3b (optional) | `cx_03_embed.py --sets nise_members` then `cx_06 --subgroups` | optional, ~1.5–2 h CPU | `E3b_nise_subgroups.csv` |
| E4 reliability | `cx_07_reliability.py` | ⏳ needs embeddings | `results/cleanx/E4/` |
| Figures | `cx_08_figures.py` | after E2–E4 | `results/cleanx/figures/` |

## 2. Step-by-step

### Step A — environment check (1 min)

```bat
.\cx.bat cx_00_check.py
```

**PASS:** last line `CX00 SUMMARY: 9/9 checks passed`. Every check must say `[PASS]`, including
"New-392 max-separation replica" (0 mismatches) and "CLEAN published metric reproduced"
(P 0.5955 R 0.4791 F1 0.4967). **FAIL on anything → stop**; the rest depends on exact reproduction.

### Step B — ESM-1b embeddings for the temporal sets (the only heavy step)

**B1 smoke / equivalence test (~5 min incl. model load):**

```bat
.\cx.bat cx_03_embed.py --verify-only
```

This re-embeds 3 New-392 proteins (shortest, median, longest = 1022 aa) and compares them with the
embeddings CLEAN already used.

**PASS:** a line starting `[PASS] recomputed ESM-1b embeddings match` with max |diff| < 1e-3 and
min cosine > 0.99999. The next line prints the CPU speed in **s/protein**; multiply by 4,099 for the
full-run time. **FAIL → do not continue** (it would mean the embedding pipeline differs from CLEAN's).

**B2 tiny run (~5 min):** 5 proteins from each set, proves writing/resume on real data.

```bat
.\cx.bat cx_03_embed.py --limit 5
```

**PASS:** ends with `CX03 DONE: 0 embeddings missing`.

**B3 full run (resumable; 839 + 1,260 + 2,000 = 4,099 proteins):**

```bat
.\cx.bat cx_03_embed.py
```

Prints `ETA` every 10 batches. If it stops (sleep, reboot), run the same command again; finished
proteins are skipped. **PASS:** `CX03 DONE: 0 embeddings missing`.

### Step C — E2 inference geometry (≈10–20 min)

Smoke (1 min):

```bat
.\cx.bat cx_05_geometry.py --smoke
```

**PASS:** prints New-392 `centroid micro-F1 0.498` and `exemplar micro-F1 0.532`, delta
`+0.035` with a 95% CI above 0.

Full (primary + legacy sets), then the second backbone as robustness check:

```bat
.\cx.bat cx_05_geometry.py
.\cx.bat cx_05_geometry.py --backbone split70
```

Primary endpoint = `results/cleanx/E2_split100/E2_deltas.csv`, row
`dataset=sp_temporal_test, subset=closed, comparison=exemplar - centroid, metric=micro_F1`.

### Step D — E3 biology of centroid failures (≈2–5 min)

```bat
.\cx.bat cx_06_ec_biology.py
```

**PASS:** summary reports ~183 NISE-like ECs (≈7% of analysed ECs) and, because E2 ran first, a
line per dataset "exemplar-minus-centroid per-protein F1 gain, true EC NISE-like … vs other …".

Optional E3b (mechanism: are non-homologous subgroups separate clusters in CLEAN space?):

```bat
.\cx.bat cx_03_embed.py --sets nise_members
.\cx.bat cx_06_ec_biology.py --subgroups
```

(The FASTA `data/cleanx/fasta/nise_members.fasta`, 3,817 proteins, was already written.)

### Step E — E4 reliability (≈5 min)

Smoke (New-392 split in halves; code path only):

```bat
.\cx.bat cx_07_reliability.py --smoke
```

**PASS:** writes `results/cleanx/E4_smoke/E4_summary.md`; ECE for `clean_gmm(raw)` on
very_rare/singleton predictions ≈ 0.49 vs ≈ 0.18 for common.

Full (calibrate on SP-temporal CAL; test on SP-temporal TEST, Price-149, New-392; non-enzymes):

```bat
.\cx.bat cx_07_reliability.py
```

### Step F — E1 on the temporal test set + figures

```bat
.\cx.bat cx_01_audit.py --datasets new392 price149 sp_temporal_test
.\cx.bat cx_08_figures.py
```

Figures land in `results/cleanx/figures/` as PNG (300 dpi), PDF, and a CSV of the plotted numbers.

## 3. Optional: MMseqs2 identity (recommended for the final paper numbers)

The identity/homology numbers now come from the built-in Python fallback (3-mer prefilter + BLOSUM62
Smith–Waterman on the top 20 candidates; identity = identical residues / shorter length). Reviewers
know MMseqs2, so re-run with it on any Linux machine or WSL:

```bat
.\cx.bat cx_04_identity.py --write-fasta
```

Then on Linux/WSL, from the folder containing `data/cleanx/fasta/`:

```bash
for s in new392 price149 sp_temporal_cal sp_temporal_test; do mmseqs easy-search fasta/$s.fasta fasta/split100.fasta $s.m8 tmp_mm -s 7.5 --max-seqs 300 -e 10 --format-output query,target,fident,alnlen,qlen,tlen,bits,evalue; done
```

Copy the `.m8` files into `data\cleanx\identity\mmseqs\` and parse them:

```bat
.\cx.bat cx_04_identity.py --tool m8 --sets new392 price149 sp_temporal_cal sp_temporal_test
```

Report the agreement between the two methods' identity bins in the Supplement.

## 4. What each output file is

| File | Contents |
|---|---|
| `results/cleanx/E1/E1_blind_spot.csv` | predictions, TP, FP, FPs invisible to the reported precision, reported vs micro P |
| `results/cleanx/E1/E1_metrics_long.csv` | every scheme × EC depth × tier |
| `results/cleanx/E1/E1_exact_duplicates.csv` | test proteins identical to a training sequence |
| `results/cleanx/E1/E1_label_drift.csv` | legacy labels vs current Swiss-Prot |
| `results/cleanx/E2_*/E2_deltas.csv` | paired-bootstrap differences vs CLEAN centroid (primary endpoint) |
| `results/cleanx/E2_*/E2_per_query.csv` | per-protein predictions and F1 for every method |
| `results/cleanx/E3_*/E3_stats.json` | NISE statistics, permutation test, test-set link |
| `results/cleanx/E4/E4_calibration.csv` | ECE/Brier by predicted and true tier |
| `results/cleanx/E4/E4_risk_control.csv` | certified thresholds: realized risk and coverage per tier |
| `results/cleanx/E4/E4_openset.csv` | abstention on unseen-EC enzymes and non-enzymes; AUROCs |
| `data/cleanx/sets/sp_temporal_manifest.json` | how the benchmark was built (filters, counts, releases) |

## 5. Troubleshooting

- `Program 'python.exe' failed to run: An Application Control policy has blocked this file` →
  you ran the venv's python directly; use `.\cx.bat`.
- `FileNotFoundError: … ESM-1b embeddings missing` → Step B not finished for that set.
- ESM load is extremely slow / PC freezes → not enough free RAM; close apps, run again (resumable).
- `enzyme.dat missing` → run `.\cx.bat cx_02_build_temporal.py` (without `--force` it reuses cached
  downloads).
- Git: `data/cleanx/` and `results/` are local artefacts; keep the GitHub repository **private**
  (see PUBLICATION_PLAN §13, CLEAN licence).
