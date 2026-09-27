# cleanx RUNBOOK — exact commands, smoke tests, pass criteria

Companion to [`../../PUBLICATION_PLAN.md`](../../PUBLICATION_PLAN.md) (the *why*). This file is the *how*.

## 0. Open a terminal in the app folder

PowerShell **or** cmd:

```bat
cd "C:\Users\pdlpr\OneDrive\Desktop\chem project\CLEAN\app"
```

Always start commands with `.\` (e.g. `.\cx_smoke.bat`); that works in both shells, including
when Windows is set not to run programs from the current folder.

`cx.bat` is the launcher every script runs through: Windows Application Control on this PC
blocks `.venv\Scripts\python.exe`, so `cx.bat` runs the same Python 3.10 base interpreter the
venv was built from, with the venv's packages. Verified environment (2026-09-27): Python 3.10.20,
numpy 1.22.3, pandas 1.4.2, scipy 1.7.3, scikit-learn 1.2.0, torch 1.11.0+cpu, fair-esm 2.0.0,
biopython 1.83 (added for the identity step with `--no-deps`; nothing else changed).

## 1. The three commands you need

### Command 1 — smoke test of everything (≈2 min)

```bat
.\cx_smoke.bat
```

**PASS** = the last line is `ALL 7 SMOKE TESTS PASSED`. It stops at the first failure and prints
`SMOKE TEST FAILED`. Expected numbers along the way (verified 2026-09-27):
- step 1: `CX00 SUMMARY: 9/9 checks passed` (exact reproduction of CLEAN: P 0.5955, R 0.4791,
  F1 0.4967; 0/392 prediction mismatches)
- step 5: `exemplar - centroid micro_F1: +0.035 [+0.003,+0.070]`
- the line `WARNING: new392 is also evaluated; its balanced* results are optimistic` is expected in
  the smoke run (λ is tuned on New-392 there only to exercise the code)

### Command 2 — ESM-1b equivalence test (≈5 min; needs ~8–10 GB free RAM)

Close browsers/IDEs first (Task Manager → Memory). Then:

```bat
.\cx.bat cx_03_embed.py --verify-only
```

**PASS** = a line `[PASS] recomputed ESM-1b embeddings match CLEAN's stored ones for 3 New-392
proteins` with max |diff| < 1e-3 and min cosine > 0.99999. The next line prints the CPU speed in
**s/protein**; the full embedding run takes about 4,099 × that. **If it says FAIL, stop and send me
the output**: the rest of the pipeline assumes identical embeddings.

### Command 3 — full pipeline (≈2–4 h, almost all of it the embedding step; resumable)

```bat
.\cx_full.bat
```

It runs, in order: embeddings (repeats the equivalence test first) → E2 split100 (primary endpoint)
→ E2 split70 (robustness) → E3 → E4 → E1 incl. the temporal test → figures. **PASS** = last line
`FULL PIPELINE FINISHED`. If the PC sleeps or the run stops, run `.\cx_full.bat` again: finished
embeddings are skipped, and every analysis step is fast.

Then send me these files (or paste their contents) and I will fill the numbers into the paper:
`results\cleanx\E2_split100\E2_summary.md`, `results\cleanx\E2_split100\E2_deltas.csv`,
`results\cleanx\E3_split100\E3_summary.md`, `results\cleanx\E4\E4_summary.md`,
`results\cleanx\E4\E4_min_certifiable_risk.csv`, `results\cleanx\E4\E4_openset.csv`.

## 2. What has already been run and verified (2026-09-27)

| Step | Script | Status | Output |
|---|---|---|---|
| E0 reproduction checks | `cx_00_check.py` | ✅ 9/9 PASS | console |
| E1 evaluation audit (New-392, Price-149) | `cx_01_audit.py` | ✅ | `results/cleanx/E1/` |
| Temporal benchmark build | `cx_02_build_temporal.py --with-pfam` | ✅ UniProt 2026_03, ENZYME 02-Sep-2026, Pfam 38.2 | `data/cleanx/sets/`, `downloads/` |
| Identity + homology baseline | `cx_04_identity.py` | ✅ (Python fallback) | `data/cleanx/identity/` |
| E2 geometry on New-392, Price-149 | `cx_05_geometry.py --datasets new392 price149` | ✅ | `results/cleanx/E2_split100/` |
| E3 NISE biology (legacy link) | `cx_06_ec_biology.py` | ✅ | `results/cleanx/E3_split100/` |
| Figures from legacy results | `cx_08_figures.py` | ✅ checked visually | `results/cleanx/figures/` |
| Smoke suite | `cx_smoke.bat` | ✅ ALL 7 PASSED | — |
| Temporal-only code paths (unseen ECs, EC-transfer mapping, λ tuning on CAL, non-enzyme open set, reliability figures) | integration test with synthetic vectors in temp folders | ✅ 5/5 PASS | — |
| Embedding orchestration (PASS/FAIL branches), MMseqs2 parser, E3b | unit tests | ✅ | — |
| **ESM-1b load itself** | `cx_03_embed.py --verify-only` | ⏳ **your Command 2** (not run here: only ~4 GB RAM free) | — |

## 3. Individual steps (only if you want to run pieces by hand)

```bat
.\cx.bat cx_00_check.py
.\cx.bat cx_03_embed.py --limit 5
.\cx.bat cx_03_embed.py
.\cx.bat cx_05_geometry.py
.\cx.bat cx_05_geometry.py --backbone split70
.\cx.bat cx_06_ec_biology.py
.\cx.bat cx_07_reliability.py
.\cx.bat cx_01_audit.py --datasets new392 price149 sp_temporal_test
.\cx.bat cx_08_figures.py
```

Optional E3b (mechanism: are non-homologous subgroups separate clusters in CLEAN space? ≈1.5–2 h):

```bat
.\cx.bat cx_03_embed.py --sets nise_members
.\cx.bat cx_06_ec_biology.py --subgroups
```

Do **not** run `cx_02_build_temporal.py --force`: it re-downloads UniProt, and a newer release
changes the benchmark. Without `--force`, cached downloads are reused.

## 4. Optional: MMseqs2 identity (recommended for final paper numbers)

Identity now comes from the built-in fallback (3-mer prefilter + BLOSUM62 Smith–Waterman on the
top 20 candidates; identity = identical residues / shorter length). To confirm with MMseqs2 on
Linux or WSL:

```bat
.\cx.bat cx_04_identity.py --write-fasta
```

Then in a Linux shell inside `app/data/cleanx/fasta`:

```bash
for s in new392 price149 sp_temporal_cal sp_temporal_test; do mmseqs easy-search $s.fasta split100.fasta $s.m8 tmp_mm -s 7.5 --max-seqs 300 -e 10 --format-output query,target,fident,alnlen,qlen,tlen,bits,evalue; done
```

Copy the four `.m8` files to `app\data\cleanx\identity\mmseqs\`, then:

```bat
.\cx.bat cx_04_identity.py --tool m8 --sets new392 price149 sp_temporal_cal sp_temporal_test
```

This overwrites the identity tables; re-run `.\cx.bat cx_05_geometry.py` afterwards so the
identity-stratified results use them.

## 5. Output files

| File | Contents |
|---|---|
| `results/cleanx/E1/E1_blind_spot.csv` | predictions, TP, FP, FPs invisible to the reported precision |
| `results/cleanx/E1/E1_overprediction.csv` | reported vs micro scores when 1–10 ECs are predicted per protein |
| `results/cleanx/E1/E1_exact_duplicates.csv`, `E1_label_drift.csv` | leakage and label drift |
| `results/cleanx/E2_*/E2_deltas.csv` | paired-bootstrap differences vs CLEAN centroid (**primary endpoint**) |
| `results/cleanx/E2_*/E2_tuning.csv/.json` | λ chosen on SP-temporal CAL |
| `results/cleanx/E2_*/E2_per_query.csv` | per-protein predictions and F1 for every method |
| `results/cleanx/E3_*/E3_stats.json` | NISE statistics, size-stratified permutation tests, test-set link |
| `results/cleanx/E4/E4_calibration.csv` | ECE/Brier by predicted and true tier |
| `results/cleanx/E4/E4_risk_control.csv`, `E4_min_certifiable_risk.csv` | certified thresholds; smallest certifiable risk per tier |
| `results/cleanx/E4/E4_openset.csv` | abstention on unseen-EC enzymes and non-enzymes; AUROCs |
| `results/cleanx/figures/*.png/.pdf/.csv` | paper figures and the numbers behind each |
| `data/cleanx/sets/sp_temporal_manifest.json` | how the benchmark was built (filters, counts, releases) |

## 6. Troubleshooting

- `An Application Control policy has blocked this file` → you ran the venv's python directly; use
  `.\cx.bat`.
- `'cx.bat' is not recognized` → you typed `cx.bat` without `.\`.
- `FileNotFoundError: … ESM-1b embeddings missing` → the embedding step is not finished for that set.
- ESM load very slow / PC freezes → not enough free RAM; close apps and re-run (resumable).
- Git: keep `data/cleanx/`, `results/` and CLEAN's own files out of any **public** repository
  (PUBLICATION_PLAN §13, CLEAN licence).
