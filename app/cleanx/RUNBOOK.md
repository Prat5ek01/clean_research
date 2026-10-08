# cleanx RUNBOOK — exact commands, smoke tests, pass criteria

Companion to [`../../PUBLICATION_PLAN.md`](../../PUBLICATION_PLAN.md) (the *why*). This file is the *how*.
Last updated 2026-09-28, after the confirmatory run: audit fixes and the BLASTp homology baseline (§1.1).

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
biopython 1.83 (added for the identity step with `--no-deps`; nothing else changed). No new packages since.

## 1. Where things stand

| Item | State |
|---|---|
| Pre-registration commit | `c8e42f704a406104a6fb61868392e73991f25698` (E2 endpoints fixed 2026-09-27; E5, E6 and the E2 sensitivity analyses 2026-09-28) — cite this hash |
| Confirmatory run of that commit (your PC, 2026-09-28) | `ALL 10 SMOKE TESTS PASSED`; `FULL PIPELINE FINISHED`; results in PUBLICATION_PLAN §2b |
| ESM-1b equivalence | PASS: max \|diff\| 1.76e-07, min cosine 1.0000000 (3 New-392 proteins) |
| Cross-platform check | E5 re-run in a Linux copy reproduces your PC's primary endpoint, rank agreement and ρ exactly |
| Final run with BLASTp (your PC, 2026-09-28) | smoke 11/11; `FULL PIPELINE FINISHED`; every pre-registered number identical to the confirmatory run; results commit `532bc1f` |
| After pulling the figure-caption fix | run only `.\cx.bat cx_08_figures.py` (≈30 s; redraws all figures, including the pilot figures figS2/figS3, from the committed results), then `git add results/cleanx/figures`, commit, push |

### 1.1 Changes after the confirmatory run (none touches a pre-registered endpoint)

1. **BLASTp homology baseline and identity** (`cx_11_blast.py`, `cx_04_identity.py --tool blast`). The old baseline
   aligned only the 30 best 3-mer candidates per query; CLEAN's paper and EC-Bench use BLASTp. This section
   (then §4) already planned to replace it with a standard aligner for the paper's numbers before the run.
   BLAST+ may be blocked on this PC, so it was run in a Linux copy; the hit tables are committed with SHA-256
   fingerprints of the exact sequences searched, and every run on this PC verifies them first.
2. **Bug fix:** in the EC-cluster bootstrap, per-protein F1 was divided by N instead of the resample's total
   weight (only the `samples_F1` rows of the `closed|cluster_bootstrap_EC` sensitivity analysis were affected;
   micro and union-macro F1 are ratios and were correct).
3. **Reporting:** every realized selective risk in E4 now carries its error count and an exact 95% CI; E6
   bootstrap bounds that fall between two non-certifiable resamples are +∞ (numpy returned NaN there), and E6
   adds `years_at_deficit_ci_lo` (years of Swiss-Prot accrual at the deficit's lower 95% bound).
4. **Figures** (`cx_08_figures.py`, reporting only): 95% CIs on every F1 (paired bootstrap, B = 2000, as in E2)
   and every risk (exact); the reliability diagrams compare the released GMM confidence and the tier-conditional
   calibration on the *same* CLEAN predictions; risk–coverage and per-tier risk control show all three
   geometries; the EC-Bench figure marks the significant reversal; the budget figure shows the deficit CI;
   new `fig10_primary_endpoint` (forest plot of the primary endpoint, its sensitivity analyses, the other
   methods and the other test sets). F1 values are recomputed from `E2_per_query.csv` and a figure is skipped if
   they differ from `E2_metrics_long.csv`.

## 2. The commands you need now, in order

### Command A — save the confirmatory results of c8e42f7 ✅ done (your commit 49a6747, pushed)

Kept here for the record; do not run it again.

```bat
git add results/cleanx data/cleanx/identity
git status
git commit -m "Confirmatory run of c8e42f7 (2026-09-28): E1-E6 results and figures"
git push
```

**PASS** = `git status` shows nothing staged under `data/esm_data`, and `git push` ends without an error.
This keeps the pre-registered run in the history exactly as it came out.

### Command B — get the audited code

```bat
git pull --no-rebase --no-edit
git log -1 --format="%H %s"
```

**PASS** = no conflict message; the last commit is the audit commit (I will tell you its hash).

### Command C — smoke test (≈3 min)

```bat
.\cx_smoke.bat
```

**PASS** = the last line is `ALL 11 SMOKE TESTS PASSED`. It stops at the first failure and prints
`SMOKE TEST FAILED`. Expected lines (times from your PC, 2026-09-28):

| Step | What | Expected line(s) | Time |
|---|---|---|---|
| 1 | environment + exact CLEAN reproduction | `CX00 SUMMARY: 9/9 checks passed` | 10 s |
| 2 | E1 audit (smoke) | `E1 DONE` | 4 s |
| 3 | temporal benchmark builder (cached downloads) | `CX02 DONE` | 3 s |
| 4 | identity search, 5 proteins (Python fallback code path) | `homology top-hit exact-EC agreement: 5/5`, `CX04 DONE` | 12 s |
| 5 | E2 geometry (smoke) | `exemplar - centroid micro_F1: +0.034 [+0.003,+0.070]`, `E2 DONE` | 13 s |
| 6 | E3 NISE biology (smoke) | `E3 DONE` | 22 s |
| 7 | E4 reliability (smoke) | `E4 DONE` | 6 s |
| 8 | unit + integration tests (synthetic data, temporary folders only) | `ALL 29 TESTS PASSED` | ≈1 min |
| 9 | E5 evaluation card (smoke) | `[PASS] EC-Bench reported weighted F1 reproduced (32 values ...) -- max \|diff\| 0.0e+00`, `[PASS] Propositions 1-2 and Lemma 3 hold ... N/N tables`, `CX09 SMOKE DONE` | 1 s |
| 10 | E6 certification budget (smoke) | three `[PASS]` lines, `CX10 SMOKE DONE` | 51 s |
| 11 | BLASTp tables | `[PASS]` for the training sequences, each of the four sets, each stored table and the TEST-vs-CAL pair; `SKIP re-search: blastp/makeblastdb not installed here`; `CX11 SMOKE DONE` | 5 s |

The line `WARNING: new392 is also evaluated; its balanced* results are optimistic` is expected in step 5 (λ is
tuned on New-392 there only to exercise the code).

### Command D — full pipeline (≈15 min: every embedding already exists and is skipped)

```bat
.\cx_full.bat
```

**PASS** = the last line is `FULL PIPELINE FINISHED`. The first header must read `1/12`.

| Step | Command inside `cx_full.bat` | PASS line | Time (your PC) |
|---|---|---|---|
| 1 | `cx_03_embed.py` (repeats the equivalence test) | `[PASS] recomputed ESM-1b embeddings match ...`, `CX03 DONE: 0 embeddings missing` | ≈1.5 min now; ≈6 h the first time |
| 2 | `cx_11_blast.py --verify` | only `[PASS]` lines, `CX11 VERIFY DONE` | 5 s |
| 3 | `cx_04_identity.py --tool blast --sets new392 price149 sp_temporal_cal sp_temporal_test` | an identity-bin table per set, `CX04 DONE` | seconds |
| 4 | `cx_04_identity.py --tool blast --pair sp_temporal_test:sp_temporal_cal` | `CX04 PAIR DONE` | seconds |
| 5 | `cx_05_geometry.py` (split100: **E2 primary endpoint**) | `E2 DONE` | ≈3 min |
| 6 | `cx_05_geometry.py --backbone split70` | `E2 DONE` | ≈2 min |
| 7 | `cx_06_ec_biology.py` | `E3 DONE` | 30 s |
| 8 | `cx_07_reliability.py` | `E4 DONE` | 50 s |
| 9 | `cx_01_audit.py --datasets new392 price149 sp_temporal_test` | `E1 DONE` | 10 s |
| 10 | `cx_09_evalcard.py` (**E5**) | `[PASS] EC-Bench reported weighted F1 reproduced (32 values ...)`, `[PASS] Propositions 1-2 and Lemma 3 hold ... N/N tables`, `E5 DONE` | 45 s |
| 11 | `cx_10_certbudget.py` (**E6**) | three `[PASS]` known-answer lines, `E6 DONE` | 3 min |
| 12 | `cx_08_figures.py` | 11 `wrote ...` lines from `fig1_blind_spot` to `fig9_certification_budget` (including the new `fig10_primary_endpoint`), no line containing `does not reproduce` or `WARNING`, `CX08 DONE` | 30 s |

The pre-registered primary endpoints must come out exactly as in the confirmatory run: E2 `balanced3 -
centroid macro_union_F1: +0.043 [+0.030,+0.060]` on `sp_temporal_test`; E5 `1 significantly reversed pairs`,
Kendall tau-b 0.900 [0.833, 0.933], Spearman rho 0.75; E6 (in `E6_summary.md`) very_rare and singleton at
alpha=0.2 `deficit CI [528, inf]` and `deficit CI [515, inf]`, rare `[CI 428, inf]`; E4 (in `E4_summary.md`) singleton risks
`(4/14, ...)`, `(5/15, ...)`, `(11/24, ...)` for centroid, exemplar, balanced3. None of these depends on the
BLASTp change. If any differs, stop and send me the output. The homology rows, the identity-stratified rows and
the CAL-disjoint subset change (BLASTp replaces the fallback aligner); that is expected.

If step 2 or 3 prints `FAIL` or `BLAST table mismatch`, the stored BLAST results do not belong to your local
sequences: stop and send me the output. If step 10 prints `E5 STOPPED: reproduction gate failed`, or any
script raises a SHA-256 error, stop likewise.

### Command E — save the final results

```bat
git add results/cleanx data/cleanx/identity
git status
git commit -m "Final run: BLASTp baseline and audit fixes"
git push
```

Then tell me; I check every table and figure and fill the numbers into PUBLICATION_PLAN.

### Command F — E7 label provenance (needs internet for UniProt; the analysis code must be committed)

The E7 pre-registration is the docstring of `cx_13_provenance.py` (drafted 2026-09-29, finalised 2026-10-07,
frozen by the commit that added it). The confirmatory run refuses to start if that code differs from the commit.

```bat
git pull
.\cx.bat cx_00_check.py
.\cx_smoke.bat
.\cx.bat tests\run_tests.py
.\cx.bat cx_12_uniprot_evidence.py --limit 5
.\cx.bat cx_12_uniprot_evidence.py
.\cx.bat cx_13_provenance.py --pilot
.\cx.bat cx_13_provenance.py
git add data/cleanx/provenance results/cleanx/E7 results/cleanx/E7_pilot results/cleanx/E7_smoke
git commit -m "E7 confirmatory run: UniProt evidence and label-provenance results"
git push
```

Expected: CX00 9/9, `ALL 13 SMOKE TESTS PASSED`, `ALL 38 TESTS PASSED`; `--limit 5` is a live check (it writes
`*_limit5` tables the analysis never reads); the full download covers New-392, SP-temporal CAL/TEST and the
EC-Bench test set (about 2,900 accessions, cached under `data/cleanx/downloads/uniprot_json/`). Send me the pilot
summary (`results/cleanx/E7_pilot/E7_summary.md`) and the confirmatory output.

## 3. What has been run and verified

| Step | Script | Status | Output |
|---|---|---|---|
| E0 reproduction checks | `cx_00_check.py` | ✅ 9/9 (your PC, 2026-09-27 and 2026-09-28) | console |
| Temporal benchmark build | `cx_02_build_temporal.py --with-pfam` | ✅ UniProt 2026_03, ENZYME 02-Sep-2026, Pfam 38.2 | `data/cleanx/sets/`, `downloads/` |
| ESM-1b embeddings (4,099 proteins) | `cx_03_embed.py` | ✅ equivalence PASS; 0 missing (your PC) | `data/esm_data/` (not in git) |
| Confirmatory run E1–E6 | `cx_full.bat` at c8e42f7 | ✅ (your PC, 2026-09-28) | `results/cleanx/` |
| E5 EC-Bench part re-run | same functions, Linux copy | ✅ identical primary endpoint, τ_b and ρ | — |
| BLASTp search | `cx_11_blast.py` | ✅ Linux copy, BLAST+ 2.12.0, ≈25 min on 4 threads | `data/cleanx/identity/blast/` |
| Unit + integration tests (29) | `tests\run_tests.py` | ✅ ALL 29 PASSED (Linux copy) | temporary folders only |
| Revised figures | `cx_08_figures.py` on your confirmatory results | ✅ Linux copy: every F1 reproduces `E2_metrics_long.csv`; all figures ≤ 7.5 in wide | `results/cleanx/figures/` (your run writes them) |
| **Final run with BLASTp** | `cx_full.bat` | ✅ your PC, commit 532bc1f (audited: only the expected rows changed) | `results/cleanx/` |

## 4. Individual steps (only if you want to run pieces by hand)

```bat
.\cx.bat cx_00_check.py
.\cx.bat cx_03_embed.py
.\cx.bat cx_11_blast.py --verify
.\cx.bat cx_04_identity.py --tool blast --sets new392 price149 sp_temporal_cal sp_temporal_test
.\cx.bat cx_04_identity.py --tool blast --pair sp_temporal_test:sp_temporal_cal
.\cx.bat cx_05_geometry.py
.\cx.bat cx_05_geometry.py --backbone split70
.\cx.bat cx_06_ec_biology.py
.\cx.bat cx_07_reliability.py
.\cx.bat cx_01_audit.py --datasets new392 price149 sp_temporal_test
.\cx.bat cx_09_evalcard.py
.\cx.bat cx_10_certbudget.py
.\cx.bat cx_08_figures.py
```

| Command | What it does | ≈ time |
|---|---|---|
| `.\cx.bat tests\run_tests.py` | all 29 tests (unit + integration on a synthetic app in a temporary folder) | 1 min |
| `.\cx.bat tests\run_tests.py --unit-only` | the 28 unit tests only | 2 s |
| `.\cx.bat cx_11_blast.py --smoke` | BLAST table verification, plus a re-search of 5 New-392 queries where BLAST+ is installed | 5 s |
| `.\cx.bat cx_09_evalcard.py --pilot` | EC-Bench Price-149 only (disclosed pilot) → `E5_pilot` | 30 s |
| `.\cx.bat cx_09_evalcard.py --skip-ecbench` | cards for our own predictions only (legacy + temporal) → `E5_local` | 5 s |
| `.\cx.bat cx_09_evalcard.py --card data\new.csv results\cleanx\preds\new392_split100_centroid_maxsep.csv` | evaluation card for any prediction file → `E5_card\` | 1 s |
| `.\cx.bat cx_10_certbudget.py --pilot` | E6 with the E4 smoke features as calibration data → `E6_pilot` | 2 min |
| `.\cx.bat cx_05_geometry.py --cap-per-ec 3` | exploratory variant of the capped sensitivity subset (pre-registered value: 5) | 3 min |

The `--card` example must print `reported (support-weighted) P/R/F1 0.596/0.479/0.497 | micro 0.518/0.479/0.498 |
union-macro F1 0.371 | per-protein F1 0.535 | 150/224 FPs invisible`.

Optional E3b (are non-homologous subgroups separate clusters in CLEAN space? ≈1.5–2 h):

```bat
.\cx.bat cx_03_embed.py --sets nise_members
.\cx.bat cx_06_ec_biology.py --subgroups
```

Do **not** run `cx_02_build_temporal.py --force`: it re-downloads UniProt, and a newer release changes the
benchmark. `cx_09_evalcard.py --force-download` is only for a corrupted EC-Bench file (pinned to commit
`5e98a2c` and fixed SHA-256 values).

## 5. The BLASTp baseline (how the committed tables were made)

`cx_11_blast.py`, run where NCBI BLAST+ is installed (here: a Linux copy, BLAST+ 2.12.0 from Ubuntu):
`makeblastdb -dbtype prot` on the split100 training sequences (no `-parse_seqids`, so accessions stay
verbatim), then per set `blastp -outfmt "6 qseqid sseqid pident length qlen slen bitscore evalue nident"
-evalue 10 -max_target_seqs 500` with default scoring (BLOSUM62, gap 11/1, composition-based statistics).
Per query the 10 best hits by bit score and the maximum-identity hit are kept
(`data/cleanx/identity/blast/<set>.tsv.gz`); `manifest.json` records the version, the exact commands, and
SHA-256 digests of every sequence set and table. Baseline = EC(s) of the best hit by bit score (ties: lower
E-value, then accession); identity = identical residues / length of the shorter sequence, maximised over hits.
To repeat it on a machine with BLAST+: `python cx_11_blast.py` (≈25 min on 4 threads), then commit
`data/cleanx/identity/blast/`.

## 6. Output files

| File | Contents |
|---|---|
| `results/cleanx/E1/E1_blind_spot.csv` | predictions, TP, FP, FPs invisible to the reported precision |
| `results/cleanx/E1/E1_overprediction.csv` | reported vs micro scores when 1–10 ECs are predicted per protein |
| `results/cleanx/E1/E1_exact_duplicates.csv`, `E1_label_drift.csv` | leakage and label drift |
| `results/cleanx/E2_*/E2_deltas.csv` | paired-bootstrap differences vs CLEAN centroid (**primary endpoint**: subset `closed`, `balanced3 - centroid`, `macro_union_F1`, `sp_temporal_test`); report `delta_obs` with `ci_lo`/`ci_hi`; tier rows carry `p_holm`; sensitivity rows `closed\|cluster_bootstrap_EC`, `closed_capped5`, `closed_cal_disjoint`, `closed_PE1` |
| `results/cleanx/E2_*/E2_bootstrap_ci.csv`, `E2_metrics_long.csv`, `E2_per_query.csv`, `E2_tuning.*` | per-method CIs, all metric families by subset, per-protein predictions, λ chosen on CAL |
| `results/cleanx/E3_*/E3_stats.json` | NISE statistics, size-stratified permutation tests (incl. `mean_d_own`, `mean_d_other`), test-set link |
| `results/cleanx/E4/E4_calibration.csv` | ECE/Brier by predicted and true tier |
| `results/cleanx/E4/E4_risk_control.csv` | realized coverage and selective risk per tier, with `errors`, `risk_ci_lo`, `risk_ci_hi` |
| `results/cleanx/E4/E4_min_certifiable_risk.csv`, `E4_openset.csv`, `E4_fallback.csv` | smallest certifiable risk per tier; abstention on unseen-EC (pure/mixed) and non-enzymes, AUROCs; EC3 fallback |
| `results/cleanx/E5/` | `E5_validation.csv` (gate), `E5_cards.csv`, `E5_pairs.csv` (**H5**), `E5_rank_agreement.csv` (τ_b, ρ), `E5_invariance.csv`, `E5_benchmark_card.csv` |
| `results/cleanx/E6/` | `E6_nstar_grid.csv`, `E6_tier_budget.csv` (**H6**: deficit and CI per tier; `years_at_deficit_ci_lo`), `E6_accrual.csv`, `E6_sgr_simulation.csv` |
| `results/cleanx/*_pilot/`, `*_smoke/`, `E5_local/`, `E5_card/` | pilot, smoke and local runs (never the paper's confirmatory numbers) |
| `results/cleanx/figures/*.png/.pdf/.csv` | figures and the numbers behind each, CIs included (`fig10_primary_endpoint` = E2 forest plot; `fig8`/`fig9` from E5/E6; `figS2`/`figS3` from the pilots); paper numbering in PUBLICATION_PLAN §8 |
| `data/cleanx/identity/blast/` | BLASTp hit tables and `manifest.json` |
| `data/cleanx/identity/<set>.csv`, `sp_temporal_test_vs_sp_temporal_cal.csv` | identity to training / TEST-vs-CAL identity (`tool` column says which aligner) |
| `data/cleanx/downloads/ec-bench/` | EC-Bench predictions at commit `5e98a2c` (MIT) with URL, date, SHA-256 |

## 7. Troubleshooting

- `An Application Control policy has blocked this file` → you ran the venv's python directly; use `.\cx.bat`.
- `'cx.bat' is not recognized` → you typed `cx.bat` without `.\`.
- `BLAST table mismatch` or a `FAIL` from `cx_11_blast.py` → your local training or set sequences differ from
  the ones searched; stop and send the output (do not edit the tables).
- `E5 STOPPED: reproduction gate failed` or a SHA-256 error → stop and send the output.
- `git pull` reports a conflict → do not resolve it by hand; send me the message.
- Never pull or edit files while a `.bat` is running: cmd.exe reads the batch file as it goes.
- `CLEANX_APP_DIR` is set only by the integration test; if it is set in your own shell, remove it.
- Keep `data/cleanx/`, `results/` and CLEAN's own files out of any **public** repository (CLEAN licence).
