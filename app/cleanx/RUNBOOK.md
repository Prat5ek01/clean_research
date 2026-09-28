# cleanx RUNBOOK — exact commands, smoke tests, pass criteria

Companion to [`../../PUBLICATION_PLAN.md`](../../PUBLICATION_PLAN.md) (the *why*). This file is the *how*.
Last updated 2026-09-28 (E5 evaluation card / EC-Bench re-scoring, E6 certification budget, sensitivity
analyses, unit and integration tests).

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
biopython 1.83 (added for the identity step with `--no-deps`; nothing else changed). The 2026-09-28
additions need no new packages.

## 1. The commands you need, in order

### Command 1 — smoke test of everything (≈4–6 min)

```bat
.\cx_smoke.bat
```

**PASS** = the last line is `ALL 10 SMOKE TESTS PASSED`. It stops at the first failure and prints
`SMOKE TEST FAILED`. Expected lines along the way:

| Step | What | Expected line(s) | Time |
|---|---|---|---|
| 1 | environment + exact CLEAN reproduction | `CX00 SUMMARY: 9/9 checks passed` (P 0.5955, R 0.4791, F1 0.4967; 0/392 prediction mismatches) | steps 1–7 together ≈2 min on your PC (2026-09-27) |
| 2 | E1 audit (smoke) | `E1 DONE` | |
| 3 | temporal benchmark builder (cached downloads) | `CX02 DONE` | 4 s * |
| 4 | identity search, 5 proteins | `homology top-hit exact-EC agreement: 5/5`, `CX04 DONE` | 19 s * |
| 5 | E2 geometry (smoke) | `exemplar - centroid micro_F1: +0.034 [+0.003,+0.070]` ¹, `E2 DONE` | |
| 6 | E3 NISE biology (smoke) | `E3 DONE` | |
| 7 | E4 reliability (smoke) | `E4 DONE` | |
| 8 | unit + integration tests (synthetic data, temporary folders only) | `ALL 22 TESTS PASSED` | 46 s * |
| 9 | E5 evaluation card (smoke) | `[PASS] EC-Bench reported weighted F1 reproduced (32 values, test set, 100% regime) -- max \|diff\| 0.0e+00` ², `[PASS] Propositions 1-2 and Lemma 3 hold numerically on every prediction table -- 29/29 tables`, `CX09 SMOKE DONE` | 2 s * |
| 10 | E6 certification budget (smoke) | three `[PASS]` lines (`n*(r=0, alpha=0.2, delta=0.025) = 17`; no budget at or above alpha; n* non-decreasing), `CX10 SMOKE DONE` | 75 s * |

\* measured in the Linux copy on 2026-09-28. Every script prints its own runtime in its `DONE in Ns` line.

¹ Since 2026-09-28 the summary prints the observed full-sample difference (+0.034) instead of the bootstrap
mean (+0.035, as before); the confidence interval is unchanged. The line
`WARNING: new392 is also evaluated; its balanced* results are optimistic` is expected in the smoke run
(λ is tuned on New-392 there only to exercise the code).
² Only when the EC-Bench files are present in `data\cleanx\downloads\ec-bench\` (they are in the
repository once the 2026-09-28 changes are committed); otherwise step 9 prints `reproduction gate skipped in smoke mode` and still passes. The table
count after `N/N` can differ if you have more prediction files; it must always read `N/N`.

Steps 3, 4, 8, 9 and 10 were run on 2026-09-28 in a Linux copy of this environment (same package versions, CUDA build of torch 1.11 run on CPU);
steps 1, 2, 5, 6 and 7 need CLEAN's weights and embeddings and have only been verified on your PC with the
2026-09-27 code. If any step-5 number other than the point estimate differs, stop and send me the output.

### Command 2 — ESM-1b equivalence test (≈5 min; needs ~8–10 GB free RAM)

Close browsers/IDEs first (Task Manager → Memory). Then:

```bat
.\cx.bat cx_03_embed.py --verify-only
```

**PASS** = a line `[PASS] recomputed ESM-1b embeddings match CLEAN's stored ones for 3 New-392
proteins` with max |diff| < 1e-3 and min cosine > 0.99999. The next line prints the CPU speed in
**s/protein**; the full embedding run takes about 4,099 × that. **If it says FAIL, stop and send me
the output**: the rest of the pipeline assumes identical embeddings.

### Command 3 — record the pre-registration commit (before Command 4)

The confirmatory run must use code that is already committed, so the Methods can cite the hash. The
pre-registrations are in the docstrings of `cx_05_geometry.py`, `cx_09_evalcard.py`, `cx_10_certbudget.py`
and in PUBLICATION_PLAN §7. Once the 2026-09-28 changes are committed (and pulled to this PC, if they were
committed elsewhere):

```bat
git status
git log -1 --format="%H %ad"
```

**PASS** = `git status` shows no modified `.py`/`.bat` files under `app\`, and you have written down the hash.
Do not edit any script between this point and the end of Command 4.

### Command 4 — full pipeline (≈2.5–4.5 h, almost all of it the embedding step; resumable)

```bat
.\cx_full.bat
```

**PASS** = the last line is `FULL PIPELINE FINISHED`. If the PC sleeps or the run stops, run `.\cx_full.bat`
again: finished embeddings are skipped, and every analysis step is fast. Steps and what each must print:

| Step | Command inside `cx_full.bat` | PASS line | Time |
|---|---|---|---|
| 1 | `cx_03_embed.py` (repeats the equivalence test first) | `[PASS] recomputed ESM-1b embeddings match ...`, `CX03 DONE: 0 embeddings missing` | 2–4 h (≈ 4,099 × the s/protein of Command 2) |
| 2 | `cx_04_identity.py --pair sp_temporal_test:sp_temporal_cal` | `CX04 PAIR DONE`; writes `data\cleanx\identity\sp_temporal_test_vs_sp_temporal_cal.csv` | ≈5 min * |
| 3 | `cx_05_geometry.py` (split100: **E2 primary endpoint**) | `E2 DONE` | not timed on real data |
| 4 | `cx_05_geometry.py --backbone split70` | `E2 DONE` | not timed on real data |
| 5 | `cx_06_ec_biology.py` | `E3 DONE` | not timed on real data |
| 6 | `cx_07_reliability.py` | `E4 DONE` | not timed on real data |
| 7 | `cx_01_audit.py --datasets new392 price149 sp_temporal_test` | `E1 DONE` | not timed on real data |
| 8 | `cx_09_evalcard.py` (**E5 confirmatory**) | `[PASS] EC-Bench reported weighted F1 reproduced (32 values ...)`, `[PASS] Propositions 1-2 and Lemma 3 hold ... N/N tables`, `E5 DONE` | ≈1–2 min (pilot on 2 of the 4 files: 26 s *) |
| 9 | `cx_10_certbudget.py` (**E6**, calibration = SP-temporal CAL) | three `[PASS]` known-answer lines, `E6 DONE` | ≈2–10 min (pilot: 116 s *) |
| 10 | `cx_08_figures.py` | `wrote fig8_ecbench_rescoring.png/.pdf/.csv`, `wrote fig9_certification_budget.png/.pdf/.csv`, `CX08 DONE` | < 1 min |

\* measured in the Linux copy on 2026-09-28; the other analysis steps were fast on 2026-09-27 and print their
runtime.

If step 8 prints `E5 STOPPED: reproduction gate failed`, or any script raises a SHA-256 error, the run stops
on purpose: send me the output and do not re-run with changed settings.

Then send me these files (or paste their contents) and I will fill the numbers into the paper:

- E2: `results\cleanx\E2_split100\E2_summary.md`, `E2_deltas.csv`; `results\cleanx\E2_split70\E2_summary.md`
- E3: `results\cleanx\E3_split100\E3_summary.md`
- E4: `results\cleanx\E4\E4_summary.md`, `E4_min_certifiable_risk.csv`, `E4_openset.csv`
- E5: `results\cleanx\E5\E5_summary.md`, `E5_validation.csv`, `E5_rank_agreement.csv`
- E6: `results\cleanx\E6\E6_summary.md`, `E6_tier_budget.csv`, `E6_sgr_simulation.csv`
- figures: `results\cleanx\figures\fig8_ecbench_rescoring.png`, `fig9_certification_budget.png`

### Optional before Command 4 — refresh the legacy E2/E3 outputs with the new columns (minutes)

The committed `E2_split100` and `E3_split100` results come from the 2026-09-27 code. To add the new columns
(observed deltas, Holm-adjusted tier p-values, sensitivity subsets, mean distance to own / other EC) without
waiting for the embeddings:

```bat
.\cx.bat cx_05_geometry.py --datasets new392 price149
.\cx.bat cx_06_ec_biology.py
```

**PASS** = the two runs end with `E2 DONE` and `E3 DONE`. Per-method scores and all confidence intervals must
not change (only columns and rows are added); the paired-delta point estimates in `E2_summary.md` now show the
observed difference instead of the bootstrap mean, so they can move in the third decimal. Anything else that
changes: send me both summaries. Command 4 re-runs both anyway.

## 2. What has been run and verified

| Step | Script | Status | Output |
|---|---|---|---|
| E0 reproduction checks | `cx_00_check.py` | ✅ 9/9 PASS (your PC, 2026-09-27) | console |
| E1 evaluation audit (New-392, Price-149) | `cx_01_audit.py` | ✅ (your PC) | `results/cleanx/E1/` |
| Temporal benchmark build | `cx_02_build_temporal.py --with-pfam` | ✅ UniProt 2026_03, ENZYME 02-Sep-2026, Pfam 38.2 | `data/cleanx/sets/`, `downloads/` |
| Identity + homology baseline | `cx_04_identity.py` | ✅ (Python fallback) | `data/cleanx/identity/` |
| E2 geometry on New-392, Price-149 | `cx_05_geometry.py --datasets new392 price149` | ✅ 2026-09-27 code; ⏳ refresh optional (above) | `results/cleanx/E2_split100/` |
| E3 NISE biology (legacy link) | `cx_06_ec_biology.py` | ✅ 2026-09-27 code; ⏳ refresh optional (above) | `results/cleanx/E3_split100/` |
| Figures from legacy results | `cx_08_figures.py` | ✅ checked visually | `results/cleanx/figures/` |
| Smoke suite, 7 original steps | `cx_smoke.bat` | ✅ ALL 7 PASSED (your PC, 2026-09-27) | — |
| Unit + integration tests (22) | `tests\run_tests.py` | ✅ ALL 22 PASSED (Linux copy, 2026-09-28) | temporary folders only |
| Smoke steps 3, 4, 8, 9, 10 | `cx_smoke.bat` | ✅ (Linux copy, 2026-09-28); steps 1, 2, 5, 6, 7 ⏳ your PC | — |
| EC-Bench reproduction gate | `cx_09_evalcard.py` | ✅ 32/32 printed values reproduced, max \|diff\| 0.0 | `E5_pilot/E5_validation.csv` |
| E5 **pilot** (EC-Bench Price-149, disclosed as pilot) | `cx_09_evalcard.py --pilot` | ✅ deterministic (identical under 3 hash seeds) | `results/cleanx/E5_pilot/` |
| E6 **pilot** (E4 smoke features as calibration data) | `cx_10_certbudget.py --pilot` | ✅ deterministic (identical on re-run) | `results/cleanx/E6_pilot/` |
| TEST vs CAL identity | `cx_04_identity.py --pair sp_temporal_test:sp_temporal_cal` | ✅ (Linux copy); Command 4 recomputes it | `data/cleanx/identity/` |
| **ESM-1b load itself** | `cx_03_embed.py --verify-only` | ⏳ **your Command 2** | — |
| **E5 confirmatory, E6 on CAL** | `cx_09_evalcard.py`, `cx_10_certbudget.py` | ⏳ **Command 4, after Command 3** | `results/cleanx/E5/`, `E6/` |

## 3. Individual steps (only if you want to run pieces by hand)

```bat
.\cx.bat cx_00_check.py
.\cx.bat cx_03_embed.py --limit 5
.\cx.bat cx_03_embed.py
.\cx.bat cx_04_identity.py --pair sp_temporal_test:sp_temporal_cal
.\cx.bat cx_05_geometry.py
.\cx.bat cx_05_geometry.py --backbone split70
.\cx.bat cx_05_geometry.py --cap-per-ec 3
.\cx.bat cx_06_ec_biology.py
.\cx.bat cx_07_reliability.py
.\cx.bat cx_01_audit.py --datasets new392 price149 sp_temporal_test
.\cx.bat cx_09_evalcard.py
.\cx.bat cx_10_certbudget.py
.\cx.bat cx_08_figures.py
```

`--cap-per-ec` changes only the capped sensitivity subset (default 5, the pre-registered value; other values
are exploratory). New modes added 2026-09-28:

| Command | What it does | ≈ time |
|---|---|---|
| `.\cx.bat tests\run_tests.py` | all 22 tests (unit + integration on a synthetic app in a temporary folder) | 1 min |
| `.\cx.bat tests\run_tests.py --unit-only` | the 21 unit tests only | 2 s |
| `.\cx.bat cx_09_evalcard.py --smoke` | E5 on synthetic EC-Bench-format data (+ the reproduction gate if the files are present) → `E5_smoke` | 5 s |
| `.\cx.bat cx_09_evalcard.py --pilot` | EC-Bench Price-149 only (disclosed pilot) → `E5_pilot` | 30 s |
| `.\cx.bat cx_09_evalcard.py --skip-ecbench` | evaluation cards for our own predictions only (legacy + temporal), no download → `E5_local` | 5 s |
| `.\cx.bat cx_09_evalcard.py --card data\new.csv results\cleanx\preds\new392_split100_centroid_maxsep.csv` | evaluation card for any prediction file (CLEAN format `id,EC:x/score,...` or `id,EC;EC`) → `E5_card\` | 1 s |
| `.\cx.bat cx_10_certbudget.py --smoke` | E6 on synthetic features + known-answer checks → `E6_smoke` | 1–2 min |
| `.\cx.bat cx_10_certbudget.py --pilot` | E6 with the E4 smoke features (New-392 half) as calibration data → `E6_pilot` | 2 min |

The `--card` example above must print `reported (support-weighted) P/R/F1 0.596/0.479/0.497 | micro
0.518/0.479/0.498 | union-macro F1 0.371 | per-protein F1 0.535 | 150/224 FPs invisible` (checked 2026-09-28).

Optional E3b (mechanism: are non-homologous subgroups separate clusters in CLEAN space? ≈1.5–2 h):

```bat
.\cx.bat cx_03_embed.py --sets nise_members
.\cx.bat cx_06_ec_biology.py --subgroups
```

Do **not** run `cx_02_build_temporal.py --force`: it re-downloads UniProt, and a newer release
changes the benchmark. Without `--force`, cached downloads are reused. Likewise, `cx_09_evalcard.py
--force-download` is only for a corrupted EC-Bench file: the files are pinned to EC-Bench commit `5e98a2c` and
checked against fixed SHA-256 values, so a new download must be byte-identical.

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
identity-stratified results use them. (The TEST-vs-CAL table from `--pair` always uses the built-in
fallback; it only defines the `closed_cal_disjoint` sensitivity subset.)

## 5. Output files

| File | Contents |
|---|---|
| `results/cleanx/E1/E1_blind_spot.csv` | predictions, TP, FP, FPs invisible to the reported precision |
| `results/cleanx/E1/E1_overprediction.csv` | reported vs micro scores when 1–10 ECs are predicted per protein |
| `results/cleanx/E1/E1_exact_duplicates.csv`, `E1_label_drift.csv` | leakage and label drift |
| `results/cleanx/E2_*/E2_deltas.csv` | paired-bootstrap differences vs CLEAN centroid (**primary endpoint**: subset `closed`, `balanced3 - centroid`, `macro_union_F1`, `sp_temporal_test`). Report `delta_obs` (observed difference) with `ci_lo`/`ci_hi`; tier rows carry `p_holm`; sensitivity rows: `closed|cluster_bootstrap_EC`, `closed_capped5`, `closed_cal_disjoint`, `closed_PE1` |
| `results/cleanx/E2_*/E2_bootstrap_ci.csv` | per-method CIs; `observed` = full-sample value |
| `results/cleanx/E2_*/E2_tuning.csv/.json` | λ chosen on SP-temporal CAL |
| `results/cleanx/E2_*/E2_per_query.csv` | per-protein predictions and F1 for every method |
| `results/cleanx/E3_*/E3_stats.json` | NISE statistics, size-stratified permutation tests (now with `mean_d_own`, `mean_d_other`), test-set link |
| `results/cleanx/E4/E4_calibration.csv` | ECE/Brier by predicted and true tier |
| `results/cleanx/E4/E4_risk_control.csv`, `E4_min_certifiable_risk.csv` | certified thresholds; smallest certifiable risk per tier |
| `results/cleanx/E4/E4_openset.csv` | abstention on unseen-EC enzymes (all, `_pure`, `_mixed`) and non-enzymes; AUROCs |
| `results/cleanx/E4/E4_fallback.csv` | EC level-3 fallback for proteins the EC-4 rule rejects (threshold calibrated on the CAL proteins it rejects) |
| `results/cleanx/E5/E5_validation.csv` | reproduction gate: our scorer vs EC-Bench's 32 printed values |
| `results/cleanx/E5/E5_cards.csv` | evaluation card per method × dataset × regime × EC level: every metric family, reported (CLEAN and EC-Bench protocols) vs micro, invisible FP share |
| `results/cleanx/E5/E5_pairs.csv` | all 120 method pairs: reported vs micro differences, bootstrap CIs, BH q-values, reversals (**H5 primary**: test / 100% / EC4) |
| `results/cleanx/E5/E5_rank_agreement.csv` | Kendall τ_b (bootstrap CI) between reported and each honest metric; Spearman ρ(inflation, ECs per protein) with one-sided p (key secondary S4) |
| `results/cleanx/E5/E5_invariance.csv` | Props. 1–2 / Lemma 3 checked on every prediction table |
| `results/cleanx/E5/E5_benchmark_card.csv` | benchmark card for SP-temporal CAL/TEST (size, unseen share pure/mixed, largest EC) |
| `results/cleanx/E6/E6_nstar_grid.csv` | calibration proteins needed to certify selective risk α with power 0.8, by true risk |
| `results/cleanx/E6/E6_tier_budget.csv` | per predicted tier: available vs needed, deficit with bootstrap CI, `status`, years of Swiss-Prot accrual |
| `results/cleanx/E6/E6_accrual.csv`, `E6_sgr_simulation.csv` | Swiss-Prot enzymes per year passing the filters; SGR certification rate at 1× and 2× the budget |
| `results/cleanx/E5_card/E5_card_<name>.csv` | output of `--card` (kept out of the confirmatory `E5` folder) |
| `results/cleanx/*_pilot/`, `*_smoke/`, `E5_local/` | pilot, smoke and local-only runs (never the paper's confirmatory numbers) |
| `results/cleanx/figures/*.png/.pdf/.csv` | paper figures and the numbers behind each; `fig8`/`fig9` from E5/E6, `figS2`/`figS3` from the pilots |
| `data/cleanx/downloads/ec-bench/*.csv(.json)` | EC-Bench prediction files at commit `5e98a2c` (MIT licence) with URL, date and SHA-256 |
| `data/cleanx/identity/sp_temporal_test_vs_sp_temporal_cal.csv` | best identity of each TEST protein to any CAL protein |
| `data/cleanx/sets/sp_temporal_manifest.json` | how the benchmark was built (filters, counts, releases) |

## 6. Troubleshooting

- `An Application Control policy has blocked this file` → you ran the venv's python directly; use
  `.\cx.bat`.
- `'cx.bat' is not recognized` → you typed `cx.bat` without `.\`.
- `FileNotFoundError: … ESM-1b embeddings missing` → the embedding step is not finished for that set.
- ESM load very slow / PC freezes → not enough free RAM; close apps and re-run (resumable).
- `E5 STOPPED: reproduction gate failed` or `SHA-256 ... is not the pinned ...` → stop and send the output.
  The first confirmatory run needs internet access to `raw.githubusercontent.com` only if the EC-Bench files
  are missing from `data\cleanx\downloads\ec-bench\`.
- `skip fig8` / `skip fig9` from `cx_08_figures.py` → E5 / E6 have not finished yet; the pilot versions
  (`figS2`, `figS3`) are drawn instead when pilot results exist.
- A test fails in step 8 → send the output. The tests write only to a temporary folder and check that
  `results\cleanx` is unchanged; they never touch `data\esm_data`.
- `CLEANX_APP_DIR` is set only by the integration test to point the scripts at a synthetic app. If it is set
  in your own shell, remove it: it redirects every path.
- Git: keep `data/cleanx/`, `results/` and CLEAN's own files out of any **public** repository
  (PUBLICATION_PLAN §13, CLEAN licence).
