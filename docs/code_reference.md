# Code reference

These module docstrings were relocated unchanged from commit `874abff997bc766e160f0b45869e21a4dbc3e187` during a
documentation-only cleanup. The original wording, dates and qualifications are preserved
below; this relocation is not a new registration or a revision of the recorded protocols.
Historical interpretations may differ from the final manuscript. Existing result digests
refer to the earlier code versions; no analyses were rerun for this cleanup.

<a id="blast"></a>

## `app/cleanx/blast.py`

```text
BLASTp homology baseline and sequence identity (NCBI BLAST+), the baseline used by CLEAN and by EC-Bench.

Search: blastp with default scoring (BLOSUM62, gap 11/1, composition-based statistics), E-value <= 10, up to
MAX_TARGET_SEQS subjects per query, against CLEAN's split100 training sequences. Per query, the hits are reduced
to the KEEP best by bit score plus the hit with the highest identity over the shorter sequence; that keeps both
quantities cx_04 needs exactly (top hit, max identity) while the stored tables stay small. Kept lines are written
exactly as BLAST printed them, gzipped without a timestamp (byte-identical on re-runs).

Top hit = highest bit score; ties broken by lower E-value, then subject accession (deterministic).
pid_short = identical residues / length of the shorter sequence, maximised over all HSPs (conservative against
short local hits; same definition as the Python fallback in cleanx/identity.py).
```

<a id="certify"></a>

## `app/cleanx/certify.py`

```text
Certification budget for selective risk control: how many calibration proteins does a group (e.g. a rarity
tier of the predicted EC) need before a rule with selective risk <= alpha can be certified?

Exact binomial facts used (one-sided Clopper-Pearson bound UCB, as in calibration.binom_ucb):
  duality:   UCB(k, n, delta) < alpha  <=>  P[X <= k] < delta  for X ~ Bin(n, alpha)
  hence      k_max(n) = min{k : P[X <= k] >= delta} - 1   is the largest certifiable number of errors
             (-1 = nothing certifiable at this n);
  power:     if the n accepted calibration proteins have true selective risk r, the rule is certified with
             probability power(n; r) = P[K <= k_max(n)], K ~ Bin(n, r);
  r = 0:     power = 1 as soon as 1 - delta**(1/n) < alpha, i.e. n > ln(delta) / ln(1 - alpha);
  r >= alpha: power -> 0 as n grows (not certifiable at any budget).
power(n; r) is not monotone in n (binomial saw-tooth), so n_first is the smallest n reaching the target and
n_stable the smallest n after which it never drops below the target again (up to n_max); budgets use n_stable.
SGR (Geifman & El-Yaniv 2017) tests ceil(log2 m) thresholds on m calibration points with delta/ceil(log2 m) each,
and the per-group (Mondrian) version splits delta over G groups, so a group certifies with per-test
delta' = delta / (G * ceil(log2 m)); group budgets are solved as a fixed point in m. The budget assumes the
accepted fraction (coverage c) is fixed in advance at its oracle threshold; cx_10 checks by simulation how
close SGR's own threshold search comes to it.
```

<a id="contrast"></a>

## `app/cleanx/contrast.py`

```text
Paired method contrasts within strata of proteins, and differences between strata (E7, E10).

delta(S) = metric(method A) - metric(method B) on the proteins of stratum S; both methods are scored on the same
proteins, so the bootstrap keeps them paired. The difference-in-differences DiD = delta(S1) - delta(S0) compares
two disjoint strata; its bootstrap resamples each stratum independently (a stratified bootstrap).

Identity standardisation (direct standardisation) re-weights the per-bin deltas of every stratum to one common bin
distribution w (the pooled distribution of the strata compared): delta_std(S) = sum_b w_b delta(S, b). Bins too
sparse to estimate are merged first by a fixed rule (merge_bins). Its bootstrap resamples every (stratum, bin) cell
independently, so no cell can be empty in a resample.
```

<a id="evalcard"></a>

## `app/cleanx/evalcard.py`

```text
EC evaluation card: every averaging scheme under both label-space conventions, the false positives each
one cannot see, and a numerical check of the invariance results below.

Notation. Proteins i = 1..n with true EC sets T_i and predicted sets P_i (EC strings truncated to the evaluated
depth). L_T = union of the T_i (the evaluated label set); L_U = union of the T_i and P_i. For a class c:
TP_c = #{i: c in T_i and P_i}, FP_c = #{i: c in P_i, not T_i}, FN_c = #{i: c in T_i, not P_i}, support
s_c = TP_c + FN_c, S = sum_c s_c. An out-of-label-set augmentation adds sets A_i with A_i disjoint from L_T.

Proposition 1 (invariance). Under any out-of-label-set augmentation, (TP_c, FP_c, FN_c) is unchanged for every
c in L_T. Hence every score that depends only on these counts for c in L_T is unchanged, in particular
  (i)  every averaging (micro, macro, weighted, samples) computed after projecting predictions onto L_T
       (the protocol of CLEAN's get_eval_metrics, which fits its MultiLabelBinarizer on the test labels), and
  (ii) the support-weighted and the true-label-macro average computed on L_U (e.g. EC-Bench's scoring),
       because every class outside L_T has s_c = 0 and therefore weight s_c / S = 0, and is excluded from the
       true-label macro mean.
Proof: A_i is disjoint from L_T, so for c in L_T membership of c in P_i is unchanged; for (ii) the classes that
change all have zero support.

Proposition 2 (sensitivity). If some A_i is non-empty and sum_c TP_c > 0, micro precision on L_U strictly
decreases; samples-averaged precision strictly decreases if some protein with a correct prediction receives an
added label; the union-label macro F1 cannot increase, and strictly decreases when it was positive and an added
class was not yet in L_U (the class enters the mean with F1 = 0; classes outside L_T already had F1 = 0).

Lemma 3. On any label space containing L_T, support-weighted recall equals micro recall
(sum_c (s_c/S)(TP_c/s_c) = sum_c TP_c / S). The whole difference between the reported and the micro scores is
therefore on the precision side.

Proposition 4 (over-prediction). Adding a correct label c to P_i raises TP_c and lowers FN_c, so it weakly
raises P_c, R_c and F1_c; adding a wrong label inside L_T weakly lowers F1_c; adding a wrong label outside L_T
changes nothing (Prop. 1). Top-k expansion therefore raises the support-weighted F1 whenever its gains in
recall outweigh the costs of the wrong labels that happen to fall inside L_T, and the smaller L_T is relative to
the predictor's vocabulary, the larger the share of wrong labels that are free.
```

<a id="freeze"></a>

## `app/cleanx/freeze.py`

```text
Pre-registration guard: a confirmatory run must use analysis code that is committed, unchanged.

code_state(paths) asks git whether any of the given files differs from HEAD or is untracked (a file outside the
repository never counts as committed). Scripts call require_frozen() before a confirmatory run: it stops with a
message unless the code is committed, and --allow-dirty (development and tests only) turns the stop into a loud
label that is written into every output of that run. digests(paths) gives the SHA-256 of every analysis file; the
scripts write them into their outputs so that a run can be matched to a commit afterwards, with or without git.
```

<a id="identity"></a>

## `app/cleanx/identity.py`

```text
Maximum sequence identity of query proteins to the CLEAN training set, and a homology-transfer
baseline (EC of the best local-alignment hit).

pid_short = identical aligned residues / length of the shorter sequence (conservative against short
local hits); pid_aln = identical / alignment columns (MMseqs2 'fident').
```

<a id="metrics"></a>

## `app/cleanx/metrics.py`

```text
Multi-label EC metrics with an explicit label space.

'weighted' reproduces the number CLEAN reports (sklearn support-weighted average over classes).
Classes that are predicted but have no true instance get weight 0 there, so false positives on EC
numbers absent from the evaluated set's labels never reduce it; 'macro_true' (mean over classes
with a true instance) shares that blind spot. Both can be raised by over-predicting. 'micro',
'macro_union' (mean over true + predicted classes) and 'samples' count every prediction.
```

<a id="provenance"></a>

## `app/cleanx/provenance.py`

```text
Label provenance (E7): how was each benchmark EC label on a UniProtKB entry asserted?

UniProtKB attaches Evidence and Conclusion Ontology (ECO) codes to every annotation. For an EC number the evidence
is the union of the codes on that EC in the protein description (recommended, alternative and submission names,
including 'includes'/'contains' domains) and on CATALYTIC ACTIVITY reactions that carry that EC:
  ECO:0000269  experimental evidence used in manual assertion
  ECO:0000250  sequence similarity evidence used in manual assertion ('by similarity'); its source names the
               protein the annotation was propagated from (source 'UniProtKB', id = accession)
  ECO:0000255 / ECO:0000256 / ECO:0000259  match to a sequence model (e.g. a HAMAP rule) used in manual /
               automatic assertion
Everything else (e.g. ECO:0000305 curator inference, ECO:0000303 author statement, imported annotations, no code)
is OTHER.

Classes (fixed in the E7 pre-registration, cx_13_provenance.py):
  EC level      EXP if any ECO:0000269; else SIM if any ECO:0000250; else RULE if any sequence-model code;
                else OTHER
  protein level LABEL_ABSENT if the entry was not found or any benchmark label is not an EC on the current entry;
                else SIM if any label is SIM; else EXP if every label is EXP; else RULE if any label is RULE;
                else OTHER. SIM and RULE together are the labels propagated by homology (PROPAGATED).

Accession resolution: an active entry whose primary accession is the query; else the single active entry that
lists the query as a secondary accession (more than one such entry: 'ambiguous_secondary', not found); else an
inactive entry for the query that was merged into exactly one active entry present in the responses ('merged');
otherwise not found ('inactive:<reason>' or 'absent').

This module parses the UniProt REST JSON (uniprotkb entries) into flat tables; it never downloads anything.
```

<a id="cx-01-audit"></a>

## `app/cx_01_audit.py`

```text
E1: evaluation validity of the CLEAN baseline.

- every averaging scheme at EC depth 1-4, with paired bootstrap CIs
- the support-weighted 'blind spot' (false positives on ECs outside the evaluated label set)
- tier x depth breakdown (reproduces app/audit_out/maxsep_tier_x_depth.csv for New-392)
- exact-sequence leakage between test sets and split100
- EC label drift of the legacy benchmarks against current Swiss-Prot (needs cx_02 --with-pfam table)

```

<a id="cx-04-identity"></a>

## `app/cx_04_identity.py`

```text
Max sequence identity of each evaluation protein to CLEAN's split100 training set, plus the
homology-transfer baseline (EC of the best alignment hit).

--tool python : 3-mer prefilter + Smith-Waterman (Biopython) on the 30 best candidates, runs anywhere
--tool blast  : NCBI BLASTp tables from cx_11_blast.py (the paper's baseline, as in CLEAN and EC-Bench); refuses to
                run unless the local sequences are exactly the ones searched (SHA-256 digests in the manifest)
--tool m8     : parse MMseqs2 results
--write-fasta : write split100 + query FASTAs for running MMseqs2 elsewhere
--pair Q:T    : max identity of every protein of set Q to set T (e.g. sp_temporal_test:sp_temporal_cal), written
                to data/cleanx/identity/Q_vs_T.csv; used by cx_05 for the CAL-disjoint TEST sensitivity analysis
                (sequences only: no prediction or label of either set is used)
```

<a id="cx-08-figures"></a>

## `app/cx_08_figures.py`

```text
Paper figures from the E1-E6 outputs. Each figure is written as PNG (300 dpi) + PDF, with a CSV
of exactly the plotted numbers next to it. Figures whose inputs do not exist yet are skipped.

Every interval drawn is a 95% interval: the paired percentile bootstrap over proteins of cx_05 (B = 2000, seed
config.SEED) for F1, exact Clopper-Pearson for accuracies and error rates, and the bootstrap interval cx_10 reports
for the calibration deficit. F1 point estimates are recomputed from E2_per_query.csv and must equal
E2_metrics_long.csv; a figure is skipped, with the reason printed, rather than drawn from numbers that disagree.
```

<a id="cx-11-blast"></a>

## `app/cx_11_blast.py`

```text
BLASTp (NCBI BLAST+) for the identity strata and the homology-transfer baseline.

Why. The Python fallback in cleanx/identity.py aligns only the 30 best 3-mer candidates per query and has no
E-value. CLEAN's paper and EC-Bench use BLASTp as the homology baseline, and RUNBOOK section 4 (committed with the
pre-registration, c8e42f7) planned to replace the fallback with a standard aligner for the paper's numbers. This
changes a comparator and the identity strata; no pre-registered endpoint uses either.

Search (cleanx/blast.py): blastp, default scoring (BLOSUM62, gap 11/1, composition-based statistics), E-value <= 10,
up to 500 subjects per query; database = CLEAN's split100 training sequences; queries = New-392, Price-149,
SP-temporal CAL and TEST; plus TEST against a CAL database for the CAL-disjoint sensitivity subset. Per query the
10 best hits by bit score and the maximum-identity hit are kept (data/cleanx/identity/blast/<set>.tsv.gz).
manifest.json records the BLAST version, the exact commands, the parameters and SHA-256 digests of every sequence
set and output file. cx_04_identity.py --tool blast refuses to use a table whose sequences differ from local data.

BLAST+ may be blocked by Application Control on the analysis PC, so the search runs wherever BLAST+ is installed
(done in a Linux copy of the environment) and the tables are committed; the PC only verifies and parses them.

Modes: default (search everything; needs blastp), --verify (digests of local data and of the stored tables against
the manifest; no BLAST needed; PASS/FAIL), --smoke (verify, then, if blastp is available, re-search 5 New-392
queries and require the same top hit, bit score and maximum identity as the stored table; otherwise SKIP).
```

<a id="cx-12-uniprot-evidence"></a>

## `app/cx_12_uniprot_evidence.py`

```text
E7 data: UniProtKB evidence behind the EC labels of the benchmark proteins (label provenance).

For every protein of New-392 (development set), SP-temporal CAL and TEST, and the EC-Bench test set, download the
current UniProtKB entry (REST API, JSON) and record, for each EC number on the entry, the ECO evidence codes and
the accessions cited by 'by similarity' (ECO:0000250) evidence (cleanx/provenance.py). No model output is read:
this script only describes the labels. The analysis is cx_13_provenance.py.

Download: accessions in batches of 100 through the UniProtKB stream endpoint (query 'accession:A OR accession:B
...', format=json, gzip); accessions a batch does not resolve are fetched one by one from the single-entry endpoint;
an inactive entry that was merged into one other entry is followed to that entry. Accessions are resolved by the
rules in cleanx/provenance.py (primary, unique secondary, merged); the rest are recorded as not found. Responses
are cached under data/cleanx/downloads/uniprot_json/<set>/ (batch files are named by the SHA-256 of their
accession list, so a changed set never reuses a stale response) and every cached file, the UniProt release it came
from and the counts are recorded with SHA-256 digests in data/cleanx/provenance/manifest.json. Runs on a machine
with internet access (the analysis PC); --offline re-parses the cache without any network.

Outputs: data/cleanx/provenance/<set>_entries.csv (one row per queried accession) and <set>_ec_evidence.csv (one
row per EC on the entry: codes, similarity sources, where the EC was found).

Modes: default (download missing responses, then parse), --offline (parse the cache only), --smoke (parser
known-answer checks on synthetic entries; no network), --limit N (first N accessions per set, for a live check;
writes <set>_limit<N> tables that the analysis never reads).
```

<a id="score-audit"></a>

## `app/score_audit.py`

```text

Stage C scorer for the CLEAN audit.

Reads:
  - ground truth        (data/new.csv        : tab-separated Entry / EC number / Sequence)
  - CLEAN predictions   (results/new_maxsep.csv)
  - EC frequency index  (data/ec_train_frequency.csv)

Reports multi-label precision / recall / F1 broken down by:
  (a) EC digit depth   1 -> 4   (does it get the family but miss the reaction?)
  (b) rarity tier      common -> singleton -> unseen
  (c) both together    (the core audit table)

Usage:
    python score_audit.py
    python score_audit.py --preds results/new_pvalue.csv --tag pvalue

```

<a id="run-tests"></a>

## `app/tests/run_tests.py`

```text
Unit and integration tests for cleanx (plain Python; pytest is not needed).

Run from app\:   .\cx.bat tests\run_tests.py            (unit + integration, ~2-4 min)
                 .\cx.bat tests\run_tests.py --unit-only
Prints [PASS]/[FAIL] per test and ALL <n> TESTS PASSED at the end; exits 1 on any failure. The integration test
builds a synthetic app in a temporary folder and never writes into data/ or results/ of this repository.
```

<a id="synth-app"></a>

## `app/tests/synth_app.py`

```text
Builds a small synthetic app/ folder (data + weights) in a temporary directory so the real cx_* scripts can run
end to end without CLEAN's weights, ESM-1b or any download. Nothing here reads or writes the real repository data.

The world is coherent: training ESM vectors are clustered by EC, 100.pt is those vectors passed through a
(randomly initialised, seeded) CLEAN LayerNormNet saved as split100.pth, and query vectors lie near their EC's
centre, so max-separation, calibration and risk control have something real to find.
```

<a id="test-integration"></a>

## `app/tests/test_integration.py`

```text
End-to-end integration test on synthetic data: builds a synthetic app in a temporary folder (tests/synth_app.py)
and runs the real scripts against it through CLEANX_APP_DIR. Covers every new or modified path: cx_04 --pair,
cx_05 (observed deltas, Holm, sensitivity subsets), cx_06 (d_own / d_other), cx_07 (pure / mixed open set, EC3
fallback), cx_09 --skip-ecbench and --card, cx_10 (certification budget) and cx_08 (figures, including the
fallback to pilot results when the confirmatory folder is incomplete). Never touches the real data/
or results/ folders (checked).
```
