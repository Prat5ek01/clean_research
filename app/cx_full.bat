@echo off
setlocal
cd /d "%~dp0"
rem Confirmatory run: commit the pre-registered code first (PUBLICATION_PLAN section 7) and cite that hash.
echo ==================== 1/14  ESM-1b equivalence test + embeddings (resumable; ~6 h the first time) ====================
call "%~dp0cx.bat" cx_03_embed.py || goto fail
echo ==================== 2/14  BLASTp tables match the local sequences (no BLAST needed) ====================
call "%~dp0cx.bat" cx_11_blast.py --verify || goto fail
echo ==================== 3/14  identity strata + BLASTp homology baseline ====================
call "%~dp0cx.bat" cx_04_identity.py --tool blast --sets new392 price149 sp_temporal_cal sp_temporal_test || goto fail
echo ==================== 4/14  TEST vs CAL identity, BLASTp (sensitivity subset) ====================
call "%~dp0cx.bat" cx_04_identity.py --tool blast --pair sp_temporal_test:sp_temporal_cal || goto fail
echo ==================== 5/14  E2 inference geometry, split100 (primary endpoint) ====================
call "%~dp0cx.bat" cx_05_geometry.py || goto fail
echo ==================== 6/14  E2 robustness: split70 backbone ====================
call "%~dp0cx.bat" cx_05_geometry.py --backbone split70 || goto fail
echo ==================== 7/14  E3 NISE biology ====================
call "%~dp0cx.bat" cx_06_ec_biology.py || goto fail
echo ==================== 8/14  E4 reliability ====================
call "%~dp0cx.bat" cx_07_reliability.py || goto fail
echo ==================== 9/14  E1 audit incl. temporal test ====================
call "%~dp0cx.bat" cx_01_audit.py --datasets new392 price149 sp_temporal_test || goto fail
echo ==================== 10/14  E5 EC-Bench re-scoring + temporal evaluation cards ====================
call "%~dp0cx.bat" cx_09_evalcard.py || goto fail
echo ==================== 11/14  E6 certification budget (SP-temporal CAL) ====================
call "%~dp0cx.bat" cx_10_certbudget.py || goto fail
echo ==================== 12/14  E7 data: UniProt evidence of every benchmark label (internet; cached) ====================
call "%~dp0cx.bat" cx_12_uniprot_evidence.py || goto fail
echo ==================== 13/14  E7 label provenance (confirmatory; refuses uncommitted analysis code) ====================
call "%~dp0cx.bat" cx_13_provenance.py || goto fail
echo ==================== 14/14  figures ====================
call "%~dp0cx.bat" cx_08_figures.py || goto fail
echo.
echo FULL PIPELINE FINISHED - results in results\cleanx\
exit /b 0
:fail
echo.
echo FAILED in the step printed above. Stop here and send the output.
exit /b 1
