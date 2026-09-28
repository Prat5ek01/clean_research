@echo off
setlocal
cd /d "%~dp0"
rem Confirmatory run: commit the pre-registered code first (PUBLICATION_PLAN section 7) and cite that hash.
echo ==================== 1/10  ESM-1b equivalence test + embeddings (resumable, 2-4 h) ====================
call "%~dp0cx.bat" cx_03_embed.py || goto fail
echo ==================== 2/10  TEST vs CAL sequence identity (sensitivity subset, ~5 min) ====================
call "%~dp0cx.bat" cx_04_identity.py --pair sp_temporal_test:sp_temporal_cal || goto fail
echo ==================== 3/10  E2 inference geometry, split100 (primary endpoint) ====================
call "%~dp0cx.bat" cx_05_geometry.py || goto fail
echo ==================== 4/10  E2 robustness: split70 backbone ====================
call "%~dp0cx.bat" cx_05_geometry.py --backbone split70 || goto fail
echo ==================== 5/10  E3 NISE biology ====================
call "%~dp0cx.bat" cx_06_ec_biology.py || goto fail
echo ==================== 6/10  E4 reliability ====================
call "%~dp0cx.bat" cx_07_reliability.py || goto fail
echo ==================== 7/10  E1 audit incl. temporal test ====================
call "%~dp0cx.bat" cx_01_audit.py --datasets new392 price149 sp_temporal_test || goto fail
echo ==================== 8/10  E5 EC-Bench re-scoring + temporal evaluation cards ====================
call "%~dp0cx.bat" cx_09_evalcard.py || goto fail
echo ==================== 9/10  E6 certification budget (SP-temporal CAL) ====================
call "%~dp0cx.bat" cx_10_certbudget.py || goto fail
echo ==================== 10/10  figures ====================
call "%~dp0cx.bat" cx_08_figures.py || goto fail
echo.
echo FULL PIPELINE FINISHED - results in results\cleanx\
exit /b 0
:fail
echo.
echo FAILED in the step printed above. Stop here and send the output.
exit /b 1
