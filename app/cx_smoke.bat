@echo off
setlocal
cd /d "%~dp0"
echo ==================== 1/13  environment + exact CLEAN reproduction ====================
call "%~dp0cx.bat" cx_00_check.py || goto fail
echo ==================== 2/13  E1 evaluation audit (smoke) ====================
call "%~dp0cx.bat" cx_01_audit.py --smoke || goto fail
echo ==================== 3/13  temporal benchmark builder (smoke, cached downloads) ====================
call "%~dp0cx.bat" cx_02_build_temporal.py --smoke || goto fail
echo ==================== 4/13  identity search (5 proteins) ====================
call "%~dp0cx.bat" cx_04_identity.py --sets new392 --limit 5 || goto fail
echo ==================== 5/13  E2 inference geometry (smoke) ====================
call "%~dp0cx.bat" cx_05_geometry.py --smoke || goto fail
echo ==================== 6/13  E3 NISE biology (smoke) ====================
call "%~dp0cx.bat" cx_06_ec_biology.py --smoke || goto fail
echo ==================== 7/13  E4 reliability (smoke) ====================
call "%~dp0cx.bat" cx_07_reliability.py --smoke || goto fail
echo ==================== 8/13  unit + integration tests (synthetic data, temp folders only) ====================
call "%~dp0cx.bat" tests\run_tests.py || goto fail
echo ==================== 9/13  E5 evaluation card + EC-Bench re-scoring (smoke) ====================
call "%~dp0cx.bat" cx_09_evalcard.py --smoke || goto fail
echo ==================== 10/13  E6 certification budget (smoke) ====================
call "%~dp0cx.bat" cx_10_certbudget.py --smoke || goto fail
echo ==================== 11/13  BLASTp tables match the local sequences (+ re-search if BLAST+ is installed) ====================
call "%~dp0cx.bat" cx_11_blast.py --smoke || goto fail
echo ==================== 12/13  E7 UniProt evidence parser (known answers, no network) ====================
call "%~dp0cx.bat" cx_12_uniprot_evidence.py --smoke || goto fail
echo ==================== 13/13  E7 label-provenance analysis (synthetic data, known answers) ====================
call "%~dp0cx.bat" cx_13_provenance.py --smoke || goto fail
echo.
echo ALL 13 SMOKE TESTS PASSED
exit /b 0
:fail
echo.
echo SMOKE TEST FAILED in the step printed above. Stop here and send the output.
exit /b 1
