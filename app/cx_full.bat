@echo off
setlocal
cd /d "%~dp0"
echo ==================== 1/7  ESM-1b equivalence test + embeddings (resumable, 2-4 h) ====================
call "%~dp0cx.bat" cx_03_embed.py || goto fail
echo ==================== 2/7  E2 inference geometry, split100 (primary endpoint) ====================
call "%~dp0cx.bat" cx_05_geometry.py || goto fail
echo ==================== 3/7  E2 robustness: split70 backbone ====================
call "%~dp0cx.bat" cx_05_geometry.py --backbone split70 || goto fail
echo ==================== 4/7  E3 NISE biology ====================
call "%~dp0cx.bat" cx_06_ec_biology.py || goto fail
echo ==================== 5/7  E4 reliability ====================
call "%~dp0cx.bat" cx_07_reliability.py || goto fail
echo ==================== 6/7  E1 audit incl. temporal test ====================
call "%~dp0cx.bat" cx_01_audit.py --datasets new392 price149 sp_temporal_test || goto fail
echo ==================== 7/7  figures ====================
call "%~dp0cx.bat" cx_08_figures.py || goto fail
echo.
echo FULL PIPELINE FINISHED - results in results\cleanx\
exit /b 0
:fail
echo.
echo FAILED in the step printed above. Stop here and send the output.
exit /b 1
