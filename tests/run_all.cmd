@echo off
setlocal EnableExtensions

set "SCRIPT_DIR=%~dp0"
for %%I in ("%SCRIPT_DIR%..") do set "PROJECT_ROOT=%%~fI"
cd /d "%PROJECT_ROOT%"

set "TOTAL_FAILED=0"
set "QA_TEST_MODE=app"

echo ============================================
echo   Quantum Analyzer - Windows Integration Tests
echo ============================================
echo.

call :run "build static core" "xmake --build lib"
call :run "build shared core" "xmake --build lib_shared"
call :run "build app" "xmake --build quantum_analyzer"

echo ==============================
echo   APP / EMBEDDED MODE
echo ==============================
call :run "app: test_chemistry_ffi.lua" ""%PROJECT_ROOT%\bin\quantum_analyzer.exe" "%PROJECT_ROOT%\tests\lib\test_chemistry_ffi.lua""
call :run "app: test_memory_lifecycle.lua" ""%PROJECT_ROOT%\bin\quantum_analyzer.exe" "%PROJECT_ROOT%\tests\lib\test_memory_lifecycle.lua""
call :run "app: test_quantum_analysis.lua" ""%PROJECT_ROOT%\bin\quantum_analyzer.exe" "%PROJECT_ROOT%\tests\lib\test_quantum_analysis.lua""
call :run "app: test_sandbox.lua" ""%PROJECT_ROOT%\bin\quantum_analyzer.exe" "%PROJECT_ROOT%\tests\app\test_sandbox.lua""

echo ==============================
echo   LUAROCKS / LIBRARY MODE
echo ==============================
where luarocks >nul 2>nul
if errorlevel 1 (
    if "%QA_TEST_LUAROCKS%"=="1" (
        echo [luarocks] missing but QA_TEST_LUAROCKS=1
        set /a TOTAL_FAILED+=1
    ) else (
        echo [skip] LuaRocks not found
    )
) else (
    call tests\run_luarocks.cmd
    if errorlevel 1 set /a TOTAL_FAILED+=1
)
echo.

echo ============================================
if "%TOTAL_FAILED%"=="0" (
    echo   ALL WINDOWS TESTS PASSED
) else (
    echo   %TOTAL_FAILED% TESTS FAILED
    exit /b 1
)
echo ============================================
exit /b 0

:run
set "NAME=%~1"
set "CMD=%~2"
echo --- %NAME% ---
%CMD%
if errorlevel 1 (
    echo   -^> FAILED
    set /a TOTAL_FAILED+=1
) else (
    echo   -^> PASSED
)
echo.
exit /b 0
