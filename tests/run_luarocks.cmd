@echo off
setlocal EnableExtensions

set "SCRIPT_DIR=%~dp0"
for %%I in ("%SCRIPT_DIR%..") do set "PROJECT_ROOT=%%~fI"
cd /d "%PROJECT_ROOT%"

where luarocks >nul 2>nul
if errorlevel 1 (
    echo [luarocks] not found
    exit /b 1
)

where luajit >nul 2>nul
if errorlevel 1 (
    echo [luarocks] luajit not found
    exit /b 1
)

set "TREE=%PROJECT_ROOT%\.cache\luarocks-test-tree-win"
set "ROCKSPEC=%PROJECT_ROOT%\quantum_analyzer-1.0-1.rockspec"

if not exist "%ROCKSPEC%" (
    echo [luarocks] rockspec not found: %ROCKSPEC%
    exit /b 1
)

if exist "%TREE%" rmdir /s /q "%TREE%"
mkdir "%TREE%"

luarocks --tree="%TREE%" make "%ROCKSPEC%"
if errorlevel 1 (
    echo [luarocks] make failed
    exit /b 1
)

set "ROCK_LIBDIR=%TREE%\lib\lua\5.1"
set "ROCK_LUADIR=%TREE%\share\lua\5.1"

if not exist "%ROCK_LIBDIR%\quantum_analyzer_core.dll" if not exist "%ROCK_LIBDIR%\libquantum_analyzer_core.dll" (
    echo [luarocks] installed native library not found
    exit /b 1
)
if not exist "%ROCK_LUADIR%\chemistry_ffi.lua" (
    echo [luarocks] chemistry_ffi.lua not found
    exit /b 1
)
if not exist "%ROCK_LUADIR%\quantum_analysis.lua" (
    echo [luarocks] quantum_analysis.lua not found
    exit /b 1
)

for /f "usebackq delims=" %%A in (`luarocks --tree="%TREE%" path --lr-path`) do set "LUA_PATH=%%A"
for /f "usebackq delims=" %%A in (`luarocks --tree="%TREE%" path --lr-cpath`) do set "LUA_CPATH=%%A"
set "LUA_PATH=%LUA_PATH%;%PROJECT_ROOT%\tests\?.lua;;"

echo --- luarocks: chemistry_ffi load/open/close ---
luajit -e "local c=require('chemistry_ffi'); local m=assert(c.Molecule.load('examples/fixtures/methane_6-31g.log')); print(m:get_basis_name()); m:close(); print('chemistry ok')"
if errorlevel 1 exit /b 1

echo --- luarocks: quantum_analysis module ---
luajit -e "local qa=require('quantum_analysis'); assert(type(qa)=='table'); print('quantum_analysis ok')"
if errorlevel 1 exit /b 1

echo --- luarocks: integrated PS workflow ---
luajit -e "local c=require('chemistry_ffi'); local qa=require('quantum_analysis'); local m=assert(c.Molecule.load('examples/fixtures/methane_6-31g.log')); local S=m:overlap(); local P=m:density(); local PS=qa.compute_PS(P,S); local tr=PS:trace(); assert(math.abs(tr-10.0)<1e-3); PS:free(); S:free(); P:free(); m:close(); print('workflow ok')"
if errorlevel 1 exit /b 1

echo luarocks library checks ok
exit /b 0
