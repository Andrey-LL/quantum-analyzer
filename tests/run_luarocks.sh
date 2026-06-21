#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
TREE="$PROJECT_ROOT/.cache/luarocks-test-tree"
ROCKSPEC="$PROJECT_ROOT/quantum_analyzer-1.0-1.rockspec"

if ! command -v luarocks >/dev/null 2>&1; then
    echo "[luarocks] not found"
    exit 1
fi

if ! command -v luajit >/dev/null 2>&1; then
    echo "[luarocks] luajit not found"
    exit 1
fi

if [ ! -f "$ROCKSPEC" ]; then
    echo "[luarocks] rockspec not found: $ROCKSPEC"
    exit 1
fi

rm -rf "$TREE"
mkdir -p "$TREE"

luarocks --tree="$TREE" make "$ROCKSPEC"

ROCK_LIBDIR="$TREE/lib/lua/5.1"
ROCK_LUADIR="$TREE/share/lua/5.1"

if [ ! -f "$ROCK_LIBDIR/libquantum_analyzer_core.so" ] || [ ! -f "$ROCK_LUADIR/chemistry_ffi.lua" ] || [ ! -f "$ROCK_LUADIR/quantum_analysis.lua" ]; then
    echo "[luarocks] installed files not found"
    exit 1
fi

BASE_LUA_PATH="$PROJECT_ROOT/tests/?.lua;;"

# Используем стандартный способ от LuaRocks: экспорт путей через `luarocks path`.
# Это убирает ручное ковыряние LUA_PATH/LUA_CPATH.
eval "$(luarocks --tree="$TREE" path)"
export LUA_PATH="${LUA_PATH:-}${BASE_LUA_PATH}"

echo "--- luarocks: chemistry_ffi load/open/close ---"
luajit -e 'local c=require("chemistry_ffi"); local m=assert(c.Molecule.load("examples/fixtures/methane_6-31g.log")); print(m:get_basis_name()); m:close(); print("chemistry ok")'

echo "--- luarocks: quantum_analysis module ---"
luajit -e 'local qa=require("quantum_analysis"); assert(type(qa)=="table"); print("quantum_analysis ok")'

echo "--- luarocks: integrated PS workflow ---"
luajit -e 'local c=require("chemistry_ffi"); local qa=require("quantum_analysis"); local m=assert(c.Molecule.load("examples/fixtures/methane_6-31g.log")); local S=m:overlap(); local P=m:density(); local PS=qa.compute_PS(P,S); local tr=PS:trace(); assert(math.abs(tr-10.0)<1e-3); PS:free(); S:free(); P:free(); m:close(); print("workflow ok")'

echo "luarocks library checks ok"
