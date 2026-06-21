#!/usr/bin/env bash
# Быстрые Lua integration tests.
# Один и тот же Lua-код гоняется в двух режимах:
#   library: luajit + build/lib/libquantum_analyzer_core.so
#   app:     bin/quantum_analyzer + встроенный static core через ffi.C

set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
TOTAL_FAILED=0

export PROJECT_ROOT

echo "============================================"
echo "  Quantum Analyzer — Lua Integration Tests"
echo "============================================"
echo ""

run_cmd() {
    local name="$1"
    shift

    echo "--- $name ---"
    "$@"
    local code=$?
    if [ "$code" -eq 0 ]; then
        echo "  -> PASSED"
    else
        echo "  -> FAILED ($code)"
        TOTAL_FAILED=$((TOTAL_FAILED + 1))
    fi
    echo ""
}

run_lua_test() {
    local mode="$1"
    local test_file="$2"
    local test_name="${mode}: $(basename "$test_file")"

    if [ "$mode" = "library" ]; then
        QA_TEST_MODE=library run_cmd "$test_name" luajit "$test_file"
    else
        QA_TEST_MODE=app run_cmd "$test_name" "$PROJECT_ROOT/bin/quantum_analyzer" "$test_file"
    fi
}

run_cmd "deps" xmake run deps
run_cmd "build static core" xmake --build lib
run_cmd "build shared core" xmake --build lib_shared
run_cmd "build app" xmake --build quantum_analyzer

# Проверка установки rockspec через LuaRocks в изолированное дерево.
run_cmd "luarocks install + smoke" bash "$SCRIPT_DIR/run_luarocks.sh"

echo "=============================="
echo "  LIBRARY MODE"
echo "=============================="
for f in "$SCRIPT_DIR"/lib/test_*.lua "$SCRIPT_DIR"/app/test_sandbox.lua "$SCRIPT_DIR"/app/test_viz_gnuplot.lua; do
    [ -f "$f" ] || continue
    run_lua_test library "$f"
done

echo "=============================="
echo "  APP / EMBEDDED MODE"
echo "=============================="
for f in "$SCRIPT_DIR"/lib/test_*.lua "$SCRIPT_DIR"/app/test_sandbox.lua "$SCRIPT_DIR"/app/test_viz_gnuplot.lua "$SCRIPT_DIR"/app/test_router.lua; do
    [ -f "$f" ] || continue
    run_lua_test app "$f"
done

echo "============================================"
if [ "$TOTAL_FAILED" -eq 0 ]; then
    echo "  ALL TESTS PASSED"
else
    echo "  $TOTAL_FAILED TEST(S) FAILED"
    exit 1
fi
echo "============================================"
exit 0
