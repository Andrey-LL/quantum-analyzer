-- Общие помощники для интеграционных Lua-тестов.

local M = {}

local function dirname(path)
    return (path:gsub("\\", "/"):match("^(.*)/[^/]*$")) or "."
end

local function detect_project_root()
    local root = os.getenv("PROJECT_ROOT")
    if root and root ~= "" then
        return root:gsub("/+$", "")
    end

    local source = debug.getinfo(1, "S").source:gsub("^@", "")
    local tests_dir = dirname(source)
    return dirname(tests_dir)
end

M.project_root = detect_project_root()
M.mode = os.getenv("QA_TEST_MODE") or "library"
M.fixture = M.project_root .. "/examples/fixtures/methane_6-31g.log"

package.path = table.concat({
    M.project_root .. "/tests/?.lua",
    M.project_root .. "/src/lib/lua_core/?.lua",
    M.project_root .. "/src/app/?.lua",
    M.project_root .. "/src/app/share/?/?.lua",
    package.path,
}, ";")

package.cpath = table.concat({
    M.project_root .. "/build/lib/?.so",
    M.project_root .. "/build/lib/?.dll",
    M.project_root .. "/build/lib/?.dylib",
    package.cpath,
}, ";")

function M.assert_true(value, message)
    if not value then
        error(message or "assertion failed", 2)
    end
end

function M.assert_eq(actual, expected, message)
    if actual ~= expected then
        error(string.format("%s: expected %s, got %s", message or "values differ", tostring(expected), tostring(actual)), 2)
    end
end

function M.assert_near(actual, expected, eps, message)
    if math.abs(actual - expected) > eps then
        error(string.format("%s: expected %.12g +/- %.3g, got %.12g", message or "values differ", expected, eps, actual), 2)
    end
end

function M.test(name, fn)
    io.write("  - " .. name .. " ... ")
    local ok, err = pcall(fn)
    if ok then
        print("ok")
        return true
    end
    print("FAILED")
    print("    " .. tostring(err))
    return false
end

function M.run_suite(name, tests)
    print(string.format("[%s] %s", M.mode, name))
    local failed = 0
    for _, case in ipairs(tests) do
        if not M.test(case[1], case[2]) then
            failed = failed + 1
        end
    end
    if failed > 0 then
        error(string.format("%s: %d failed", name, failed), 0)
    end
end

return M
