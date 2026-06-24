local t = dofile((os.getenv("PROJECT_ROOT") or ".") .. "/tests/support.lua")
local is_windows = os.getenv("QA_TEST_WINDOWS") == "1"
    or package.config:sub(1, 1) == "\\"
    or (os.getenv("OS") or ""):lower():find("windows", 1, true) ~= nil
    or os.getenv("MSYSTEM") ~= nil
local app_bin = t.project_root .. (is_windows and "/bin/quantum_analyzer.exe" or "/bin/quantum_analyzer")

local function shell_quote(s)
    if is_windows then
        return '"' .. tostring(s):gsub('"', '\\"') .. '"'
    end
    return "'" .. s:gsub("'", "'\\''") .. "'"
end

local function capture(command)
    local pipe = assert(io.popen(command .. " 2>&1"))
    local output = pipe:read("*a")
    local ok, _, code = pipe:close()
    return output, ok, code
end

local function remove_tree(path)
    if is_windows then
        os.execute('if exist ' .. shell_quote(path) .. ' rmdir /s /q ' .. shell_quote(path))
    else
        os.execute("rm -rf " .. shell_quote(path))
    end
end

t.run_suite("binary router", {
    {"binary exists", function()
        local f = io.open(app_bin, "rb")
        t.assert_true(f ~= nil, "quantum_analyzer binary exists")
        if f then f:close() end
    end},

    {"-e executes code", function()
        local cmd = shell_quote(app_bin) .. " -e " .. shell_quote([[print('qa-app-ok')]])
        local output, ok = capture(cmd)
        t.assert_true(ok, "command exits successfully")
        t.assert_true(output:find("qa%-app%-ok") ~= nil, "prints marker")
    end},

    {"embedded ffi can open fixture", function()
        local fixture = t.fixture:gsub("\\", "/")
        local code = "local chemistry = require('chemistry_ffi')\n" ..
            "local mol = chemistry.Molecule.load([[" .. fixture .. "]])\n" ..
            "print(mol:get_num_atoms())\n" ..
            "mol:close()"
        local cmd = shell_quote(app_bin) .. " -e " .. shell_quote(code)
        local output, ok = capture(cmd)
        t.assert_true(ok, "command exits successfully")
        t.assert_true(output:find("5", 1, true) ~= nil, "fixture parsed")
    end},

    {"batch routes through sandbox", function()
        local cmd = shell_quote(app_bin) ..
            " --batch --template standard_analysis --files " .. shell_quote(t.fixture)
        local output, ok = capture(cmd)
        t.assert_true(ok, "batch exits successfully")
        t.assert_true(output:find("Completed: 1 successful, 0 failed", 1, true) ~= nil, "batch summary")
        t.assert_true(output:find("Детальный блочный анализ", 1, true) ~= nil, "template output")
    end},

    {"batch accepts multiple files", function()
        local cmd = shell_quote(app_bin) ..
            " --batch --template standard_analysis --files " ..
            shell_quote(t.fixture) .. " " .. shell_quote(t.project_root .. "/examples/fixtures/methane_sto-3g.log")
        local output, ok = capture(cmd)
        t.assert_true(ok, "batch exits successfully")
        t.assert_true(output:find("Completed: 2 successful, 0 failed", 1, true) ~= nil, "two files processed")
    end},

    {"plain-code summary template runs through sandbox", function()
        local cmd = shell_quote(app_bin) ..
            " --batch --template summary_analysis --files " ..
            shell_quote(t.fixture) .. " " .. shell_quote(t.project_root .. "/examples/fixtures/methane_sto-3g.log")
        local output, ok = capture(cmd)
        t.assert_true(ok, "summary batch exits successfully")
        t.assert_true(output:find("Completed: 2 successful, 0 failed", 1, true) ~= nil, "summary processed")
        t.assert_true(output:find("Заряды и валентность", 1, true) ~= nil, "summary table")
        t.assert_true(output:find("Индексы связей", 1, true) ~= nil, "bond table")
    end},

    {"single molecule overview template runs", function()
        local cmd = shell_quote(app_bin) ..
            " --batch --template single_molecule_overview --files " .. shell_quote(t.fixture)
        local output, ok = capture(cmd)
        t.assert_true(ok, "single molecule overview exits successfully")
        t.assert_true(output:find("Completed: 1 successful, 0 failed", 1, true) ~= nil, "overview processed")
        t.assert_true(output:find("Краткий анализ молекулы", 1, true) ~= nil, "overview heading")
        t.assert_true(output:find("Атомная сводка", 1, true) ~= nil, "overview atom table")
    end},

    {"method comparison template groups basis reports", function()
        local cmd = shell_quote(app_bin) ..
            " --batch --template method_comparison_analysis --files " ..
            shell_quote(t.fixture) .. " " .. shell_quote(t.project_root .. "/examples/fixtures/methane_sto-3g.log")
        local output, ok = capture(cmd)
        t.assert_true(ok, "method comparison exits successfully")
        t.assert_true(output:find("Completed: 2 successful, 0 failed", 1, true) ~= nil, "comparison processed")
        t.assert_true(output:find("Сравнение методов зарядов", 1, true) ~= nil, "comparison heading")
        t.assert_true(output:find("Сводка расхождения методов", 1, true) ~= nil, "comparison summary")
    end},

    {"summary groups two molecules", function()
        local cmd = shell_quote(app_bin) ..
            " --batch --template summary_analysis --files " ..
            shell_quote(t.fixture) .. " " ..
            shell_quote(t.project_root .. "/examples/fixtures/methane_sto-3g.log") .. " " ..
            shell_quote(t.project_root .. "/examples/fixtures/water_sto-3g.log")
        local output, ok = capture(cmd)
        t.assert_true(ok, "summary grouped command exits successfully")
        t.assert_true(output:find("Completed: 3 successful, 0 failed", 1, true) ~= nil, "three files processed")
        t.assert_true(output:find("# Молекула 1:4;6:1", 1, true) ~= nil, "methane group output")
        t.assert_true(output:find("# Молекула 1:2;8:1", 1, true) ~= nil, "water group output")
    end},

    {"batch saves grouped reports to output directory", function()
        local out_dir = t.project_root .. "/.cache/test_outputs/router_reports"
        remove_tree(out_dir)

        local cmd = shell_quote(app_bin) ..
            " --batch --template summary_analysis --out-dir " .. shell_quote(out_dir) .. " --files " ..
            shell_quote(t.fixture) .. " " ..
            shell_quote(t.project_root .. "/examples/fixtures/methane_sto-3g.log") .. " " ..
            shell_quote(t.project_root .. "/examples/fixtures/water_sto-3g.log")
        local output, ok = capture(cmd)
        t.assert_true(ok, "batch exits successfully")
        t.assert_true(output:find("Saved report:", 1, true) ~= nil, "saved report message")

        local water = out_dir .. "/01_molecule_1_2_8_1.md"
        local methane = out_dir .. "/02_molecule_1_4_6_1.md"
        local wf = io.open(water, "rb")
        local mf = io.open(methane, "rb")
        t.assert_true(wf ~= nil, "water report saved")
        t.assert_true(mf ~= nil, "methane report saved")
        if wf then
            local data = wf:read("*a")
            wf:close()
            t.assert_true(data:find("# Молекула 1:2;8:1", 1, true) ~= nil, "water report content")
        end
        if mf then
            local data = mf:read("*a")
            mf:close()
            t.assert_true(data:find("# Молекула 1:4;6:1", 1, true) ~= nil, "methane report content")
        end
    end},

    {"help routes through router", function()
        local cmd = shell_quote(app_bin) .. " --help"
        local output, ok = capture(cmd)
        t.assert_true(ok, "help exits successfully")
        t.assert_true(output:find("Usage:", 1, true) ~= nil, "usage text")
        t.assert_true(output:find("--batch", 1, true) ~= nil, "batch option")
    end},

    {"no args prints standalone banner", function()
        local cmd = shell_quote(app_bin)
        local output, ok = capture(cmd)
        t.assert_true(ok, "command exits successfully")
        t.assert_true(output:find("Standalone", 1, true) ~= nil, "standalone banner")
    end},
})
