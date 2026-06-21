#!/usr/bin/env luajit
-- Запуск отдельного интеграционного Lua-теста.

local test_file = arg[1]
if not test_file then
    print("Usage: luajit tests/runner.lua <test_file>")
    os.exit(1)
end

-- Алиасы позволяют запускать тесты из структуры LuaRocks и из дерева исходников.
local ok, mod = pcall(require, "quantum_analyzer.chemistry_ffi")
if ok then package.loaded["chemistry_ffi"] = mod end

local ok, mod = pcall(require, "quantum_analyzer.quantum_analysis")
if ok then package.loaded["quantum_analysis"] = mod end

local ok, mod = pcall(require, "quantum_analyzer.sandbox")
if ok then package.loaded["sandbox"] = mod end

-- PROJECT_ROOT передаётся тесту через переменную окружения как абсолютный путь.
local script_dir = arg[0]:match("(.*/)") or "./"
local project_root = script_dir:gsub("tests/$", "")
if project_root:sub(1, 1) ~= "/" then
    project_root = os.getenv("PWD") and (os.getenv("PWD") .. "/" .. project_root:gsub("^%./", "")) or project_root
end
os.setenv("PROJECT_ROOT", project_root)

local ok, err = pcall(dofile, test_file)
if not ok then
    print("TEST ERROR: " .. tostring(err))
    os.exit(1)
end
