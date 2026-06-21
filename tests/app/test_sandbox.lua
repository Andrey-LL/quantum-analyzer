local t = dofile((os.getenv("PROJECT_ROOT") or ".") .. "/tests/support.lua")

local sandbox = require("sandbox")

t.run_suite("sandbox", {
    {"module exports", function()
        t.assert_eq(type(sandbox.run), "function", "run")
        t.assert_eq(type(sandbox.execute), "function", "execute")
        t.assert_eq(type(sandbox.create_base_env), "function", "create_base_env")
        t.assert_eq(type(sandbox.group_molecules), "function", "group_molecules")
        t.assert_eq(type(sandbox.matrix_to_md), "function", "matrix_to_md")
        t.assert_eq(type(sandbox.table_to_md), "function", "table_to_md")
    end},

    {"base environment", function()
        local env = sandbox.create_base_env()
        t.assert_true(env.chemistry ~= nil, "chemistry in env")
        t.assert_true(env.Molecule ~= nil, "Molecule in env")
        t.assert_eq(type(env.emit), "function", "emit")
        t.assert_eq(type(env.h1), "function", "h1")
        t.assert_eq(type(env.matrix_to_md), "function", "matrix_to_md")
        t.assert_eq(type(env.atom_table), "function", "atom_table")
        t.assert_eq(type(env.bond_table), "function", "bond_table")
        t.assert_eq(type(env.Table), "function", "Table")
        t.assert_eq(type(env.col), "function", "col")
        t.assert_eq(type(env.val), "function", "val")
        t.assert_eq(type(env.summary), "function", "summary")
        t.assert_eq(type(env.delta), "function", "delta")
        t.assert_eq(type(env.cached), "function", "cached")
    end},

    {"safe environment hides unsafe globals", function()
        local env = sandbox.create_base_env()
        t.assert_eq(env.os, nil, "os hidden")
        t.assert_eq(env.io, nil, "io hidden")
        t.assert_eq(env.debug, nil, "debug hidden")
        t.assert_eq(env.require, nil, "require hidden")
        t.assert_eq(getmetatable(env), nil, "safe env has no global fallback")

        local outputs = sandbox.execute({
            [[assert(os == nil); assert(require == nil); print("safe")]],
            [[os.execute("true")]],
        }, nil)
        t.assert_true(outputs[1]:find("safe", 1, true) ~= nil, "safe code runs")
        t.assert_true(outputs[1]:find("Ошибка выполнения", 1, true) ~= nil, "os.execute is unavailable")
    end},

    {"unsafe environment exposes Lua globals", function()
        _G.__qa_sandbox_probe = "unsafe-fallback"
        local env = sandbox.create_base_env({unsafe = true})
        t.assert_eq(type(env.os.execute), "function", "os.execute available")
        t.assert_eq(type(env.io.open), "function", "io.open available")
        t.assert_eq(type(env.debug.traceback), "function", "debug available")
        t.assert_eq(type(env.require), "function", "require available")
        t.assert_eq(env.__qa_sandbox_probe, "unsafe-fallback", "global fallback")
        _G.__qa_sandbox_probe = nil
    end},

    {"plugin exports are explicit in safe environment", function()
        local env = sandbox.create_base_env()
        t.assert_eq(env.os, nil, "plugin does not expose os")
    end},

    {"execute without files keeps compatibility", function()
        local results = sandbox.execute({
            [[print("hello")]],
            [[emit({headers={"A","B"}, rows={{1,2}}})]],
        }, nil)

        t.assert_eq(#results, 1, "one document without files")
        t.assert_true(results[1]:find("hello", 1, true) ~= nil, "print captured")
        t.assert_true(results[1]:find("| A | B |", 1, true) ~= nil, "table rendered")
    end},

    {"single mode runs all blocks per file", function()
        local run = sandbox.run({
            [[
                h1("Single")
                print(mol:get_num_atoms())
            ]],
            [[metric("basis_size", "Размер базиса", basis_size)]],
        }, {t.fixture}, {mode = "single"})

        t.assert_eq(run.mode, "single", "single mode")
        t.assert_eq(#run.records, 1, "one record")
        t.assert_eq(#run.outputs, 1, "one output")
        t.assert_true(run.outputs[1]:find("Single", 1, true) ~= nil, "heading rendered")
        t.assert_true(run.outputs[1]:find("5", 1, true) ~= nil, "molecule available")
        t.assert_eq(run.records[1].molecule_key, "1:4;6:1", "composition key")
    end},

    {"grouped atom table", function()
        local run = sandbox.run({
            [[
                h1("Grouped")
                atom_table("atoms", {title = "Атомные номера", value_name = "Z"})
                for A = 1, mol:get_num_atoms() do
                    atom_value("atoms", A, mol:get_atomic_number(A))
                end
            ]],
        }, {t.fixture, t.project_root .. "/examples/fixtures/methane_sto-3g.log"}, {mode = "auto"})

        t.assert_eq(run.mode, "grouped", "auto grouped mode")
        t.assert_eq(#run.outputs, 1, "same molecule grouped into one output")
        t.assert_true(run.outputs[1]:find("Атомные номера", 1, true) ~= nil, "atom table title")
        t.assert_true(run.outputs[1]:find("| Атом | Z |", 1, true) ~= nil, "grouped table")
    end},

    {"pivot Table DSL groups basis columns", function()
        local run = sandbox.run({
            [[
                Table("DSL flat", {
                    rows = atoms,
                    left = {
                        col("Атом", atom_label),
                    },
                    across = basis,
                    values = {
                        val("Z", function(mol, atom) return atom.Z end, {fmt = "%.0f"}),
                    },
                    derived = {
                        summary("Сумма", {Z = "sum"}),
                    },
                    header = "flat",
                })

                Table("DSL spanned", {
                    rows = atoms,
                    left = {
                        col("Атом", atom_label),
                    },
                    across = basis,
                    values = {
                        val("Номер", function(mol, atom) return atom.id end, {fmt = "%.0f"}),
                    },
                    header = {"metric", "basis"},
                })

                Table("DSL derived", {
                    rows = atoms,
                    left = {
                        col("Атом", atom_label),
                    },
                    across = basis,
                    values = {
                        val("A", function(mol, atom) return atom.id end, {fmt = "%.0f"}),
                        val("B", function(mol, atom) return atom.id + 1 end, {fmt = "%.0f"}),
                    },
                    derived = {
                        delta("Средняя дельта", "A", "B", {fmt = "%.1f", mode = "mean_abs"}),
                    },
                    header = {"basis", "metric"},
                })

                Table("DSL transposed", {
                    rows = atoms,
                    left = {
                        col("Атом", atom_label),
                    },
                    across = basis,
                    values = {
                        val("Z", function(mol, atom) return atom.Z end, {fmt = "%.0f"}),
                    },
                    header = "flat",
                    transpose = true,
                })
            ]],
        }, {t.fixture, t.project_root .. "/examples/fixtures/methane_sto-3g.log"}, {mode = "auto"})

        t.assert_eq(run.mode, "grouped", "auto grouped mode")
        t.assert_eq(#run.outputs, 1, "same molecule grouped into one output")
        t.assert_true(run.outputs[1]:find("Z (6-31G)", 1, true) ~= nil, "flat basis column")
        t.assert_true(run.outputs[1]:find("Z (STO-3G)", 1, true) ~= nil, "flat second basis column")
        t.assert_true(run.outputs[1]:find("| Сумма | 10 | 10 |", 1, true) ~= nil, "summary row")
        local has_spanned_1 = run.outputs[1]:find("|  | 6-31G | STO-3G |", 1, true) ~= nil
        local has_spanned_2 = run.outputs[1]:find("|  | STO-3G | 6-31G |", 1, true) ~= nil
        t.assert_true(has_spanned_1 or has_spanned_2, "spanned basis row")
        t.assert_true(run.outputs[1]:find("| Средняя дельта | 1.0 | — | 1.0 | — |", 1, true) ~= nil, "derived delta row")
        t.assert_true(run.outputs[1]:find("## DSL transposed", 1, true) ~= nil, "transposed table title")
    end},

    {"different molecule keys produce different groups", function()
        local water = t.project_root .. "/examples/fixtures/water_sto-3g.log"
        local run = sandbox.run({
            [[
                atom_table("atoms", {title = "Атомы", value_name = "Z"})
                for A = 1, mol:get_num_atoms() do atom_value("atoms", A, mol:get_atomic_number(A)) end
            ]],
        }, {t.fixture, water}, {mode = "grouped"})

        t.assert_true(#run.outputs >= 2, "two molecule outputs")
        t.assert_eq(#run.errors, 0, "no file errors")
    end},

    {"plugin figures disable grouped mode", function()
        if t.mode ~= "app" then return end

        local out = t.project_root .. "/.cache/test_outputs/sandbox/no_group_bar.svg"
        os.remove(out)

        local run = sandbox.run({
            [[
                Table("Grouped table", {
                    rows = atoms(),
                    left = {col("Atom", atom_label)},
                    values = {val("Z", function(_, atom) return atom.Z end, {fmt = "%.0f"})},
                    header = "flat",
                })

                bar({
                    {label = basis_name, value = basis_size},
                    {label = "zero", value = 0},
                }, {
                    caption = "Basis size",
                    output = "]] .. out .. [[",
                    link = "no_group_bar.svg",
                })
            ]],
        }, {t.fixture, t.project_root .. "/examples/fixtures/methane_sto-3g.log"}, {mode = "auto"})

        t.assert_eq(run.mode, "single", "figure output disables grouping")
        t.assert_eq(#run.outputs, 2, "one output per input file")
        t.assert_true(run.outputs[1]:find("![Basis size](no_group_bar.svg)", 1, true) ~= nil, "figure markdown")
    end},

    {"broken file does not abort batch", function()
        local run = sandbox.run({
            [[text("ok")]],
        }, {t.fixture, t.project_root .. "/examples/does-not-exist.log"}, {mode = "single"})

        t.assert_eq(#run.records, 2, "success and error records")
        t.assert_eq(#run.errors, 1, "one error")
        t.assert_true(#run.outputs >= 2, "success output plus error output")
        t.assert_true(run.outputs[#run.outputs]:find("Ошибки обработки", 1, true) ~= nil, "error section")
    end},

    {"group molecules smoke", function()
        local groups = sandbox.group_molecules({t.fixture})
        t.assert_eq(type(groups), "table", "groups result")
        t.assert_true(groups["1:4;6:1"] ~= nil, "composition key present")
    end},
})
