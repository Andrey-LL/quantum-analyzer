local t = dofile((os.getenv("PROJECT_ROOT") or ".") .. "/tests/support.lua")

local qa = require("quantum_analysis")
local chemistry = require("chemistry_ffi")

local function with_molecule(fn)
    local mol = assert(chemistry.Molecule.load(t.fixture))
    local density = mol:density()
    local overlap = mol:overlap()
    local ok, err = pcall(fn, mol, density, overlap)
    density:free()
    overlap:free()
    mol:close()
    if not ok then error(err, 0) end
end

t.run_suite("quantum_analysis", {
    {"module exports", function()
        t.assert_eq(type(qa.compute_PS), "function", "compute_PS")
        t.assert_eq(type(qa.mulliken_populations), "function", "mulliken_populations")
        t.assert_eq(type(qa.mulliken_charges), "function", "mulliken_charges")
        t.assert_eq(type(qa.wiberg_bond_orders), "function", "wiberg_bond_orders")
        t.assert_eq(type(qa.mayer_bond_orders), "function", "mayer_bond_orders")
    end},

    {"P*S and electron count", function()
        with_molecule(function(_, density, overlap)
            local PS = qa.compute_PS(density, overlap)
            local rows, cols = PS:size()
            t.assert_eq(rows, cols, "PS is square")
            t.assert_near(PS:trace(), 10.0, 1e-3, "PS trace")
            PS:free()
        end)
    end},

    {"Mulliken populations and charges", function()
        with_molecule(function(mol, density, overlap)
            local natoms = mol:get_num_atoms()
            local pops = qa.mulliken_populations(mol, density, overlap)
            local charges = qa.mulliken_charges(mol, density, overlap)

            t.assert_eq(#pops, natoms, "population count")
            t.assert_eq(#charges, natoms, "charge count")

            local pop_sum = 0
            for _, value in ipairs(pops) do pop_sum = pop_sum + value end
            t.assert_near(pop_sum, 10.0, 1e-3, "population sum")

            local charge_sum = 0
            for _, value in ipairs(charges) do charge_sum = charge_sum + value end
            t.assert_near(charge_sum, 0.0, 1e-3, "neutral molecule charge")
        end)
    end},

    {"bond order tables", function()
        with_molecule(function(mol, density, overlap)
            local natoms = mol:get_num_atoms()
            local wiberg = qa.wiberg_bond_orders(mol, density, overlap)
            local mayer = qa.mayer_bond_orders(mol, density, overlap)

            t.assert_eq(#wiberg, natoms, "Wiberg row count")
            t.assert_eq(#mayer, natoms, "Mayer row count")
            t.assert_true(type(wiberg[1][2]) == "number", "Wiberg C-H value")
            t.assert_true(type(mayer[1][2]) == "number", "Mayer C-H value")
        end)
    end},

    {"advanced analysis smoke", function()
        with_molecule(function(mol, density, overlap)
            local lowdin = qa.loewdin_charges(mol, density, overlap)
            local valences = qa.atomic_valence_indices(mol, density, overlap)
            local PS = qa.compute_PS(density, overlap)
            local ao2atom = qa.get_ao_atom_mapping(mol)
            local lone_pairs = qa.analyze_lone_pairs(mol, PS, ao2atom)

            t.assert_eq(#lowdin, mol:get_num_atoms(), "Lowdin charge count")
            t.assert_eq(#valences, mol:get_num_atoms(), "valence count")
            t.assert_eq(type(lone_pairs), "table", "lone pair result")
            PS:free()
        end)
    end},
})

do return end

local project_root = os.getenv("PROJECT_ROOT") or "./"

print("=== quantum_analysis Tests ===\n")

local qa = require("quantum_analysis")
local chemistry = require("chemistry_ffi")

local mol_ok, mol = pcall(chemistry.Molecule.load, project_root .. "examples/fixtures/methane_6-31g.log")
assert(mol_ok, "methane_6-31g.log загружается для анализа")

if not mol_ok then
    print("  Не могу загрузить молекулу — остальные тесты пропущены")
    print(string.format("\n=== Результат: %d/%d пройдено, %d провалено ===", passed, total, failed))
    os.exit(failed > 0 and 1 or 0)
end

print("[1] Базовые свойства")
local natoms = mol:get_num_atoms()
assert(natoms == 5, string.format("Метан: 5 атомов (получено: %d)", natoms))
local nbasis = mol:get_basis_size()
assert(nbasis > 0, string.format("Базис > 0 (%d)", nbasis))

print("\n[2] PS матрица")
local density = mol:density()
local overlap = mol:overlap()
local PS = qa.compute_PS(density, overlap)
assert(PS ~= nil, "PS = P*S вычисляется")
if PS then local r, c = PS:size(); assert(r == c, string.format("PS квадратная (%dx%d)", r, c)); PS:free() end

print("\n[3] Mulliken популяции")
local pops = qa.mulliken_populations(mol, density, overlap)
assert(pops ~= nil, "Mulliken populations возвращает результат")
assert(#pops == natoms, string.format("Популяции для всех атомов (%d)", #pops))
local total_pop = 0
for _, p in ipairs(pops) do total_pop = total_pop + p end
assert(total_pop > 0, string.format("Суммарная популяция > 0 (%.2f)", total_pop))

print("\n[4] Mulliken заряды")
local charges = qa.mulliken_charges(mol, density, overlap)
assert(charges ~= nil, "Mulliken charges возвращает результат")
assert(#charges == natoms, string.format("Заряды для всех атомов (%d)", #charges))
local total_charge = 0
for _, q in ipairs(charges) do total_charge = total_charge + q end
assert(math.abs(total_charge) < 0.5, string.format("Сумма зарядов ≈ 0 (%.4f)", total_charge))

print("\n[5] Löwdin заряды")
local lowdin = qa.loewdin_charges(mol, density, overlap)
assert(lowdin ~= nil, "Löwdin charges возвращает результат")

print("\n[6] Wiberg bond orders")
local wiberg = qa.wiberg_bond_orders(mol, density, overlap)
assert(wiberg ~= nil, "Wiberg bond orders вычисляются")

print("\n[7] Mayer bond orders")
local mayer = qa.mayer_bond_orders(mol, density, overlap)
assert(mayer ~= nil, "Mayer bond orders вычисляются")

print("\n[8] Atomic valence indices")
local valences = qa.atomic_valence_indices(mol, density, overlap)
assert(valences ~= nil, "Valence indices вычисляются")

print("\n[9] Lone pair analysis")
local lp = qa.analyze_lone_pairs(mol)
assert(lp ~= nil, "Lone pair analysis возвращает результат")

density:free()
overlap:free()
mol:close()

print(string.format("\n=== Результат: %d/%d пройдено, %d провалено ===", passed, total, failed))
os.exit(failed > 0 and 1 or 0)
