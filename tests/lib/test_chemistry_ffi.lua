local t = dofile((os.getenv("PROJECT_ROOT") or ".") .. "/tests/support.lua")

local chemistry = require("chemistry_ffi")

local function with_molecule(fn)
    local mol = assert(chemistry.Molecule.load(t.fixture))
    local ok, err = pcall(fn, mol)
    mol:close()
    if not ok then error(err, 0) end
end

t.run_suite("chemistry_ffi", {
    {"module exports", function()
        t.assert_eq(type(chemistry.Molecule), "table", "Molecule export")
        t.assert_eq(type(chemistry.Matrix), "table", "Matrix export")
        t.assert_eq(type(chemistry.Group), "table", "Group export")
        t.assert_true(chemistry.gaussian ~= nil, "FFI namespace")
    end},

    {"gaussian parser metadata", function()
        with_molecule(function(mol)
            t.assert_eq(mol:get_num_atoms(), 5, "methane atom count")
            t.assert_true(mol:get_basis_size() > 0, "basis size")
            t.assert_eq(mol:get_atomic_number(1), 6, "carbon atomic number")
            t.assert_eq(mol:get_atomic_number(2), 1, "hydrogen atomic number")
            t.assert_eq(type(mol:get_nuclear_repulsion()), "number", "nuclear repulsion")
            local mapping = mol:ao_atom_mapping()
            t.assert_eq(#mapping, mol:get_basis_size(), "AO atom mapping size")
        end)
    end},

    {"overlap and density matrices", function()
        with_molecule(function(mol)
            local S = mol:overlap()
            local P = mol:density()
            local sr, sc = S:size()
            local pr, pc = P:size()

            t.assert_eq(sr, sc, "overlap is square")
            t.assert_eq(pr, pc, "density is square")
            t.assert_eq(sr, mol:get_basis_size(), "overlap size matches basis")
            t.assert_true(S:is_symmetric(), "overlap is symmetric")
            t.assert_true(P:is_symmetric(), "density is symmetric")

            for i = 1, math.min(sr, 5) do
                t.assert_near(S:get(i, i), 1.0, 1e-6, "overlap diagonal")
            end

            local PS = P * S
            t.assert_near(PS:trace(), 10.0, 1e-3, "Tr(P*S) electron count")
            PS:free()
            S:free()
            P:free()
        end)
    end},

    {"matrix operations", function()
        with_molecule(function(mol)
            local S = mol:overlap()
            local P = mol:density()

            local sum = S + P
            local diff = S - P
            local scaled = S * 0.5
            local transposed = S:transpose()

            local rows = sum:size()
            t.assert_eq(rows, mol:get_basis_size(), "matrix add")
            rows = diff:size()
            t.assert_eq(rows, mol:get_basis_size(), "matrix subtract")
            rows = scaled:size()
            t.assert_eq(rows, mol:get_basis_size(), "matrix scale")
            rows = transposed:size()
            t.assert_eq(rows, mol:get_basis_size(), "matrix transpose")

            sum:free()
            diff:free()
            scaled:free()
            transposed:free()
            S:free()
            P:free()
        end)
    end},

    {"manual matrix primitives", function()
        with_molecule(function(mol)
            local A = chemistry.Matrix.from_table(mol, 2, 2, {2, 0, 0, 8})
            local B = chemistry.Matrix.from_table(mol, 2, 2, {1, 2, 0, 1})

            t.assert_true(A:is_symmetric(), "manual symmetric matrix")
            t.assert_true(not B:is_symmetric(), "manual non-symmetric matrix")
            t.assert_near(A:trace(), 10.0, 1e-12, "manual trace")
            t.assert_near(A:condition_number(), 4.0, 1e-12, "symmetric condition number")

            local eig = A:eigenvalues()
            t.assert_near(eig[1], 2.0, 1e-12, "first eigenvalue")
            t.assert_near(eig[2], 8.0, 1e-12, "second eigenvalue")

            local sqrt_a = A:pow(0.5)
            t.assert_near(sqrt_a:get(1, 1), math.sqrt(2), 1e-12, "matrix_power diagonal")
            t.assert_near(sqrt_a:get(2, 2), math.sqrt(8), 1e-12, "matrix_power diagonal")

            local sqrt_lapack = A:symm_pow(0.5)
            t.assert_near(sqrt_lapack:get(1, 1), math.sqrt(2), 1e-12, "matrix_symm_pow diagonal")
            t.assert_near(sqrt_lapack:get(2, 2), math.sqrt(8), 1e-12, "matrix_symm_pow diagonal")

            local ok = pcall(function() return B:pow(0.5) end)
            t.assert_true(not ok, "matrix_power rejects non-symmetric matrix")

            local had = A:hadamard(A)
            t.assert_near(had:get(1, 1), 4.0, 1e-12, "hadamard")
            t.assert_near(had:get(2, 2), 64.0, 1e-12, "hadamard")

            local cwise = A:cwise_pow(2)
            t.assert_near(cwise:get(1, 1), 4.0, 1e-12, "cwise pow")
            t.assert_near(cwise:get(2, 2), 64.0, 1e-12, "cwise pow")

            local C = chemistry.Matrix.from_table(mol, 2, 2, {0.1, -2.0, 3.0, 0.05})
            local threshold = C:threshold(1.0)
            t.assert_near(threshold:get(1, 1), 0.0, 1e-12, "threshold")
            t.assert_near(threshold:get(1, 2), -2.0, 1e-12, "threshold")
            t.assert_near(threshold:get(2, 1), 3.0, 1e-12, "threshold")
            t.assert_near(threshold:get(2, 2), 0.0, 1e-12, "threshold")

            local clamp = C:clamp(-1.0, 1.0)
            t.assert_near(clamp:get(1, 1), 0.1, 1e-12, "clamp")
            t.assert_near(clamp:get(1, 2), -1.0, 1e-12, "clamp")
            t.assert_near(clamp:get(2, 1), 1.0, 1e-12, "clamp")
            t.assert_near(clamp:get(2, 2), 0.05, 1e-12, "clamp")

            t.assert_near(A:norm_fro(), math.sqrt(68), 1e-12, "frobenius norm")
            t.assert_near(A:max_abs(), 8.0, 1e-12, "max abs")
            t.assert_near(A:min_abs_nonzero(), 2.0, 1e-12, "min abs nonzero")

            local diag = A:diagonal()
            local dr, dc = diag:size()
            t.assert_eq(dr, 2, "diagonal rows")
            t.assert_eq(dc, 1, "diagonal cols")
            t.assert_near(diag:get(1, 1), 2.0, 1e-12, "diagonal")
            t.assert_near(diag:get(2, 1), 8.0, 1e-12, "diagonal")

            local triple = B:triple_product_symm(A)
            t.assert_true(triple:is_symmetric(), "triple product symmetry")
            t.assert_near(triple:get(1, 1), 34.0, 1e-12, "triple product")
            t.assert_near(triple:get(1, 2), 16.0, 1e-12, "triple product")
            t.assert_near(triple:get(2, 1), 16.0, 1e-12, "triple product")
            t.assert_near(triple:get(2, 2), 8.0, 1e-12, "triple product")

            local g1 = chemistry.Group.create(2, mol)
            local g2 = chemistry.Group.create(2, mol)
            local gall = chemistry.Group.create_full(2, mol)
            g1:set(1, true)
            g2:set(2, true)
            t.assert_true(g1:get(1), "group set/get")
            t.assert_true(not g1:get(2), "group set/get")

            local block = A:extract_block(g1, g2)
            local br, bc = block:size()
            t.assert_eq(br, 1, "manual block rows")
            t.assert_eq(bc, 1, "manual block cols")
            t.assert_near(block:get(1, 1), 0.0, 1e-12, "manual block value")

            t.assert_near(A:block_trace(gall), 10.0, 1e-12, "block trace")
            t.assert_near(A:block_sum_squares(gall, gall), 68.0, 1e-12, "block sum squares")
            t.assert_near(A:block_norm_fro(gall, gall), math.sqrt(68), 1e-12, "block norm")
            t.assert_near(A:block_mayer_pair(A, gall, gall), 68.0, 1e-12, "block mayer pair")

            block:free()
            g1:free()
            g2:free()
            gall:free()
            triple:free()
            diag:free()
            clamp:free()
            threshold:free()
            C:free()
            cwise:free()
            had:free()
            sqrt_lapack:free()
            sqrt_a:free()
            B:free()
            A:free()
        end)
    end},

    {"ao to mo matrix transform", function()
        with_molecule(function(mol)
            local C = mol:mo_coefficients()
            t.assert_true(C ~= nil, "MO coefficients available")

            local cr, cc = C:size()
            t.assert_eq(cr, mol:get_basis_size(), "MO coefficient rows")

            local S = mol:overlap()
            local S_mo = S:ao_to_mo(C)
            local mr, mc = S_mo:size()
            t.assert_eq(mr, cc, "AO to MO rows")
            t.assert_eq(mc, cc, "AO to MO cols")

            local energies, total = mol:orbital_energies()
            t.assert_true(energies ~= nil, "orbital energies available")
            t.assert_true(total >= #energies, "orbital energies count")
            t.assert_true(#energies > 0, "orbital energies not empty")

            S_mo:free()
            S:free()
            C:free()
        end)
    end},

    {"temporary resource lifecycle", function()
        with_molecule(function(mol)
            local S = mol:overlap()
            local P = mol:density()
            local PS = P * S
            local id = PS._id

            t.assert_eq(mol._registry.matrix[id].refcount, 1, "temporary matrix has one Lua owner")
            PS:free()
            t.assert_eq(mol._registry.matrix[id].refcount, 0, "released temporary matrix has no Lua owners")

            local ok = pcall(function()
                return PS:trace()
            end)
            t.assert_true(not ok, "released matrix rejects further use")

            S:free()
            P:free()
        end)
    end},

    {"symmetric-only matrix power rejects non-symmetric matrices", function()
        with_molecule(function(mol)
            local S = mol:overlap()
            local P = mol:density()
            local PS = P * S

            t.assert_true(not PS:is_symmetric(), "P*S is treated as non-symmetric")
            local ok = pcall(function()
                return PS:pow(0.5)
            end)
            t.assert_true(not ok, "matrix power requires a symmetric matrix")

            PS:free()
            S:free()
            P:free()
        end)
    end},

    {"groups and blocks", function()
        with_molecule(function(mol)
            local carbon = chemistry.Group.from_atom(mol, 1)
            local hydrogen = chemistry.Group.from_atom(mol, 2)
            local both = carbon:union(hydrogen)
            local intersection = carbon:intersection(hydrogen)

            t.assert_true(carbon:count() > 0, "carbon AO group")
            t.assert_true(hydrogen:count() > 0, "hydrogen AO group")
            t.assert_true(both:count() >= carbon:count(), "group union")
            t.assert_eq(intersection:count(), 0, "different atoms do not intersect")

            local S = mol:overlap()
            local block = S:extract_block(carbon, hydrogen)
            local rows, cols = block:size()
            t.assert_eq(rows, carbon:count(), "block rows")
            t.assert_eq(cols, hydrogen:count(), "block cols")

            block:free()
            S:free()
            intersection:free()
            both:free()
            hydrogen:free()
            carbon:free()
        end)
    end},
})

do return end

-- Корень проекта для test-файлов (выставляется runner.lua)
local project_root = os.getenv("PROJECT_ROOT") or "./"

print("=== chemistry_ffi Tests ===\n")

-- ---------------------------------------------------------------------------
-- 1) Загрузка модуля
-- ---------------------------------------------------------------------------
print("[1] Загрузка модуля")
local ok, chemistry = pcall(require, "chemistry_ffi")
assert(ok, "Модуль загружается")
assert(chemistry.Molecule ~= nil, "chemistry.Molecule существует")
assert(chemistry.Matrix ~= nil, "chemistry.Matrix существует")
assert(chemistry.Group ~= nil, "chemistry.Group существует")

-- ---------------------------------------------------------------------------
-- 2) Загрузка молекулы
-- ---------------------------------------------------------------------------
print("\n[2] Загрузка молекулы")
local mol_ok, mol = pcall(chemistry.Molecule.load, project_root .. "examples/fixtures/methane_6-31g.log")
assert(mol_ok, "methane_6-31g.log загружается")

if mol_ok then
    assert(mol ~= nil, "Молекула не nil")
    assert(type(mol.get_num_atoms) == "function", "get_num_atoms — функция")
    assert(type(mol.get_basis_size) == "function", "get_basis_size — функция")
    assert(type(mol.close) == "function", "close — функция")

    local natoms = mol:get_num_atoms()
    assert(natoms == 5, string.format("Число атомов = 5 (получено: %d)", natoms))

    local nbasis = mol:get_basis_size()
    assert(nbasis > 0, string.format("Базис > 0 (получено: %d)", nbasis))

    local basis_name = mol:get_basis_name()
    assert(basis_name ~= nil and basis_name ~= "", "Имя базиса не пустое")

    mol:close()
    assert(mol:is_closed(), "Молекула закрыта")
end

-- ---------------------------------------------------------------------------
-- 3) Матрицы (overlap, density)
-- ---------------------------------------------------------------------------
print("\n[3] Матрицы (overlap, density)")
local mol_ok2, mol2 = pcall(chemistry.Molecule.load, project_root .. "examples/fixtures/methane_6-31g.log")
if mol_ok2 then
    local ov_ok, overlap = pcall(function() return mol2:overlap() end)
    assert(ov_ok, "Overlap загружается")

    if ov_ok then
        local rows, cols = overlap:size()
        assert(rows > 0, string.format("Overlap rows > 0 (%d)", rows))
        assert(cols > 0, string.format("Overlap cols > 0 (%d)", cols))
        assert(rows == cols, "Overlap квадратная")

        local tr = overlap:trace()
        assert(tr ~= nil, "Trace вычисляется")

        local sym = overlap:is_symmetric()
        assert(type(sym) == "boolean", "is_symmetric возвращает boolean")

        local scaled = overlap:scale(2.0)
        assert(scaled ~= nil, "scale(2.0) работает")
        scaled:free()

        local transposed = overlap:transpose()
        assert(transposed ~= nil, "transpose() работает")
        transposed:free()

        overlap:free()
    end

    local den_ok, density = pcall(function() return mol2:density() end)
    assert(den_ok, "Density загружается")
    if den_ok then density:free() end

    mol2:close()
end

-- ---------------------------------------------------------------------------
-- 4) Group
-- ---------------------------------------------------------------------------
print("\n[4] Group")
local mol_ok3, mol3 = pcall(chemistry.Molecule.load, project_root .. "examples/fixtures/methane_6-31g.log")
if mol_ok3 then
    local nbasis = mol3:get_basis_size()
    local g = chemistry.Group.create(nbasis, mol3)
    assert(g ~= nil, "Group.create работает")
    assert(g:nbasis() == nbasis, "Group nbasis совпадает")
    g:free()
    assert(g:is_freed(), "Group освобоён")

    local g_atom_ok, g_atom = pcall(chemistry.Group.from_atom, mol3, 1)
    assert(g_atom_ok, "Group.from_atom работает")
    if g_atom_ok then
        assert(g_atom:count() > 0, "Group count > 0 для атома")
        g_atom:free()
    end

    mol3:close()
end

-- ---------------------------------------------------------------------------
-- 5) Операторы матриц
-- ---------------------------------------------------------------------------
print("\n[5] Операторы матриц (*, +, -)")
local mol_ok4, mol4 = pcall(chemistry.Molecule.load, project_root .. "examples/fixtures/methane_6-31g.log")
if mol_ok4 then
    local overlap = mol4:overlap()
    local density = mol4:density()

    local prod = overlap * density
    assert(prod ~= nil, "overlap * density работает")
    if prod then prod:free() end

    local sum = overlap + density
    assert(sum ~= nil, "overlap + density работает")
    if sum then sum:free() end

    local diff = overlap - density
    assert(diff ~= nil, "overlap - density работает")
    if diff then diff:free() end

    local sc = overlap * 0.5
    assert(sc ~= nil, "overlap * 0.5 работает")
    if sc then sc:free() end

    local neg = -overlap
    assert(neg ~= nil, "-overlap работает")
    if neg then neg:free() end

    overlap:free()
    density:free()
    mol4:close()
end

-- ---------------------------------------------------------------------------
-- 6) Автоочистка (GC)
-- ---------------------------------------------------------------------------
print("\n[6] GC-безопасность")
collectgarbage("collect")
local mol_ok5, mol5 = pcall(chemistry.Molecule.load, project_root .. "examples/fixtures/methane_6-31g.log")
if mol_ok5 then
    do
        local ov = mol5:overlap()
    end
    collectgarbage("collect")
    assert(true, "GC не крашится после orphan матрицы")
    mol5:close()
end

-- ---------------------------------------------------------------------------
-- Итог
-- ---------------------------------------------------------------------------
print(string.format("\n=== Результат: %d/%d пройдено, %d провалено ===", passed, total, failed))
os.exit(failed > 0 and 1 or 0)
