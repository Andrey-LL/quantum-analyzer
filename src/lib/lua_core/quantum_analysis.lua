-- ============================================================================
-- КВАНТОВО-ХИМИЧЕСКИЙ АНАЛИЗ
-- ============================================================================
-- Molecule владеет C-ресурсами, а Matrix и Group являются refcount/LRU-дескрипторами.
-- Временные матрицы, возвращаемые анализом, освобождаются через :free().
-- ============================================================================

local chemistry = require("chemistry_ffi")
local qa = {}

-- ============================================================================
-- БАЗОВЫЕ ФУНКЦИИ АНАЛИЗА
-- ============================================================================

--- Вычисляет индекс валентности атома по формуле Дмитриева-Семёнова
-- Формула: V_A = Tr(2·(PS)_AA - (PS)_AA²) = Σ_ā n_ā·(2 - n_ā)
-- где n_ā — заселённости гибридных орбиталей (собственные значения блока (PS)_AA).
-- Источник: Дмитриев, Семёнов, "Квантовая химия", формула валентности.
-- @param PS_AA Matrix: Блок матрицы PS для атома A.
-- @return number: Индекс валентности атома.
function qa.dmitriev_semenov_valence(PS_AA)
    -- V_A = 2·Tr(PS) - Tr(PS²)
    local trace_PS = PS_AA:trace()
    
    local PS_sq = PS_AA * PS_AA
    local trace_PS_sq = PS_sq:trace()
    
    -- Освобождаем временную матрицу PS_sq
    PS_sq:free()
    
    return 2.0 * trace_PS - trace_PS_sq
end

--- Диагонализует блок (PS)_AA и возвращает заселённости гибридных орбиталей
-- Формула: (PS)_AA = U · Λ · U⁺, где Λ — диагональная матрица заселённостей n_ā.
-- Собственные значения n_ā интерпретируются как заселённости естественных гибридных орбиталей.
-- @param PS_AA Matrix: Блок матрицы (PS)_AA.
-- @return table: eigenvalues — массив заселённостей {n_1, n_2, ...}.
-- @return Matrix: eigenvectors — матрица собственных векторов (столбцы = орбитали).
function qa.diagonalize_atomic_block(PS_AA)
    return PS_AA:eigensystem()
end

--- Классифицирует гибридные орбитали по заселённостям
-- Пороги классификации:
--   n ≥ 1.80 → неподелённая пара (lone pair)
--   0.85 ≤ n < 1.80 → валентная орбиталь (valence)
--   n < 0.15 → вакантная орбиталь (vacant)
-- Промежуточные значения игнорируются.
-- @param eigenvalues table: Заселённости орбиталей.
-- @param eigenvectors Matrix: Матрица собственных векторов.
-- @return table: lone_pairs — список {population, orbital_index, type}.
-- @return table: valence_orbitals — список валентных орбиталей.
-- @return table: vacant_orbitals — список вакантных орбиталей.
function qa.classify_hybrid_orbitals(eigenvalues, eigenvectors)
    local lone_pairs = {}
    local valence_orbitals = {}
    local vacant_orbitals = {}

    local TH_LONE = 1.80
    local TH_VALENCE = 0.85
    local TH_VACANT = 0.15

    for i, n in ipairs(eigenvalues) do
        if n >= TH_LONE then
            table.insert(lone_pairs, { population = n, orbital_index = i, type = "lone_pair" })
        elseif n >= TH_VALENCE then
            table.insert(valence_orbitals, { population = n, orbital_index = i, type = "valence" })
        elseif n < TH_VACANT then
            table.insert(vacant_orbitals, { population = n, orbital_index = i, type = "vacant" })
        end
        -- Промежуточные орбитали пропускаются
    end

    return lone_pairs, valence_orbitals, vacant_orbitals
end

--- Вычисляет индекс связи между атомами A и B
-- Формула: I_AB = Tr((PS)_AB · (PS)_BA)
-- Индекс связи характеризует кратность и прочность связи.
-- @param PS_AB Matrix: Блок (PS) от атома A к B.
-- @param PS_BA Matrix: Блок (PS) от атома B к A.
-- @return number: Индекс связи I_AB.
function qa.bond_index(PS_AB, PS_BA)
    local product = PS_AB * PS_BA
    local idx = product:trace()
    
    -- Освобождаем временную матрицу произведения
    product:free()
    
    return idx
end

--- Вычисляет матрицу произведения плотности и перекрытия: PS = P · S
-- Матрица PS используется во многих анализах (Малликен, Майер, валентности).
-- @param density Matrix: Матрица плотности P.
-- @param overlap Matrix: Матрица перекрытия S.
-- @return Matrix: PS = P · S (временная, требует :free()).
function qa.compute_PS(density, overlap)
    return density * overlap
end

--- Получает маппинг базисных функций (AO) на атомы
-- @param molecule Molecule: Объект молекулы.
-- @return table: ao2atom[ao_index] = atom_index.
function qa.get_ao_atom_mapping(molecule)
    return molecule:ao_atom_mapping()
end

--- Внутренняя функция строит списки AO индексов для каждого атома
-- @param molecule Molecule: Объект молекулы.
-- @param ao2atom table: Маппинг AO → атом.
-- @return table: atom_ao_lists[atom_idx] = {ao_indices...}.
local function build_atom_ao_lists(molecule, ao2atom)
    local natoms = molecule:get_num_atoms()
    local nbasis = molecule:get_basis_size()
    local lists = {}
    for i = 1, natoms do lists[i] = {} end
    
    for mu = 1, nbasis do
        local a = ao2atom[mu]
        if a and a >= 1 and a <= natoms then
            table.insert(lists[a], mu)
        end
    end
    return lists
end

-- ============================================================================
-- ПОПУЛЯЦИОННЫЙ АНАЛИЗ И ЗАРЯДЫ
-- Параметры density/overlap опциональны: извлекаются из molecule, если nil.
-- ============================================================================

--- Mulliken популяционный анализ по атомам
-- Формула: G_A = Σ_{μ∈A} (PS)_{μμ}
-- Популяция атома равна сумме диагональных элементов PS для его базисных функций.
-- Источник: Mulliken, J. Chem. Phys. 23, 1833 (1955).
-- @param molecule Molecule: Объект молекулы.
-- @param density Matrix: Матрица плотности P (опционально).
-- @param overlap Matrix: Матрица перекрытия S (опционально).
-- @return table: populations[atom_idx] = популяция.
function qa.mulliken_populations(molecule, density, overlap)
    -- Автоматическое извлечение матриц из молекулы, если не переданы
    density = density or molecule:density()
    overlap = overlap or molecule:overlap()
    
    -- PS = P * S (используем __mul)
    local PS = density * overlap
    
    local ao2atom = molecule:ao_atom_mapping()
    local lists = build_atom_ao_lists(molecule, ao2atom)
    local natoms = molecule:get_num_atoms()
    
    local pops = {}
    for i = 1, natoms do pops[i] = 0.0 end

    -- G_A = Σ_{μ∈A} (PS)_{μμ}
    for A = 1, natoms do
        for _, mu in ipairs(lists[A]) do
            pops[A] = pops[A] + PS:get(mu, mu)
        end
    end
    
    -- Освобождаем временную PS
    PS:free()
    
    return pops
end

--- Вычисление зарядов Малликена
-- Формула: q_A = Z_A - G_A
-- @param molecule Molecule: Объект молекулы.
-- @param density Matrix: Матрица плотности P (опционально).
-- @param overlap Matrix: Матрица перекрытия S (опционально).
-- @return table: charges[atom_idx] = заряд.
function qa.mulliken_charges(molecule, density, overlap)
    local pops = qa.mulliken_populations(molecule, density, overlap)
    local natoms = molecule:get_num_atoms()
    local charges = {}
    
    for A = 1, natoms do
        charges[A] = molecule:get_atomic_number(A) - pops[A]
    end
    
    return charges
end

--- Löwdin популяционный анализ
-- Формула: P_L = S^(1/2) · P · S^(1/2), затем G_A = Σ_{μ∈A} (P_L)_{μμ}
-- Преимущество: меньшая зависимость от базиса по сравнению с Малликеном.
-- Источник: Löwdin, Adv. Quantum Chem. 5, 185 (1970).
-- @param molecule Molecule: Объект молекулы.
-- @param density Matrix: Матрица плотности P (опционально).
-- @param overlap Matrix: Матрица перекрытия S (опционально).
-- @return table: charges[atom_idx] = заряд.
function qa.loewdin_charges(molecule, density, overlap)
    density = density or molecule:density()
    overlap = overlap or molecule:overlap()
    
    local natoms = molecule:get_num_atoms()
    local nbasis = molecule:get_basis_size()
    local ao2atom = molecule:ao_atom_mapping()
    
    -- S^(1/2) через матричную степень
    local sqrt_S = overlap:pow(0.5)
    
    -- P_L = S^(1/2) * P * S^(1/2) (используем __mul)
    local temp = sqrt_S * density
    local P_loew = temp * sqrt_S
    
    -- Освобождаем промежуточные матрицы
    temp:free()
    sqrt_S:free()
    
    local pops = {}
    for i = 1, natoms do pops[i] = 0.0 end
    
    for mu = 1, nbasis do
        local a = ao2atom[mu]
        if a and a >= 1 and a <= natoms then
            pops[a] = pops[a] + P_loew:get(mu, mu)
        end
    end
    
    P_loew:free()
    
    local charges = {}
    for A = 1, natoms do
        charges[A] = molecule:get_atomic_number(A) - pops[A]
    end
    
    return charges
end

-- ============================================================================
-- ПОРЯДКИ СВЯЗЕЙ
-- ============================================================================

--- Матрица порядков связи Виберга
-- Формула: W_AB = Σ_{μ∈A} Σ_{ν∈B} (PS)_{μν}²
-- Источник: Wiberg, Tetrahedron 24, 1083 (1968).
-- @param molecule Molecule: Объект молекулы.
-- @param density Matrix: Матрица плотности P (опционально).
-- @param overlap Matrix: Матрица перекрытия S (опционально).
-- @return table: bond_orders[A][B] = порядок связи.
function qa.wiberg_bond_orders(molecule, density, overlap)
    density = density or molecule:density()
    overlap = overlap or molecule:overlap()
    
    local PS = density * overlap
    local ao2atom = molecule:ao_atom_mapping()
    local lists = build_atom_ao_lists(molecule, ao2atom)
    local natoms = molecule:get_num_atoms()
    
    local bo = {}
    for i = 1, natoms do 
        bo[i] = {}
        for j = 1, natoms do bo[i][j] = 0.0 end 
    end

    for A = 1, natoms do
        for B = 1, natoms do
            for _, mu in ipairs(lists[A]) do
                for _, nu in ipairs(lists[B]) do
                    local v = PS:get(mu, nu)
                    bo[A][B] = bo[A][B] + v * v
                end
            end
        end
    end
    
    PS:free()
    return bo
end

--- Матрица порядков связи Майера
-- Формула: M_AB = Σ_{μ∈A} Σ_{ν∈B} (PS)_{μν} · (PS)_{νμ}
-- Источник: Mayer, Chem. Phys. Lett. 97, 270 (1983).
-- @param molecule Molecule: Объект молекулы.
-- @param density Matrix: Матрица плотности P (опционально).
-- @param overlap Matrix: Матрица перекрытия S (опционально).
-- @return table: bond_orders[A][B] = порядок связи.
function qa.mayer_bond_orders(molecule, density, overlap)
    density = density or molecule:density()
    overlap = overlap or molecule:overlap()
    
    local PS = density * overlap
    local ao2atom = molecule:ao_atom_mapping()
    local lists = build_atom_ao_lists(molecule, ao2atom)
    local natoms = molecule:get_num_atoms()
    
    local bo = {}
    for i = 1, natoms do 
        bo[i] = {}
        for j = 1, natoms do bo[i][j] = 0.0 end 
    end

    for A = 1, natoms do
        for B = 1, natoms do
            for _, mu in ipairs(lists[A]) do
                for _, nu in ipairs(lists[B]) do
                    local v1 = PS:get(mu, nu)
                    local v2 = PS:get(nu, mu)
                    bo[A][B] = bo[A][B] + v1 * v2
                end
            end
        end
    end
    
    PS:free()
    return bo
end

--- Порядок связи Виберга для конкретной пары атомов
-- @param molecule Molecule: Объект молекулы.
-- @param density Matrix: Матрица плотности P (опционально).
-- @param overlap Matrix: Матрица перекрытия S (опционально).
-- @param atom_A number: Индекс атома A.
-- @param atom_B number: Индекс атома B.
-- @return number: Порядок связи или nil.
function qa.wiberg_bond_order_pair(molecule, density, overlap, atom_A, atom_B)
    local natoms = molecule:get_num_atoms()
    -- Проверка границ.
    if atom_A < 0 or atom_A >= natoms or atom_B < 0 or atom_B >= natoms then 
        return nil 
    end
    
    local bo = qa.wiberg_bond_orders(molecule, density, overlap)
    return bo[atom_A + 1][atom_B + 1]
end

--- Индексы валентности атомов
-- Формула: V_A = Σ_{B≠A} M_AB
-- Валентность атома равна сумме порядков всех его связей.
-- @param molecule Molecule: Объект молекулы.
-- @param density Matrix: Матрица плотности P (опционально).
-- @param overlap Matrix: Матрица перекрытия S (опционально).
-- @return table: valences[atom_idx] = индекс валентности.
function qa.atomic_valence_indices(molecule, density, overlap)
    local bo = qa.mayer_bond_orders(molecule, density, overlap)
    local natoms = molecule:get_num_atoms()
    local vals = {}
    
    for A = 1, natoms do
        local s = 0.0
        for B = 1, natoms do 
            if A ~= B then s = s + bo[A][B] end 
        end
        vals[A] = s
    end
    
    return vals
end

-- ============================================================================
-- БЛОЧНЫЙ АНАЛИЗ И ИНДЕКСЫ
-- Параметры PS/ao2atom опциональны: вычисляются/извлекаются из molecule.
-- ============================================================================

--- Вычисляет все индексы связей для молекулы
-- @param molecule Molecule: Объект молекулы.
-- @param PS Matrix: Матрица PS (опционально).
-- @param ao2atom table: Маппинг AO→атом (опционально).
-- @param threshold number: Порог значимости связи (по умолчанию 0.1).
-- @return table: bonds = {{atom_a, atom_b, index}, ...}.
function qa.compute_all_bond_indices(molecule, PS, ao2atom, threshold)
    threshold = threshold or 0.1
    
    -- Если PS не передана, вычисляем из молекулы
    local PS_owned = false
    if not PS then
        local d = molecule:density()
        local s = molecule:overlap()
        PS = d * s
        PS_owned = true
    end
    
    ao2atom = ao2atom or molecule:ao_atom_mapping()
    local natoms = molecule:get_num_atoms()
    local lists = build_atom_ao_lists(molecule, ao2atom)
    local bonds = {}

    for A = 1, natoms do
        for B = A + 1, natoms do
            if #lists[A] > 0 and #lists[B] > 0 then
                local gA = chemistry.Group.from_indices(lists[A], molecule)
                local gB = chemistry.Group.from_indices(lists[B], molecule)

                local PS_AB = PS:extract_block(gA, gB)
                local PS_BA = PS:extract_block(gB, gA)

                local idx = qa.bond_index(PS_AB, PS_BA)
                
                if idx > threshold then
                    table.insert(bonds, { atom_a = A, atom_b = B, index = idx })
                end

                gA:free(); gB:free(); PS_AB:free(); PS_BA:free()
            end
        end
    end
    
    -- Освобождаем PS только если она была создана внутри функции
    if PS_owned then PS:free() end
    
    return bonds
end

--- Вычисляет индексы валентности всех атомов по формуле Дмитриева-Семёнова
-- @param molecule Molecule: Объект молекулы.
-- @param PS Matrix: Матрица PS (опционально).
-- @param ao2atom table: Маппинг AO→атом (опционально).
-- @return table: valences[atom_idx] = индекс валентности.
function qa.compute_all_atomic_valences(molecule, PS, ao2atom)
    local PS_owned = false
    if not PS then
        local d = molecule:density()
        local s = molecule:overlap()
        PS = d * s
        PS_owned = true
    end
    
    ao2atom = ao2atom or molecule:ao_atom_mapping()
    local natoms = molecule:get_num_atoms()
    local lists = build_atom_ao_lists(molecule, ao2atom)
    local vals = {}

    for A = 1, natoms do
        if #lists[A] > 0 then
            local gA = chemistry.Group.from_indices(lists[A], molecule)
            local PS_AA = PS:extract_block(gA, gA)
            
            vals[A] = qa.dmitriev_semenov_valence(PS_AA)
            
            gA:free(); PS_AA:free()
        else
            vals[A] = 0.0
        end
    end
    
    if PS_owned then PS:free() end
    return vals
end

--- Анализ неподелённых электронных пар для всех атомов
-- @param molecule Molecule: Объект молекулы.
-- @param PS Matrix: Матрица PS (опционально).
-- @param ao2atom table: Маппинг AO→атом (опционально).
-- @return table: results[atom_idx] = {lone_pairs, valence, ...}.
function qa.analyze_lone_pairs(molecule, PS, ao2atom)
    local PS_owned = false
    if not PS then
        local d = molecule:density()
        local s = molecule:overlap()
        PS = d * s
        PS_owned = true
    end
    
    ao2atom = ao2atom or molecule:ao_atom_mapping()
    local natoms = molecule:get_num_atoms()
    local lists = build_atom_ao_lists(molecule, ao2atom)
    local res = {}

    for A = 1, natoms do
        if #lists[A] > 0 then
            local gA = chemistry.Group.from_indices(lists[A], molecule)
            local PS_AA = PS:extract_block(gA, gA)
            
            local valence = qa.dmitriev_semenov_valence(PS_AA)
            local evs, evecs = qa.diagonalize_atomic_block(PS_AA)
            local lp, vo, vac = qa.classify_hybrid_orbitals(evs, evecs)

            res[A] = {
                atom_index = A,
                ao_count = #lists[A],
                valence = valence,
                trace_PS = PS_AA:trace(),
                lone_pairs = lp,
                valence_orbitals = vo,
                vacant_orbitals = vac,
                eigenvalues = evs,
                eigenvectors = evecs
            }
            gA:free(); PS_AA:free()
        else
            res[A] = { atom_index = A, ao_count = 0, valence = 0.0, lone_pairs = {}, valence_orbitals = {}, vacant_orbitals = {} }
        end
    end
    
    if PS_owned then PS:free() end
    return res
end

-- ============================================================================
-- НИЗКОУРОВНЕВЫЕ МАТРИЧНЫЕ ПРЕОБРАЗОВАНИЯ
-- ============================================================================

--- Вычисляет эрмитову матрицу плотности Лёвдина
-- Формула: P_L = S^(-1/2) · P · S^(-1/2)
-- @param density Matrix: Матрица плотности P.
-- @param overlap Matrix: Матрица перекрытия S.
-- @return Matrix: P_L (требует :free()).
function qa.compute_loewdin_density(density, overlap)
    -- X = S^(-1/2)
    local X = overlap:pow(-0.5)
    
    -- P_L = X * P * X (используем __mul)
    local XT_P = X * density
    local P_loew = XT_P * X
    
    X:free(); XT_P:free()
    return P_loew
end

--- Вычисляет неэрмитову матрицу плотности PS
-- Формула: PS = P · S
-- @param density Matrix: Матрица плотности P.
-- @param overlap Matrix: Матрица перекрытия S.
-- @return Matrix: PS (требует :free()).
function qa.compute_nonhermitian_density(density, overlap)
    return density * overlap
end

--- SVD анализ связи между атомами A и B
-- Формула: SVD блоков перекрытия S_AB и S_BA.
-- Анализирует сингулярные числа для оценки кратности связи.
-- @param overlap Matrix: Матрица перекрытия S.
-- @param ao2atom table: Маппинг AO→атом.
-- @param atom_A number: Индекс атома A.
-- @param atom_B number: Индекс атома B.
-- @return Matrix: U_A, Matrix: U_B, table: singular_values.
function qa.svd_bond_analysis(overlap, ao2atom, atom_A, atom_B)
    local molecule = overlap.molecule
    local nbasis = molecule:get_basis_size()
    
    -- Сбор AO индексов для атомов.
    local tA, tB = atom_A + 1, atom_B + 1
    local ao_A = {}
    local ao_B = {}
    
    for mu = 1, nbasis do
        local a = ao2atom[mu]
        if a == tA then table.insert(ao_A, mu)
        elseif a == tB then table.insert(ao_B, mu) end
    end

    if #ao_A == 0 or #ao_B == 0 then return nil, nil, nil end

    local gA = chemistry.Group.from_indices(ao_A, molecule)
    local gB = chemistry.Group.from_indices(ao_B, molecule)
    
    local S_AB = overlap:extract_block(gA, gB)
    local S_BA = overlap:extract_block(gB, gA)
    
    -- M1 = S_AB * S_BA, M2 = S_BA * S_AB
    local M1 = S_AB * S_BA
    local M2 = S_BA * S_AB

    local U_A, S_A, Vt_A = M1:svd()
    local U_B, S_B, Vt_B = M2:svd()

    local svals = {}
    local r, c = S_A:size()
    for i = 1, math.min(r, c) do
        local v = S_A:get(i, i)
        if v > 1e-6 then table.insert(svals, v) end
    end

    -- Освобождаем все временные ресурсы.
    gA:free(); gB:free(); S_AB:free(); S_BA:free()
    M1:free(); M2:free()
    S_A:free(); Vt_A:free(); S_B:free(); Vt_B:free()
    
    return U_A, U_B, svals
end

-- ============================================================================
-- ЭКСПОРТ: АНАЛИЗ + FFI-ОБЁРТКА
-- ============================================================================
local exports = {}

for k, v in pairs(qa) do exports[k] = v end

-- Реэкспортируем chemistry_ffi: Molecule, Matrix, Group и низкоуровневый C API.
for k, v in pairs(chemistry) do exports[k] = v end

return exports
