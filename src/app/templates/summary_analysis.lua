-- Сводный анализ по всем базисам одной молекулы.
-- Файл шаблона является кодовым блоком: router читает его как текст и выполняет в sandbox.

h1("Сводный анализ по базисам")

local mulliken = cached(function(mol)
    return qa.mulliken_charges(mol, mol:density(), mol:overlap())
end)

local lowdin = cached(function(mol)
    return qa.loewdin_charges(mol, mol:density(), mol:overlap())
end)

local valence = cached(function(mol)
    local S = mol:overlap()
    local P = mol:density()
    local PS = qa.compute_PS(P, S)
    local ao2atom = qa.get_ao_atom_mapping(mol)
    return qa.compute_all_atomic_valences(mol, PS, ao2atom)
end)

local bonds = cached(function(mol)
    local S = mol:overlap()
    local P = mol:density()
    local PS = qa.compute_PS(P, S)
    local ao2atom = qa.get_ao_atom_mapping(mol)
    local list = qa.compute_all_bond_indices(mol, PS, ao2atom, 0.05)
    local by_key = {}
    for _, bond in ipairs(list) do
        by_key[bond.atom_a .. "-" .. bond.atom_b] = bond
    end
    return {list = list, by_key = by_key}
end)

local lone_pairs = cached(function(mol)
    local S = mol:overlap()
    local P = mol:density()
    local PS = qa.compute_PS(P, S)
    local ao2atom = qa.get_ao_atom_mapping(mol)
    return qa.analyze_lone_pairs(mol, PS, ao2atom)
end)

Table("Заряды и валентность", {
    rows = atoms,
    left = {
        col("Атом", atom_index),
        col("Тип", atom_symbol),
    },
    across = basis,
    values = {
        val("$q_{Mulliken}$", function(mol, atom) return mulliken(mol)[atom.id] end,
            {fmt = "%.4f"}),
        val("$q_{Lowdin}$", function(mol, atom) return lowdin(mol)[atom.id] end,
            {fmt = "%.4f"}),
        val("$V$", function(mol, atom) return valence(mol)[atom.id] end,
            {fmt = "%.3f"}),
    },
    derived = {
        summary("Сумма", {
            ["$q_{Mulliken}$"] = "sum",
            ["$q_{Lowdin}$"] = "sum",
            ["$V$"] = "—",
        }),
        delta("Средняя |qM - qL|", "$q_{Mulliken}$", "$q_{Lowdin}$",
            {fmt = "%.4f", mode = "mean_abs"}),
    },
    header = "flat",
})

Table("Валентность по базисам", {
    rows = atoms,
    left = {
        col("Атом", atom_label),
        col("Z", atom_z),
    },
    across = basis,
    values = {
        val("Валентность", function(mol, atom) return valence(mol)[atom.id] end,
            {fmt = "%.4f"}),
    },
    header = {"metric", "basis"},
})

Table("Индексы связей", {
    rows = function(mol)
        local rows = {}
        for _, bond in ipairs(bonds(mol).list) do
            rows[#rows + 1] = {
                key = bond.atom_a .. "-" .. bond.atom_b,
                label = bond.atom_a .. "-" .. bond.atom_b,
            }
        end
        return rows
    end,
    left = {
        col("Связь", function(_, row) return row.label end),
    },
    across = basis,
    values = {
        val("$I_{AB}$", function(mol, row)
            local bond = bonds(mol).by_key[row.key]
            return bond and bond.index or nil
        end, {fmt = "%.4f"}),
    },
    header = {"metric", "basis"},
})

Table("Заселённости атомных орбиталей", {
    rows = function(mol)
        local rows = {}
        for atom, data in pairs(lone_pairs(mol) or {}) do
            for orbital, _ in ipairs(data.eigenvalues or {}) do
                rows[#rows + 1] = {
                    key = atom .. ":" .. orbital,
                    atom = atom,
                    orbital = orbital,
                    label = "A" .. atom .. " / " .. orbital,
                }
            end
        end
        table.sort(rows, function(a, b) return a.key < b.key end)
        return rows
    end,
    left = {
        col("Орбиталь", function(_, row) return row.label end),
    },
    across = basis,
    values = {
        val("n", function(mol, row)
            local data = lone_pairs(mol)[row.atom]
            return data and data.eigenvalues and data.eigenvalues[row.orbital] or nil
        end, {fmt = "%.4f"}),
    },
    header = {"metric", "basis"},
})
