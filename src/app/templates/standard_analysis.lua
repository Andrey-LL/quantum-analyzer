-- Детальный блочный анализ

h1("Детальный блочный анализ")

local data = cached(function(mol)
    local S = mol:overlap()
    local P = mol:density()
    local PS = qa.compute_PS(P, S)
    local ao2atom = qa.get_ao_atom_mapping(mol)
    return {
        S = S,
        P = P,
        PS = PS,
        ao2atom = ao2atom,
        electron_count = PS:trace(),
        valence = qa.compute_all_atomic_valences(mol, PS, ao2atom),
        mulliken = qa.mulliken_charges(mol, P, S),
        lowdin = qa.loewdin_charges(mol, P, S),
        bonds = qa.compute_all_bond_indices(mol, PS, ao2atom, 0.05),
        lone_pairs = qa.analyze_lone_pairs(mol, PS, ao2atom),
    }
end)

metric("basis_size", "Размер базиса", basis_size, {grouped = false})
metric("electrons", "Tr(PS)", data(mol).electron_count, {grouped = false})
matrix_preview("PS", data(mol).PS, {title = "Матрица PS = P*S", max_size = 10})

Table("Атомные показатели", {
    rows = atoms,
    left = {
        col("Атом", atom_label),
        col("Z", atom_z),
    },
    across = basis,
    values = {
        val("$V$", function(mol, atom) return data(mol).valence[atom.id] end,
            {fmt = "%.4f"}),
        val("$q_{Mulliken}$", function(mol, atom) return data(mol).mulliken[atom.id] end,
            {fmt = "%.4f"}),
        val("$q_{Lowdin}$", function(mol, atom) return data(mol).lowdin[atom.id] end,
            {fmt = "%.4f"}),
    },
    derived = {
        summary("Сумма", {
            ["$V$"] = "—",
            ["$q_{Mulliken}$"] = "sum",
            ["$q_{Lowdin}$"] = "sum",
        }),
        delta("Средняя |qM - qL|", "$q_{Mulliken}$", "$q_{Lowdin}$",
            {fmt = "%.4f", mode = "mean_abs"}),
    },
    header = {"metric", "basis"},
})

Table("Индексы связей", {
    rows = function(mol)
        local rows = {}
        for _, bond in ipairs(data(mol).bonds) do
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
            for _, bond in ipairs(data(mol).bonds) do
                if row.key == bond.atom_a .. "-" .. bond.atom_b then
                    return bond.index
                end
            end
            return nil
        end, {fmt = "%.4f"}),
    },
    header = {"metric", "basis"},
})

Table("Заселённости атомных орбиталей", {
    rows = function(mol)
        local rows = {}
        for atom, atom_data in pairs(data(mol).lone_pairs or {}) do
            for orbital, _ in ipairs(atom_data.eigenvalues or {}) do
                rows[#rows + 1] = {
                    key = atom .. ":" .. orbital,
                    atom = atom,
                    orbital = orbital,
                    label = "A" .. atom .. " / " .. orbital,
                }
            end
        end
        return rows
    end,
    left = {
        col("Орбиталь", function(_, row) return row.label end),
    },
    across = basis,
    values = {
        val("n", function(mol, row)
            local atom_data = data(mol).lone_pairs[row.atom]
            return atom_data and atom_data.eigenvalues and atom_data.eigenvalues[row.orbital] or nil
        end, {fmt = "%.4f"}),
    },
    header = {"metric", "basis"},
})
