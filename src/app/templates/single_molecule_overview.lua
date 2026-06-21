-- Compact one-molecule report with the main atomic, charge and bond indicators.

h1("Краткий анализ молекулы")

local data = cached(function(mol)
    local S = mol:overlap()
    local P = mol:density()
    local PS = qa.compute_PS(P, S)
    local ao2atom = qa.get_ao_atom_mapping(mol)
    local mulliken = qa.mulliken_charges(mol, P, S)
    local lowdin = qa.loewdin_charges(mol, P, S)
    return {
        electrons = PS:trace(),
        mulliken = mulliken,
        lowdin = lowdin,
        valence = qa.compute_all_atomic_valences(mol, PS, ao2atom),
        bonds = qa.compute_all_bond_indices(mol, PS, ao2atom, 0.03),
    }
end)

metric("basis_size", "Размер базиса", basis_size, {grouped = false})
metric("electrons", "Tr(P*S)", data(mol).electrons, {grouped = false})

Table("Атомная сводка", {
    rows = atoms,
    left = {
        col("Атом", atom_label),
        col("Z", atom_z),
    },
    across = basis,
    values = {
        val("$q_{Mulliken}$", function(mol, atom)
            return data(mol).mulliken[atom.id]
        end, {fmt = "%.4f"}),
        val("$q_{Lowdin}$", function(mol, atom)
            return data(mol).lowdin[atom.id]
        end, {fmt = "%.4f"}),
        val("$\\Delta q$", function(mol, atom)
            return data(mol).mulliken[atom.id] - data(mol).lowdin[atom.id]
        end, {fmt = "%.4f"}),
        val("$V$", function(mol, atom)
            return data(mol).valence[atom.id]
        end, {fmt = "%.4f"}),
    },
    derived = {
        summary("Итог", {
            ["$q_{Mulliken}$"] = "sum",
            ["$q_{Lowdin}$"] = "sum",
            ["$\\Delta q$"] = "—",
            ["$V$"] = "—",
        }),
        delta("Средняя abs(qM - qL)", "$q_{Mulliken}$", "$q_{Lowdin}$",
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
