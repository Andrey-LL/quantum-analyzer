-- Анализ неподелённых электронных пар по формулам Дмитриева-Семёнова.

h1("Анализ неподелённых электронных пар")

local ps_data = cached(function(mol)
    local S = mol:overlap()
    local P = mol:density()
    local PS = qa.compute_PS(P, S)
    local ao2atom = qa.get_ao_atom_mapping(mol)
    return {
        PS = PS,
        ao2atom = ao2atom,
        valence = qa.compute_all_atomic_valences(mol, PS, ao2atom),
        lone_pairs = qa.analyze_lone_pairs(mol, PS, ao2atom),
        bonds = qa.compute_all_bond_indices(mol, PS, ao2atom, 0.1),
    }
end)

Table("Валентность атомов", {
    rows = atoms,
    left = {
        col("Атом", atom_label),
        col("Z", atom_z),
    },
    across = basis,
    values = {
        val("$V$", function(mol, atom)
            return ps_data(mol).valence[atom.id]
        end, {fmt = "%.4f"}),
    },
    header = {"metric", "basis"},
})

Table("Неподелённые пары", {
    rows = function(mol)
        local rows = {}
        for atom, data in pairs(ps_data(mol).lone_pairs or {}) do
            for pair, _ in ipairs(data.lone_pairs or {}) do
                rows[#rows + 1] = {
                    key = atom .. ":" .. pair,
                    atom = atom,
                    pair = pair,
                    label = "A" .. atom .. " / LP" .. pair,
                }
            end
        end
        return rows
    end,
    left = {
        col("Пара", function(_, row) return row.label end),
    },
    across = basis,
    values = {
        val("n", function(mol, row)
            local atom_data = ps_data(mol).lone_pairs[row.atom]
            local lp = atom_data and atom_data.lone_pairs and atom_data.lone_pairs[row.pair]
            return lp and lp.population or nil
        end, {fmt = "%.4f"}),
    },
    header = {"metric", "basis"},
})

Table("Индексы связей", {
    rows = function(mol)
        local rows = {}
        for _, bond in ipairs(ps_data(mol).bonds) do
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
            for _, bond in ipairs(ps_data(mol).bonds) do
                if row.key == bond.atom_a .. "-" .. bond.atom_b then
                    return bond.index
                end
            end
            return nil
        end, {fmt = "%.4f"}),
    },
    header = {"metric", "basis"},
})
