-- Charge-method comparison report for one molecule across available basis sets.

h1("Сравнение методов зарядов")

local charges = cached(function(mol)
    local S = mol:overlap()
    local P = mol:density()
    return {
        mulliken = qa.mulliken_charges(mol, P, S),
        lowdin = qa.loewdin_charges(mol, P, S),
    }
end)

local function q_delta(mol, atom)
    local q = charges(mol)
    return q.mulliken[atom.id] - q.lowdin[atom.id]
end

Table("Заряды по атомам", {
    rows = atoms,
    left = {
        col("Атом", atom_label),
        col("Z", atom_z),
    },
    across = basis,
    values = {
        val("$q_{Mulliken}$", function(mol, atom)
            return charges(mol).mulliken[atom.id]
        end, {fmt = "%.4f"}),
        val("$q_{Lowdin}$", function(mol, atom)
            return charges(mol).lowdin[atom.id]
        end, {fmt = "%.4f"}),
        val("$q_M - q_L$", q_delta, {fmt = "%.4f"}),
        val("$abs(q_M - q_L)$", function(mol, atom)
            return math.abs(q_delta(mol, atom))
        end, {fmt = "%.4f"}),
    },
    derived = {
        summary("Сводка", {
            ["$q_{Mulliken}$"] = "sum",
            ["$q_{Lowdin}$"] = "sum",
            ["$q_M - q_L$"] = "—",
            ["$abs(q_M - q_L)$"] = "—",
        }),
        delta("Средняя abs(qM - qL)", "$q_{Mulliken}$", "$q_{Lowdin}$",
            {fmt = "%.4f", mode = "mean_abs"}),
    },
    header = {"metric", "basis"},
})

Table("Сводка расхождения методов", {
    rows = function()
        return {
            {
                key = "sum_mulliken",
                label = "Сумма Mulliken",
                fn = function(q)
                    local s = 0
                    for _, v in ipairs(q.mulliken) do s = s + v end
                    return s
                end,
            },
            {
                key = "sum_lowdin",
                label = "Сумма Lowdin",
                fn = function(q)
                    local s = 0
                    for _, v in ipairs(q.lowdin) do s = s + v end
                    return s
                end,
            },
            {
                key = "mean_abs_delta",
                label = "Среднее abs(qM - qL)",
                fn = function(q)
                    local s, n = 0, 0
                    for atom, mulliken in ipairs(q.mulliken) do
                        if q.lowdin[atom] ~= nil then
                            s = s + math.abs(mulliken - q.lowdin[atom])
                            n = n + 1
                        end
                    end
                    return n > 0 and s / n or nil
                end,
            },
        }
    end,
    left = {
        col("Показатель", function(_, row) return row.label end),
    },
    across = basis,
    values = {
        val("Значение", function(mol, row)
            return row.fn(charges(mol))
        end, {fmt = "%.4f"}),
    },
    header = {"metric", "basis"},
})
