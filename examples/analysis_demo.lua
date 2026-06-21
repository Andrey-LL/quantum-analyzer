-- Демонстрационный шаблон анализа для публичных примеров.
-- Запуск:
-- bin/quantum_analyzer --batch --template examples/analysis_demo.lua --files examples/fixtures/methane_6-31g.log examples/fixtures/methane_sto-3g.log

h1("Quantum Analyzer demo report")

text("This report is generated from trimmed Gaussian fixtures. It demonstrates grouped pivot tables for the same molecule in different basis sets.")

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

Table("Atomic charges and valence", {
    rows = atoms,
    left = {
        col("Atom", atom_index),
        col("Type", atom_symbol),
    },
    across = basis,
    values = {
        val("q Mulliken", function(mol, atom) return mulliken(mol)[atom.id] end, {fmt = "%.4f"}),
        val("q Lowdin", function(mol, atom) return lowdin(mol)[atom.id] end, {fmt = "%.4f"}),
        val("V", function(mol, atom) return valence(mol)[atom.id] end, {fmt = "%.3f"}),
    },
    derived = {
        summary("Sum", {
            ["q Mulliken"] = "sum",
            ["q Lowdin"] = "sum",
            ["V"] = "—",
        }),
        delta("Mean |qM - qL|", "q Mulliken", "q Lowdin", {fmt = "%.4f", mode = "mean_abs"}),
    },
    header = "flat",
})
