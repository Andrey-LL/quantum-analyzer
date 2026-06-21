---
title: Quantum Analyzer Pandoc Demo
---

<style>
body{max-width:980px;font-family:Arial,Helvetica,sans-serif;line-height:1.55}
h1,h2,h3{color:#1f2937}table{border-collapse:collapse;margin:1rem 0}
th,td{border-bottom:1px solid #d8dee6;padding:.35rem .65rem}th{background:#f3f6f8}
img{display:block;max-width:100%;margin:1rem 0 1.7rem}.demo-note{color:#4b5563}
</style>

<p class="demo-note">
This document is rendered by Pandoc. The Gaussian code block is executed by
Quantum Analyzer, and the generated Markdown replaces the block in the final
report.
</p>

The example uses one methane calculation and produces two figures: an
atom-grouped overlap heatmap and a Mulliken charge bar chart.

```{.qa}
local file = "examples/fixtures/methane_6-31g.log"
local mol = Molecule.load(file)
local S, P = mol:overlap(), mol:density()
local PS = qa.compute_PS(P, S)
local charges = qa.mulliken_charges(mol, P, S)

h2("Methane calculation summary")
text("Input file: `" .. file .. "`")
text("Basis: " .. mol:get_basis_name() .. ", basis size: " .. mol:get_basis_size())
text(string.format("Tr(P*S): %.4f", PS:trace()))

local rows = {}
for atom = 1, mol:get_num_atoms() do
    rows[#rows + 1] = {
        tostring(atom),
        tostring(mol:get_atomic_number(atom)),
        string.format("%.4f", charges[atom]),
    }
end

h2("Atomic charges")
emit({headers = {"Atom", "Z", "Mulliken charge"}, rows = rows})

h2("Overlap matrix")
heatmap_atoms(S, mol:ao_atom_mapping(), {
    caption = "Overlap matrix grouped by atoms",
    output = "examples/generated/pandoc_overlap_atoms_heatmap.svg",
    link = "pandoc_overlap_atoms_heatmap.svg",
    overwrite = true,
})

local bars = {}
for atom = 1, mol:get_num_atoms() do
    bars[#bars + 1] = {label = atom .. "", value = charges[atom]}
end
h2("Charge distribution")
bar(bars, {
    caption = "Mulliken charges by atom",
    xlabel = "Atom", ylabel = "Charge",
    output = "examples/generated/pandoc_mulliken_charges.svg",
    link = "pandoc_mulliken_charges.svg",
    overwrite = true,
})

PS:free(); S:free(); P:free(); mol:close()
```
