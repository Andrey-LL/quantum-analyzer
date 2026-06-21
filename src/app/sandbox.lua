-- Среда выполнения Lua-блоков анализа Gaussian.

local function dirname(path)
  return (path or "."):gsub("\\", "/"):match("^(.*)/[^/]*$") or "."
end

local function add_unique(list, seen, value)
  if value and value ~= "" and not seen[value] then
    seen[value] = true
    list[#list + 1] = value
  end
end

local base_dir = "."
if arg and arg[0] then
  base_dir = dirname(arg[0])
end

local module_dir = dirname((debug.getinfo(1, "S").source or ""):gsub("^@", ""))
local plugin_dirs, seen_plugin_dirs = {}, {}
add_unique(plugin_dirs, seen_plugin_dirs, base_dir .. "/plugins")
add_unique(plugin_dirs, seen_plugin_dirs, base_dir .. "/../plugins")
add_unique(plugin_dirs, seen_plugin_dirs, module_dir .. "/plugins")
add_unique(plugin_dirs, seen_plugin_dirs, module_dir .. "/../plugins")

for _, plugin_dir in ipairs(plugin_dirs) do
  package.path = plugin_dir .. "/?.lua;" .. plugin_dir .. "/*/?.lua;" .. package.path
end

local chemistry = require("chemistry_ffi")
local qa = require("quantum_analysis")

local plugins = {}
do
  local loaded = {}
  local is_windows = package.config:sub(1, 1) == "\\"
  local function shquote(path)
    return '"' .. tostring(path):gsub('"', '\\"') .. '"'
  end
  local function list_command(dir)
    if is_windows then
      return "dir /b " .. shquote((dir:gsub("/", "\\"))) .. " 2>nul"
    end
    return "ls -1 " .. shquote(dir) .. " 2>/dev/null"
  end

  for _, plugin_dir in ipairs(plugin_dirs) do
    local p = io.popen(list_command(plugin_dir))
    if p then
      for dir in p:lines() do
        if not loaded[dir] then
          package.path = plugin_dir .. "/" .. dir .. "/?.lua;" .. package.path
          local mod_path = plugin_dir .. "/" .. dir .. "/" .. dir .. ".lua"
          local f = io.open(mod_path, "r")
          if f then
            f:close()
            local ok, mod = pcall(require, dir)
            if ok and type(mod) == "table" then
              loaded[dir] = true
              table.insert(plugins, mod)
            end
          end
        end
      end
      p:close()
    end
  end
end

local function install_plugin_exports(env, ctx)
  for _, plugin in ipairs(plugins) do
    if type(plugin.install) == "function" then
      plugin.install(env, ctx)
    elseif type(plugin.exports) == "table" then
      local exports = plugin.exports
      for name, value in pairs(exports) do
        if type(name) == "string" and env[name] == nil then
          env[name] = value
        end
      end
    else
      for name, value in pairs(plugin) do
        if type(name) == "string"
            and name:sub(1, 1) ~= "_"
            and name ~= "install"
            and type(value) == "function"
            and env[name] == nil then
          env[name] = value
        end
      end
    end
  end
end

_LOG_FILES = {}

local function fmt_value(value)
  if type(value) == "number" then
    return string.format("%.4f", value)
  end
  if value == nil then return "" end
  return tostring(value)
end

local function basename(path)
  return (path or ""):gsub("\\", "/"):match("([^/]+)$") or tostring(path)
end

local atom_symbols = {
  [1] = "H", [2] = "He", [3] = "Li", [4] = "Be", [5] = "B",
  [6] = "C", [7] = "N", [8] = "O", [9] = "F", [10] = "Ne",
  [11] = "Na", [12] = "Mg", [13] = "Al", [14] = "Si", [15] = "P",
  [16] = "S", [17] = "Cl", [18] = "Ar", [30] = "Zn",
}

local function chemical_symbol(z)
  return atom_symbols[z] or ("Z" .. tostring(z or "?"))
end

local function table_to_md(t)
  if not t or not t.headers then return "" end
  local res = {}
  local function line(cells)
    local out = {}
    for i = 1, #cells do out[i] = tostring(cells[i] or "") end
    return "| " .. table.concat(out, " | ") .. " |"
  end

  table.insert(res, line(t.headers))
  table.insert(res, "|" .. string.rep("---|", #t.headers))
  for _, r in ipairs(t.rows or {}) do
    table.insert(res, line(r))
  end
  return table.concat(res, "\n")
end

local function atom_rows(mol)
  local rows = {}
  for atom = 1, mol:get_num_atoms() do
    local z = mol:get_atomic_number(atom)
    rows[#rows + 1] = {
      key = tostring(atom),
      id = atom,
      atom = atom,
      index = atom,
      Z = z,
      z = z,
      symbol = chemical_symbol(z),
      label = chemical_symbol(z) .. " (" .. tostring(atom) .. ")",
    }
  end
  return rows
end

local function normalize_columns(items)
  local cols = {}
  for _, item in ipairs(items or {}) do
    cols[#cols + 1] = {
      label = item.label or item.name or item[1],
      fn = item.fn or item.value or item[2],
      opts = item.opts or item[3] or {},
    }
  end
  return cols
end

local function normalize_derived(items)
  local rows = {}
  for _, item in ipairs(items or {}) do
    rows[#rows + 1] = item
  end
  return rows
end

local function normalize_header(header)
  if type(header) == "table" then
    return {mode = "levels", levels = header}
  end
  if header == "spanned" then
    return {mode = "levels", levels = {"metric", "basis"}}
  end
  if header == "basis_metric" then
    return {mode = "levels", levels = {"basis", "metric"}}
  end
  return {mode = "flat", levels = {"metric", "basis"}}
end

local function call_table_fn(fn, record, row)
  if type(fn) == "function" then
    return fn(record.mol, row, record)
  end
  if type(fn) == "string" and row then
    return row[fn]
  end
  return fn
end

local function table_pivot_value(spec, record)
  local across = spec.across or spec.pivot or "basis"
  local value
  if type(across) == "function" then
    value = across(record.mol, record)
  elseif across == "basis" then
    value = record.basis_name
  elseif across == "basis_size" then
    value = record.basis_size
  else
    value = across
  end

  local label = tostring(value or record.basis_name or basename(record.file))
  return {
    key = label,
    label = label,
    order = record.basis_size or label,
  }
end

local function resolve_table_rows(spec, record)
  local rows = spec.rows or "atoms"
  if rows == "atoms" then
    return atom_rows(record.mol)
  end
  if type(rows) == "function" then
    return rows(record.mol, record)
  end
  return rows or {}
end

local function build_pivot_table_block(title, spec, record)
  spec = spec or {}
  local left = normalize_columns(spec.left or spec.stub or {
    {"Атом", function(_, row) return row.label or row.id end},
  })
  local values = normalize_columns(spec.values or spec.metrics)
  local pivot = table_pivot_value(spec, record)
  local block = {
    type = "pivot_table",
    id = spec.id or title,
    title = title,
    header = normalize_header(spec.header or spec.header_layout or spec.header_mode or "flat"),
    transpose = spec.transpose == true,
    fill = spec.fill or "—",
    summary_label = spec.summary_label or spec.footer,
    left = left,
    values = values,
    derived = normalize_derived(spec.derived or spec.after or {}),
    pivot = pivot,
    rows = {},
    row_order = {},
    cells = {},
  }

  for index, row in ipairs(resolve_table_rows(spec, record)) do
    local row_key = tostring(row.key or row.id or row.atom or index)
    block.row_order[#block.row_order + 1] = row_key
    block.rows[row_key] = {key = row_key, left = {}, source = row}

    for i, col in ipairs(left) do
      block.rows[row_key].left[i] = call_table_fn(col.fn, record, row)
    end

    block.cells[row_key] = {}
    for value_index, value_spec in ipairs(values) do
      block.cells[row_key][value_index] = call_table_fn(value_spec.fn, record, row)
    end
  end

  return block
end

local function matrix_to_md(m, max_size)
  max_size = max_size or 10
  if not m then return "" end

  local rows, cols = m:size()
  local res = {}
  table.insert(res, string.format("Размер: %d×%d  |  Tr = %.4f  |  Симметрична: %s",
    rows, cols, m:trace() or 0, m:is_symmetric() and "да" or "нет"))
  table.insert(res, "")

  local show_rows = math.min(rows, max_size)
  local show_cols = math.min(cols, max_size)
  local header = {" "}
  for j = 1, show_cols do header[#header + 1] = j end
  local body = {headers = header, rows = {}}
  for i = 1, show_rows do
    local row = {i}
    for j = 1, show_cols do
      row[#row + 1] = string.format("%.3f", m:get(i, j))
    end
    body.rows[#body.rows + 1] = row
  end
  table.insert(res, table_to_md(body))
  if rows > show_rows or cols > show_cols then
    table.insert(res, "")
    table.insert(res, string.format("*Показан preview %d×%d из %d×%d.*", show_rows, show_cols, rows, cols))
  end
  return table.concat(res, "\n")
end

local function make_molecule_signature(mol)
  local counts = {}
  local ordered = {}
  local atomic_numbers = {}
  for atom = 1, mol:get_num_atoms() do
    local z = mol:get_atomic_number(atom)
    atomic_numbers[atom] = z
    ordered[atom] = tostring(z)
    counts[z] = (counts[z] or 0) + 1
  end

  local zs = {}
  for z in pairs(counts) do zs[#zs + 1] = z end
  table.sort(zs)

  local parts = {}
  for _, z in ipairs(zs) do
    parts[#parts + 1] = tostring(z) .. ":" .. tostring(counts[z])
  end

  return table.concat(parts, ";"), table.concat(ordered, ","), atomic_numbers
end

local DocumentBuilder = {}
DocumentBuilder.__index = DocumentBuilder

function DocumentBuilder.new(record)
  return setmetatable({
    record = record,
    blocks = {},
    block_index = {},
    wants_grouped = false,
    disable_grouping = false,
  }, DocumentBuilder)
end

function DocumentBuilder:_add(block)
  self.blocks[#self.blocks + 1] = block
  if block.id then self.block_index[block.id] = block end
  return block
end

function DocumentBuilder:h(level, text)
  return self:_add({type = "heading", level = level, text = tostring(text or "")})
end

function DocumentBuilder:h1(text) return self:h(1, text) end
function DocumentBuilder:h2(text) return self:h(2, text) end

function DocumentBuilder:text(text, opts)
  opts = opts or {}
  return self:_add({type = "paragraph", text = tostring(text or ""), merge = opts.merge or "first"})
end

function DocumentBuilder:metric(id, label, value, opts)
  opts = opts or {}
  if opts.grouped ~= false then self.wants_grouped = true end
  return self:_add({
    type = "metric",
    id = id,
    label = label or id,
    value = value,
    grouped = opts.grouped ~= false,
  })
end

function DocumentBuilder:atom_table(id, opts)
  opts = opts or {}
  self.wants_grouped = true
  local block = self.block_index[id]
  if not block then
    block = self:_add({
      type = "atom_table",
      id = id,
      title = opts.title or id,
      value_name = opts.value_name or "Значение",
      values = {},
      opts = opts,
    })
  end
  return block
end

function DocumentBuilder:atom_value(table_id, atom_index, value, opts)
  local block = self.block_index[table_id] or self:atom_table(table_id, opts)
  block.values[atom_index] = value
end

function DocumentBuilder:bond_table(id, opts)
  opts = opts or {}
  self.wants_grouped = true
  local block = self.block_index[id]
  if not block then
    block = self:_add({
      type = "bond_table",
      id = id,
      title = opts.title or id,
      value_name = opts.value_name or "Значение",
      values = {},
      opts = opts,
    })
  end
  return block
end

function DocumentBuilder:bond_value(table_id, atom_a, atom_b, value, opts)
  local block = self.block_index[table_id] or self:bond_table(table_id, opts)
  local a, b = atom_a, atom_b
  if a > b then a, b = b, a end
  block.values[tostring(a) .. "-" .. tostring(b)] = {atom_a = a, atom_b = b, value = value}
end

function DocumentBuilder:table(id, opts)
  opts = opts or {}
  if opts.grouped then self.wants_grouped = true end
  local block = self.block_index[id]
  if not block then
    block = self:_add({
      type = "table",
      id = id,
      title = opts.title or id,
      rows = {},
      row_order = {},
      col_order = {},
      grouped = opts.grouped == true,
      opts = opts,
    })
  end
  return block
end

function DocumentBuilder:set_cell(table_id, row_key, col_key, value, opts)
  local block = self.block_index[table_id] or self:table(table_id, opts)
  if not block.rows[row_key] then
    block.rows[row_key] = {}
    block.row_order[#block.row_order + 1] = row_key
  end
  if block.rows[row_key][col_key] == nil then
    local seen = false
    for _, key in ipairs(block.col_order) do
      if key == col_key then seen = true; break end
    end
    if not seen then block.col_order[#block.col_order + 1] = col_key end
  end
  block.rows[row_key][col_key] = value
end

function DocumentBuilder:figure(id, path, caption, opts)
  opts = opts or {}
  return self:_add({type = "figure", id = id, path = path, caption = caption, merge = opts.merge or "first"})
end

function DocumentBuilder:matrix_preview(id, matrix, opts)
  opts = opts or {}
  local rows, cols = matrix:size()
  return self:_add({
    type = "matrix_preview",
    id = id,
    title = opts.title or id,
    markdown = matrix_to_md(matrix, opts.max_size or 6),
    rows = rows,
    cols = cols,
  })
end

function DocumentBuilder:document()
  return {
    blocks = self.blocks,
    wants_grouped = self.wants_grouped,
    disable_grouping = self.disable_grouping,
  }
end

local function format_table_value(value, opts, fill)
  if value == nil then return fill or "" end
  if type(value) == "number" then
    return string.format((opts and opts.fmt) or "%.4f", value)
  end
  return tostring(value)
end

local function markdown_cell(value)
  return tostring(value or ""):gsub("\n", " "):gsub("|", "\\|")
end

local function collect_pivot_table(blocks)
  local first = blocks[1]
  local state = {
    title = first.title,
    header = first.header,
    transpose = first.transpose,
    fill = first.fill,
    summary_label = first.summary_label,
    left = first.left,
    values = first.values,
    derived = first.derived or {},
    rows = {},
    row_order = {},
    pivots = {},
    pivot_order = {},
    cells = {},
  }

  local function add_pivot(pivot)
    if state.pivots[pivot.key] then return end
    state.pivots[pivot.key] = pivot
    state.pivot_order[#state.pivot_order + 1] = pivot.key
  end

  local function add_row(row_key, row)
    if state.rows[row_key] then return end
    state.rows[row_key] = row
    state.row_order[#state.row_order + 1] = row_key
  end

  for _, block in ipairs(blocks) do
    add_pivot(block.pivot)
    for _, row_key in ipairs(block.row_order or {}) do
      add_row(row_key, block.rows[row_key])
      state.cells[row_key] = state.cells[row_key] or {}
      for value_index, value in pairs(block.cells[row_key] or {}) do
        state.cells[row_key][value_index] = state.cells[row_key][value_index] or {}
        state.cells[row_key][value_index][block.pivot.key] = value
      end
    end
  end

  table.sort(state.pivot_order, function(a, b)
    local pa, pb = state.pivots[a], state.pivots[b]
    if pa.order == pb.order then
      return tostring(pa.label) < tostring(pb.label)
    end
    if type(pa.order) == "number" and type(pb.order) == "number" then
      return pa.order < pb.order
    end
    return tostring(pa.order) < tostring(pb.order)
  end)

  table.sort(state.row_order, function(a, b)
    local function parts(key)
      local out = {}
      for p in tostring(key):gmatch("[^:%-]+") do
        out[#out + 1] = tonumber(p) or p
      end
      return out
    end
    local pa, pb = parts(a), parts(b)
    for i = 1, math.max(#pa, #pb) do
      if pa[i] == nil then return true end
      if pb[i] == nil then return false end
      if pa[i] ~= pb[i] then return pa[i] < pb[i] end
    end
    return tostring(a) < tostring(b)
  end)

  return state
end

local function pivot_metric_index(state)
  local index = {}
  for i, value_spec in ipairs(state.values) do
    index[i] = i
    if value_spec.label then index[value_spec.label] = i end
    if value_spec.name then index[value_spec.name] = i end
  end
  return index
end

local function pivot_columns(state)
  local levels = state.header.levels or {"metric", "basis"}
  local cols = {}
  local function add(value_index, pivot_key)
    local value_spec = state.values[value_index]
    local pivot = state.pivots[pivot_key]
    cols[#cols + 1] = {
      value_index = value_index,
      pivot_key = pivot_key,
      opts = value_spec.opts or {},
      labels = {
        metric = value_spec.label or tostring(value_index),
        basis = pivot.label or pivot_key,
      },
    }
  end

  if levels[1] == "basis" then
    for _, pivot_key in ipairs(state.pivot_order) do
      for value_index = 1, #state.values do add(value_index, pivot_key) end
    end
  else
    for value_index = 1, #state.values do
      for _, pivot_key in ipairs(state.pivot_order) do add(value_index, pivot_key) end
    end
  end
  return cols
end

local function pivot_value(state, row_key, col)
  return state.cells[row_key]
      and state.cells[row_key][col.value_index]
      and state.cells[row_key][col.value_index][col.pivot_key]
end

local function sum_column(state, col)
  local sum, seen = 0, false
  for _, row_key in ipairs(state.row_order) do
    local value = pivot_value(state, row_key, col)
    if type(value) == "number" then
      sum = sum + value
      seen = true
    end
  end
  return seen and sum or nil
end

local function delta_column(state, row_spec, col)
  local metric_index = pivot_metric_index(state)
  local left_index = metric_index[row_spec.left or row_spec[2]]
  local right_index = metric_index[row_spec.right or row_spec[3]]
  if not left_index or not right_index then return nil end

  local sum, count = 0, 0
  for _, row_key in ipairs(state.row_order) do
    local left = state.cells[row_key] and state.cells[row_key][left_index] and state.cells[row_key][left_index][col.pivot_key]
    local right = state.cells[row_key] and state.cells[row_key][right_index] and state.cells[row_key][right_index][col.pivot_key]
    if type(left) == "number" and type(right) == "number" then
      local d = left - right
      if (row_spec.mode or row_spec.opts and row_spec.opts.mode or "mean_abs") == "mean_abs" then
        d = math.abs(d)
      end
      sum = sum + d
      count = count + 1
    end
  end

  if count == 0 then return nil end
  local mode = row_spec.mode or (row_spec.opts and row_spec.opts.mode) or "mean_abs"
  if mode == "sum" or mode == "sum_abs" then return sum end
  return sum / count
end

local function derived_row_value(state, row_spec, col)
  if row_spec.kind == "summary" then
    local metric = col.labels.metric
    local rules = row_spec.rules or row_spec.values or {}
    local rule = rules[metric] or rules[col.value_index] or row_spec.rule
    if rule == "sum" then return sum_column(state, col), col.opts end
    if type(rule) == "function" then return rule(state, col), row_spec.opts end
    if rule ~= nil then return rule, row_spec.opts end
    return nil, row_spec.opts
  end

  if row_spec.kind == "delta" then
    if col.labels.metric ~= (row_spec.target or row_spec.left or row_spec[2]) then
      return nil, row_spec.opts
    end
    return delta_column(state, row_spec, col), row_spec.opts
  end

  if type(row_spec.fn) == "function" then
    return row_spec.fn(state, col), row_spec.opts
  end

  local rule = col.opts.summary or col.opts.bottom
  if row_spec.legacy_summary and rule == "sum" then return sum_column(state, col), col.opts end
  if row_spec.legacy_summary and rule then return rule, col.opts end
  return nil, row_spec.opts
end

local function pivot_derived_rows(state, cols)
  local rows = {}
  local specs = {}
  if state.summary_label then
    specs[#specs + 1] = {kind = "summary", legacy_summary = true, label = state.summary_label}
  end
  for _, row_spec in ipairs(state.derived or {}) do specs[#specs + 1] = row_spec end

  for _, row_spec in ipairs(specs) do
    local row = {}
    row[#row + 1] = row_spec.label or row_spec[1] or "Итог"
    for _ = 2, #state.left do row[#row + 1] = "—" end
    for _, col in ipairs(cols) do
      local value, opts = derived_row_value(state, row_spec, col)
      row[#row + 1] = format_table_value(value, opts or row_spec.opts or col.opts, state.fill)
    end
    rows[#rows + 1] = row
  end

  return rows
end

local function render_grid(title, header_rows, prefix_rows, body_rows, left_count)
  local lines = {}
  lines[#lines + 1] = "## " .. (title or "Таблица")
  lines[#lines + 1] = ""

  local column_count = #header_rows[1]
  local sep = {}
  for i = 1, column_count do
    sep[i] = i <= left_count and ":--" or "--:"
  end

  local function line(cells)
    local out = {}
    for i = 1, column_count do out[i] = markdown_cell(cells[i]) end
    return "| " .. table.concat(out, " | ") .. " |"
  end

  for _, row in ipairs(header_rows) do lines[#lines + 1] = line(row) end
  lines[#lines + 1] = "| " .. table.concat(sep, " | ") .. " |"
  for _, row in ipairs(prefix_rows or {}) do lines[#lines + 1] = line(row) end
  for _, row in ipairs(body_rows or {}) do lines[#lines + 1] = line(row) end
  return table.concat(lines, "\n")
end

local function transpose_grid(header_rows, prefix_rows, body_rows)
  local source = {}
  for _, row in ipairs(header_rows or {}) do source[#source + 1] = row end
  for _, row in ipairs(prefix_rows or {}) do source[#source + 1] = row end
  for _, row in ipairs(body_rows or {}) do source[#source + 1] = row end

  local width = 0
  for _, row in ipairs(source) do width = math.max(width, #row) end

  local transposed = {}
  for col = 1, width do
    local row = {}
    for line_index = 1, #source do row[#row + 1] = source[line_index][col] or "" end
    transposed[#transposed + 1] = row
  end

  return {transposed[1] or {" "}}, {}, {unpack(transposed, 2)}
end

local function render_pivot_table(blocks)
  local state = collect_pivot_table(blocks)
  local cols = pivot_columns(state)
  local header_rows = {}
  local prefix_rows = {}

  if state.header.mode == "levels" then
    local levels = state.header.levels
    for level_index, level in ipairs(levels) do
      local row = {}
      for _, col in ipairs(state.left) do
        row[#row + 1] = level_index == 1 and (col.label or "") or ""
      end
      local previous = nil
      for _, col in ipairs(cols) do
        local label = col.labels[level] or ""
        row[#row + 1] = label == previous and "" or label
        previous = label
      end
      if level_index == 1 then
        header_rows[#header_rows + 1] = row
      else
        prefix_rows[#prefix_rows + 1] = row
      end
    end
  else
    local header = {}
    for _, col in ipairs(state.left) do header[#header + 1] = col.label or "" end
    for _, col in ipairs(cols) do
      header[#header + 1] = string.format("%s (%s)", col.labels.metric, col.labels.basis)
    end
    header_rows[#header_rows + 1] = header
  end

  local body_rows = {}
  for _, row_key in ipairs(state.row_order) do
    local row_data = state.rows[row_key]
    local row = {}
    for i = 1, #state.left do row[#row + 1] = row_data.left[i] end
    for _, col in ipairs(cols) do
      row[#row + 1] = format_table_value(pivot_value(state, row_key, col), col.opts, state.fill)
    end
    body_rows[#body_rows + 1] = row
  end

  for _, row in ipairs(pivot_derived_rows(state, cols)) do body_rows[#body_rows + 1] = row end

  local left_count = #state.left
  if state.transpose then
    header_rows, prefix_rows, body_rows = transpose_grid(header_rows, prefix_rows, body_rows)
    left_count = 1
  end

  return render_grid(state.title, header_rows, prefix_rows, body_rows, left_count)
end

local function render_block_single(block, record)
  if block.type == "heading" then
    return string.rep("#", block.level or 1) .. " " .. block.text
  elseif block.type == "paragraph" then
    return block.text
  elseif block.type == "metric" then
    return string.format("**%s:** %s", block.label or block.id, fmt_value(block.value))
  elseif block.type == "atom_table" then
    local rows = {}
    for atom = 1, #(record.atomic_numbers or {}) do
      rows[#rows + 1] = {atom, record.atomic_numbers[atom] or "", fmt_value(block.values[atom])}
    end
    return "## " .. (block.title or block.id) .. "\n\n" .. table_to_md({
      headers = {"Атом", "Z", block.value_name or "Значение"},
      rows = rows,
    })
  elseif block.type == "bond_table" then
    local rows = {}
    local keys = {}
    for key in pairs(block.values or {}) do keys[#keys + 1] = key end
    table.sort(keys)
    for _, key in ipairs(keys) do
      local value = block.values[key]
      rows[#rows + 1] = {value.atom_a, value.atom_b, fmt_value(value.value)}
    end
    return "## " .. (block.title or block.id) .. "\n\n" .. table_to_md({
      headers = {"Атом A", "Атом B", block.value_name or "Значение"},
      rows = rows,
    })
  elseif block.type == "table" then
    local headers = {block.opts.row_header or "Строка"}
    for _, col in ipairs(block.col_order) do headers[#headers + 1] = col end
    local rows = {}
    for _, row_key in ipairs(block.row_order) do
      local row = {row_key}
      for _, col in ipairs(block.col_order) do row[#row + 1] = fmt_value(block.rows[row_key][col]) end
      rows[#rows + 1] = row
    end
    return "## " .. (block.title or block.id) .. "\n\n" .. table_to_md({headers = headers, rows = rows})
  elseif block.type == "figure" then
    return string.format("![%s](%s)", block.caption or "", block.path or "")
  elseif block.type == "matrix_preview" then
    return "## " .. (block.title or block.id) .. "\n\n" .. block.markdown
  elseif block.type == "pivot_table" then
    return render_pivot_table({block})
  end
  return nil
end

local function render_record(record)
  local out = {}
  out[#out + 1] = "# " .. basename(record.file)
  out[#out + 1] = ""
  out[#out + 1] = "**Файл:** `" .. tostring(record.file) .. "`"
  out[#out + 1] = "**Базис:** " .. tostring(record.basis_name or "unknown")
  out[#out + 1] = "**Размер базиса:** " .. tostring(record.basis_size or "?")
  out[#out + 1] = ""

  if not record.ok then
    out[#out + 1] = "## Ошибка обработки"
    out[#out + 1] = tostring(record.error)
    return table.concat(out, "\n")
  end

  for _, block in ipairs(record.document.blocks or {}) do
    local rendered = render_block_single(block, record)
    if rendered and rendered ~= "" then
      out[#out + 1] = rendered
      out[#out + 1] = ""
    end
  end
  return table.concat(out, "\n")
end

local function find_block(record, id, typ)
  for _, block in ipairs(record.document.blocks or {}) do
    if block.id == id and block.type == typ then return block end
  end
end

local function first_groupable_blocks(record)
  local blocks = {}
  for _, block in ipairs(record.document.blocks or {}) do
    if block.type == "atom_table" or block.type == "bond_table" or block.type == "pivot_table" or (block.type == "table" and block.grouped) or block.type == "metric" then
      blocks[#blocks + 1] = block
    elseif block.type == "heading" or block.type == "paragraph" then
      blocks[#blocks + 1] = block
    end
  end
  return blocks
end

local function render_group(group)
  table.sort(group.records, function(a, b)
    if (a.basis_size or 0) == (b.basis_size or 0) then
      return tostring(a.basis_name or "") < tostring(b.basis_name or "")
    end
    return (a.basis_size or 0) < (b.basis_size or 0)
  end)

  local first = group.records[1]
  local out = {"# Молекула " .. group.molecule_key, ""}
  if group.warning and #group.warning > 0 then
    out[#out + 1] = "## Предупреждения"
    for _, warning in ipairs(group.warning) do out[#out + 1] = "- " .. warning end
    out[#out + 1] = ""
  end

  for _, block in ipairs(first_groupable_blocks(first)) do
    if block.type == "heading" then
      out[#out + 1] = string.rep("#", block.level or 1) .. " " .. block.text
      out[#out + 1] = ""
    elseif block.type == "paragraph" then
      out[#out + 1] = block.text
      out[#out + 1] = ""
    elseif block.type == "metric" then
      local headers = {"Показатель"}
      local row = {block.label or block.id}
      for _, record in ipairs(group.records) do
        headers[#headers + 1] = record.basis_name or basename(record.file)
        local b = find_block(record, block.id, "metric")
        row[#row + 1] = b and fmt_value(b.value) or ""
      end
      out[#out + 1] = table_to_md({headers = headers, rows = {row}})
      out[#out + 1] = ""
    elseif block.type == "atom_table" then
      local headers = {"Атом", "Z"}
      for _, record in ipairs(group.records) do headers[#headers + 1] = record.basis_name or basename(record.file) end
      local rows = {}
      for atom = 1, #(first.atomic_numbers or {}) do
        local row = {atom, first.atomic_numbers[atom] or ""}
        for _, record in ipairs(group.records) do
          local b = find_block(record, block.id, "atom_table")
          row[#row + 1] = b and fmt_value(b.values[atom]) or ""
        end
        rows[#rows + 1] = row
      end
      out[#out + 1] = "## " .. (block.title or block.id)
      out[#out + 1] = ""
      out[#out + 1] = table_to_md({headers = headers, rows = rows})
      out[#out + 1] = ""
    elseif block.type == "bond_table" then
      local keyset = {}
      for _, record in ipairs(group.records) do
        local b = find_block(record, block.id, "bond_table")
        if b then
          for key in pairs(b.values) do keyset[key] = true end
        end
      end
      local keys = {}
      for key in pairs(keyset) do keys[#keys + 1] = key end
      table.sort(keys)

      local headers = {"Связь"}
      for _, record in ipairs(group.records) do headers[#headers + 1] = record.basis_name or basename(record.file) end
      local rows = {}
      for _, key in ipairs(keys) do
        local row = {key}
        for _, record in ipairs(group.records) do
          local b = find_block(record, block.id, "bond_table")
          local value = b and b.values[key]
          row[#row + 1] = value and fmt_value(value.value) or ""
        end
        rows[#rows + 1] = row
      end
      out[#out + 1] = "## " .. (block.title or block.id)
      out[#out + 1] = ""
      out[#out + 1] = table_to_md({headers = headers, rows = rows})
      out[#out + 1] = ""
    elseif block.type == "pivot_table" then
      local blocks = {}
      for _, record in ipairs(group.records) do
        local b = find_block(record, block.id, "pivot_table")
        if b then blocks[#blocks + 1] = b end
      end
      if #blocks > 0 then
        out[#out + 1] = render_pivot_table(blocks)
        out[#out + 1] = ""
      end
    elseif block.type == "table" and block.grouped then
      out[#out + 1] = render_block_single(block, first)
      out[#out + 1] = ""
    end
  end
  return table.concat(out, "\n")
end

local function compile_blocks(code_blocks)
  local chunks = {}
  for index, code in ipairs(code_blocks or {}) do
    if type(code) == "function" then
      chunks[#chunks + 1] = {kind = "function", fn = code, index = index}
    else
      chunks[#chunks + 1] = {kind = "code", code = tostring(code or ""), index = index}
    end
  end
  return chunks
end

local function make_env(record, builder, opts)
  opts = opts or {}
  local function emit(x) builder:text(type(x) == "table" and table_to_md(x) or tostring(x or "")) end
  local function table_col(label, fn, opts) return {label, fn, opts or {}} end
  local function table_value(label, fn, opts) return {label, fn, opts or {}} end
  local function table_summary(label, rules, opts)
    return {kind = "summary", label = label, rules = rules or {}, opts = opts or {}}
  end
  local function table_delta(label, left, right, opts)
    opts = opts or {}
    return {kind = "delta", label = label, left = left, right = right, opts = opts, mode = opts.mode}
  end
  local function basis_value(mol)
    mol = mol or record.mol
    return (mol and mol.get_basis_name and mol:get_basis_name()) or record.basis_name
  end
  local function cached(fn)
    local cache_mol, cache_value
    return function(mol, ...)
      mol = mol or record.mol
      if cache_mol ~= mol then
        cache_mol = mol
        cache_value = fn(mol)
      end
      return cache_value
    end
  end

  local env = {
    mol = record.mol,
    molecule = record.mol,
    file = record.file,
    filename = record.file,
    basis = record.basis_name,
    basis_name = record.basis_name,
    basis_size = record.basis_size,
    atomic_numbers = record.atomic_numbers,
    molecule_key = record.molecule_key,
    ordered_key = record.ordered_key,
    out = builder,
    result = builder,

    chemistry = chemistry,
    qa = qa,
    Molecule = chemistry.Molecule,
    Matrix = chemistry.Matrix,
    Group = chemistry.Group,

    print = function(...)
      local parts = {}
      for i = 1, select("#", ...) do parts[i] = tostring(select(i, ...)) end
      builder:text(table.concat(parts, " "))
    end,
    emit = emit,
    h = function(level, text) builder:h(level, text) end,
    h1 = function(text) builder:h1(text) end,
    h2 = function(text) builder:h2(text) end,
    h3 = function(text) builder:h(3, text) end,
    text = function(value, opts) builder:text(value, opts) end,
    metric = function(id, label, value, opts) builder:metric(id, label, value, opts) end,
    atom_table = function(id, opts) builder:atom_table(id, opts) end,
    atom_value = function(id, atom, value, opts) builder:atom_value(id, atom, value, opts) end,
    bond_table = function(id, opts) builder:bond_table(id, opts) end,
    bond_value = function(id, a, b, value, opts) builder:bond_value(id, a, b, value, opts) end,
    table_block = function(id, opts) builder:table(id, opts) end,
    set_cell = function(id, row, col, value, opts) builder:set_cell(id, row, col, value, opts) end,
    figure = function(id, path, caption, opts) builder:figure(id, path, caption, opts) end,
    matrix_preview = function(id, matrix, opts) builder:matrix_preview(id, matrix, opts) end,
    Table = function(title, spec)
      builder.wants_grouped = true
      return builder:_add(build_pivot_table_block(title, spec or {}, record))
    end,
    col = table_col,
    val = table_value,
    summary = table_summary,
    delta = table_delta,
    atoms = function() return atom_rows(record.mol) end,
    basis = basis_value,
    cached = cached,
    atom_index = function(_, row) return row.id end,
    atom_symbol = function(_, row) return row.symbol end,
    atom_label = function(_, row) return row.label end,
    atom_z = function(_, row) return row.Z end,

    table_to_md = table_to_md,
    matrix_to_md = matrix_to_md,
    pairs = pairs, ipairs = ipairs, next = next, select = select,
    tostring = tostring, tonumber = tonumber, type = type,
    assert = assert, error = error, pcall = pcall,
    math = math, string = string, table = table,
    LOG_FILES = _LOG_FILES,
  }

  install_plugin_exports(env, {
    emit = function(value)
      builder:text(type(value) == "table" and table_to_md(value) or tostring(value or ""))
    end,
    text = function(value)
      builder:text(tostring(value or ""))
    end,
    output = function(value)
      builder.disable_grouping = true
      builder:text(type(value) == "table" and table_to_md(value) or tostring(value or ""))
    end,
  })

  if opts.unsafe == true then
    env.debug = debug
    env.io = io
    env.os = os
    env.require = require
    setmetatable(env, {__index = _G})
  end

  return env
end

local function run_one_file(chunks, file, opts)
  local record = {
    ok = false,
    file = file,
    warnings = {},
  }

  local ok_open, mol_or_err = pcall(chemistry.Molecule.load, file)
  if not ok_open then
    record.error = tostring(mol_or_err)
    return record
  end

  local mol = mol_or_err
  record.mol = mol
  record.filename = file
  record.basis_name = mol:get_basis_name() or ("basis-" .. tostring(mol:get_basis_size()))
  record.basis_size = mol:get_basis_size()
  record.num_atoms = mol:get_num_atoms()
  record.molecule_key, record.ordered_key, record.atomic_numbers = make_molecule_signature(mol)

  local builder = DocumentBuilder.new(record)
  local env = make_env(record, builder, opts)

  for _, chunk in ipairs(chunks) do
    local ok, err
    if chunk.kind == "function" then
      ok, err = pcall(chunk.fn, {mol = mol, out = builder, qa = qa, chemistry = chemistry, file = file, record = record})
    else
      local loaded
      loaded, err = load(chunk.code, "sandbox_block_" .. tostring(chunk.index), "t", env)
      if loaded then ok, err = pcall(loaded) else ok = false end
    end
    if not ok then
      record.error = tostring(err)
      record.traceback = debug.traceback(tostring(err), 2)
      break
    end
  end

  record.document = builder:document()
  record.wants_grouped = record.document.wants_grouped
  record.disable_grouping = record.document.disable_grouping
  record.ok = record.error == nil

  mol:close()
  record.mol = nil
  return record
end

local function group_records(records)
  local by_key = {}
  local groups = {}

  for _, record in ipairs(records or {}) do
    if record.ok then
      local key = record.molecule_key
      local group = by_key[key]
      if not group then
        group = {molecule_key = key, records = {}, warning = {}, ordered_key = record.ordered_key}
        by_key[key] = group
        groups[#groups + 1] = group
      elseif group.ordered_key ~= record.ordered_key then
        group.warning[#group.warning + 1] = "Одинаковый состав, но разный порядок атомов: " .. tostring(record.file)
      end
      group.records[#group.records + 1] = record
    end
  end

  table.sort(groups, function(a, b) return tostring(a.molecule_key) < tostring(b.molecule_key) end)
  return groups
end

local function render_errors(errors)
  if #errors == 0 then return nil end
  local rows = {}
  for _, record in ipairs(errors) do rows[#rows + 1] = {record.file or "", record.error or "unknown error"} end
  return "## Ошибки обработки\n\n" .. table_to_md({headers = {"Файл", "Ошибка"}, rows = rows})
end

local function render_outputs(records, groups, errors, opts)
  opts = opts or {}
  local mode = opts.mode or "auto"
  local wants_grouped = false
  local disable_grouping = false
  for _, record in ipairs(records or {}) do
    if record.ok and record.wants_grouped then wants_grouped = true end
    if record.ok and record.disable_grouping then disable_grouping = true end
  end

  local use_grouped = not disable_grouping and (mode == "grouped" or (mode == "auto" and wants_grouped))
  local outputs = {}

  if use_grouped then
    for _, group in ipairs(groups) do
      outputs[#outputs + 1] = render_group(group)
    end
  else
    for _, record in ipairs(records or {}) do
      if record.ok then outputs[#outputs + 1] = render_record(record) end
    end
  end

  local err_md = render_errors(errors)
  if err_md then outputs[#outputs + 1] = err_md end
  return outputs, use_grouped and "grouped" or "single"
end

local function run(code_blocks, log_files, opts)
  opts = opts or {}
  _LOG_FILES = log_files or {}
  local chunks = compile_blocks(code_blocks or {})
  local records = {}
  local errors = {}

  if #_LOG_FILES == 0 then
    local record = {
      ok = true,
      file = "(no file)",
      basis_name = "",
      basis_size = 0,
      num_atoms = 0,
      molecule_key = "none",
      ordered_key = "",
      atomic_numbers = {},
      document = DocumentBuilder.new({}):document(),
    }
    local builder = DocumentBuilder.new(record)
    local env = make_env(record, builder, opts)
    for _, chunk in ipairs(chunks) do
      if chunk.kind == "code" then
        local loaded, err = load(chunk.code, "sandbox_block_" .. tostring(chunk.index), "t", env)
        if loaded then
          local ok, runerr = pcall(loaded)
          if not ok then builder:text("**Ошибка выполнения:** " .. tostring(runerr)) end
        else
          builder:text("**Ошибка компиляции:** " .. tostring(err))
        end
      end
    end
    record.document = builder:document()
    records[#records + 1] = record
  else
    for _, file in ipairs(_LOG_FILES) do
      local record = run_one_file(chunks, file, opts)
      records[#records + 1] = record
      if not record.ok then errors[#errors + 1] = record end
    end
  end

  local groups = group_records(records)
  local outputs, mode = render_outputs(records, groups, errors, opts)
  return {
    ok = #errors == 0,
    mode = mode,
    records = records,
    groups = groups,
    outputs = outputs,
    errors = errors,
  }
end

local function execute(code_blocks, log_files, opts)
  return run(code_blocks, log_files, opts).outputs
end

local function group_molecules(log_files)
  local groups = {}
  for _, file in ipairs(log_files or {}) do
    local ok, mol = pcall(chemistry.Molecule.load, file)
    if ok then
      local key = make_molecule_signature(mol)
      groups[key] = groups[key] or {}
      groups[key][#groups[key] + 1] = file
      mol:close()
    end
  end
  return groups
end

local function create_base_env(opts)
  local record = {
    file = nil,
    basis_name = nil,
    basis_size = nil,
    atomic_numbers = {},
  }
  local builder = DocumentBuilder.new(record)
  local env = make_env(record, builder, opts)
  env.clear_buf = function() builder.blocks = {}; builder.block_index = {} end
  env.result = function() return render_record({ok = true, file = "(env)", basis_name = "", basis_size = 0, atomic_numbers = {}, document = builder:document()}) end
  return env
end

return {
  run = run,
  execute = execute,
  create_base_env = create_base_env,
  group_molecules = group_molecules,
  make_molecule_signature = make_molecule_signature,
  matrix_to_md = matrix_to_md,
  table_to_md = table_to_md,
  DocumentBuilder = DocumentBuilder,
  _LOG_FILES = _LOG_FILES,
}
