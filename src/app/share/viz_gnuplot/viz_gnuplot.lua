-- Модуль визуализации через gnuplot

local bit = require("bit")

local M = {}

-- Пути вычисляются из расположения текущего модуля.
local function get_project_root()
  local info = debug.getinfo(1, "S")
  local source = info.source
  if source:sub(1, 1) == "@" then source = source:sub(2) end
  return source:match("^(.*)/plugins/viz_gnuplot/viz_gnuplot%.lua$")
      or source:match("^(.*)/src/app/share/viz_gnuplot/viz_gnuplot%.lua$")
      or "."
end

local MODULE_SOURCE = (debug.getinfo(1, "S").source or ""):gsub("^@", "")
local PLUGIN_DIR = MODULE_SOURCE:gsub("/[^/]+$", "")
local PROJECT_ROOT = get_project_root()

-- Хэш FNV-1a 32-bit используется для имён файлов кэша.
local function fnv1a32_init()
  return 2166136261
end

local function fnv1a32_update(hash, s)
  local prime = 16777619
  for i = 1, #s do
    hash = bit.bxor(hash, s:byte(i))
    hash = bit.band(hash * prime, 0xFFFFFFFF)
  end
  return hash
end

local function fnv1a32_hex(hash)
  return string.format("%08x", hash)
end

local function file_exists(p)
  local f = io.open(p, "rb")
  if f then f:close(); return true end
  return false
end

local function ensure_dir(dir)
  if package.config:sub(1, 1) == "\\" then
    os.execute(string.format('if not exist %q mkdir %q', dir, dir))
  else
    os.execute(string.format("mkdir -p %q", dir))
  end
end

local function parent_dir(path)
  return tostring(path or ""):match("^(.*)/[^/]+$")
end

local function run_gnuplot_script(script_text)
  local gpfile = os.tmpname() .. ".gp"
  local f = assert(io.open(gpfile, "w"))
  f:write(script_text)
  f:close()
  -- %q корректно обрабатывает пробелы в путях, но не является защитой для недоверенного ввода.
  local ok = os.execute(string.format("gnuplot %q", gpfile))
  os.remove(gpfile)
  return ok == true or ok == 0
end

local function read_text_file(path)
  local f = io.open(path, "rb")
  if not f then
    return nil
  end
  local s = f:read("*a")
  f:close()
  return s
end

local function write_text_file(path, content)
  local f = io.open(path, "wb")
  if not f then
    return false
  end
  f:write(content)
  f:close()
  return true
end

local function normalize_svg_tmp_paths(svg_path)
  local content = read_text_file(svg_path)
  if not content then
    return
  end
  content = content:gsub("/tmp/lua_[%w_]+%.dat", "@DATA@")
  content = content:gsub("\\\\[^\"<>]-lua_[%w_]+%.dat", "@DATA@")
  write_text_file(svg_path, content)
end

local function read_all(path)
  local f = assert(io.open(path, "rb"))
  local s = f:read("*a")
  f:close()
  return s
end

-- Каталог кэша SVG/PNG для gnuplot.
local function cache_dir()
    return PROJECT_ROOT .. "/.cache/gaussian"
end

local function apply_placeholders(tpl, dict)
  return (tpl:gsub("@([A-Z0-9_]+)@", function(key)
    local v = dict[key]
    if v == nil then return "" end
    return tostring(v)
  end))
end

function M.render(template_name, write_dat_fn, opts, extra)
    opts = opts or {}
    extra = extra or {}

    local cdir = opts.cache_dir or cache_dir()
    ensure_dir(cdir)

    local dat = os.tmpname() .. ".dat"
    local f = assert(io.open(dat, "w"))

    local h = fnv1a32_init()
    h = fnv1a32_update(h, template_name .. "\n")
    h = fnv1a32_update(h, (opts.caption or "") .. "\n")

    -- Данные графика участвуют в хэше кэша.
    h = write_dat_fn(f, h, fnv1a32_update)
    f:close()

    -- Дополнительные параметры, влияющие на изображение, тоже входят в хэш.
    for k, v in pairs(extra) do
        h = fnv1a32_update(h, k .. "=" .. tostring(v) .. "\n")
    end

    local hex = fnv1a32_hex(h)
    local out = opts.output or (cdir .. "/" .. template_name .. "_" .. hex .. ".svg")
    local link = opts.link or out
    local out_dir = parent_dir(out)
    if out_dir and out_dir ~= "" then
        ensure_dir(out_dir)
    end

    if file_exists(out) and not opts.overwrite then
        os.remove(dat)
        return string.format("![%s](%s)", opts.caption or "", link)
    end

    local tpl_path = PLUGIN_DIR .. "/" .. template_name .. ".gp"
    local tpl = read_all(tpl_path)

    local ph = {
        DATA = dat,
        OUTPUT = out,
        CAPTION = opts.caption or "",
    }
    -- Дополнительные плейсхолдеры передаются в шаблон gnuplot.
    for k, v in pairs(extra) do ph[k] = v end

    tpl = apply_placeholders(tpl, ph)

    assert(run_gnuplot_script(tpl), "gnuplot failed: " .. template_name)
    os.remove(dat)

    if out:sub(-4):lower() == ".svg" then
      normalize_svg_tmp_paths(out)
    end

    return string.format("![%s](%s)", opts.caption or "", link)
end

local atom_symbols = {
  [1]="H", [2]="He", [6]="C", [7]="N", [8]="O", [9]="F", [15]="P", [16]="S", [17]="Cl"
}
local function get_symbol(z)
  return atom_symbols[z] or ("Z"..z)
end


local function atom_boundaries(ao2atom)
    local bounds = {}
    for i = 2, #ao2atom do
        if ao2atom[i] ~= ao2atom[i - 1] then
            bounds[#bounds + 1] = i - 0.5
        end
    end
    return bounds
end

local function atom_labels(ao2atom, molecule)
  local labels = {}
  local start = 1
  for i = 2, #ao2atom + 1 do
    if i > #ao2atom or ao2atom[i] ~= ao2atom[i-1] then
      local finish = i - 1
      local center = (start + finish) / 2 -- Центр блока атома в 1-based координатах.

      local atom_idx = ao2atom[start] -- Индекс атома 1-based.

      local name = get_symbol(molecule:get_atomic_number(atom_idx))

      table.insert(labels, string.format('set label "%s" at first %g, graph 0 offset 0,-1.5 center font "Arial-Bold,16"', name, center))

      table.insert(labels, string.format('set label "%s" at graph 0, first %g offset -1.5,0 right font "Arial-Bold,16"', name, center))

      start = i
    end
  end
  return table.concat(labels, "\n")
end


local function gnuplot_grid(bounds, nrows, ncols)
  local t = {}
  for _, b in ipairs(bounds) do
    -- Вертикальная линия границы атомного блока.
    t[#t+1] = string.format("set arrow from %g, 0.5 to %g, %g nohead lw 2 lc rgb 'black' front", b, b, nrows + 0.5)
    -- Горизонтальная линия границы атомного блока.
    t[#t+1] = string.format("set arrow from 0.5, %g to %g, %g nohead lw 2 lc rgb 'black' front", b, ncols + 0.5, b)
  end
  return table.concat(t, "\n")
end

function M.heatmap_atoms(matrix, ao2atom, opts)
  local rows, cols = matrix:size()
  local bounds = atom_boundaries(ao2atom)
  local grid = gnuplot_grid(bounds, rows, cols)
  local mol = matrix.molecule
  local labels_cmd = ""
  if mol then
    labels_cmd = atom_labels(ao2atom, mol)
  end
  local ext = (rows > 100 or cols > 100) and "png" or "svg"

  return M.render("heatmap_atoms", function(f, h, upd)
    h = upd(h, string.format("%d %d\n", rows, cols))
    for i=1,rows do
      for j=1,cols do
        local line = string.format("%d %d %.17g\n", i, j, matrix:get(i,j))
        f:write(line); h = upd(h, line)
      end
      f:write("\n"); h = upd(h, "\n")
    end
    return h

  end, opts, { GRID = grid , LABELS = labels_cmd, FORMAT = ext})
end


function M.heatmap(matrix, opts)
  local rows, cols = matrix:size()
  return M.render("heatmap", function(f, h, upd)
    h = upd(h, string.format("%d %d\n", rows, cols))
    for i = 1, rows do
      for j = 1, cols do
        local line = string.format("%d %d %.17g\n", i, j, matrix:get(i, j))
        f:write(line)
        h = upd(h, line)
      end
      f:write("\n")
      h = upd(h, "\n")
    end
    return h
  end, opts, {})
end

function M.bar(values, opts)
  opts = opts or {}
  return M.render("bar", function(f, h, upd)
    for i, item in ipairs(values or {}) do
      local label = item.label or item[1] or tostring(i)
      local value = item.value or item[2] or 0
      local line = string.format("%s %.17g\n", tostring(label), tonumber(value) or 0)
      f:write(line)
      h = upd(h, line)
    end
    return h
  end, opts, {
    XLABEL = opts.xlabel or "",
    YLABEL = opts.ylabel or "",
  })
end

function M.install(env, ctx)
  local function output(name, fn)
    env[name] = function(...)
      local markdown = fn(...)
      if markdown ~= nil and ctx then
        if ctx.output then
          ctx.output(markdown)
        elseif ctx.emit then
          ctx.emit(markdown)
        end
      end
      return markdown
    end
  end

  output("heatmap_atoms", M.heatmap_atoms)
  output("heatmap", M.heatmap)
  output("bar", M.bar)
end

return M
