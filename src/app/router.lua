local arg = arg or {}
local mode = arg[1]

local function dirname(path)
    return (path or ""):gsub("\\", "/"):match("^(.*)/[^/]*$") or "."
end

local function read_all(path)
    local f = io.open(path, "rb")
    if not f then return nil end
    local data = f:read("*a")
    f:close()
    return data
end

local function write_all(path, data)
    local f = assert(io.open(path, "wb"))
    f:write(data or "")
    f:close()
end

local function ensure_dir(path)
    local is_windows = package.config:sub(1, 1) == "\\"
    if is_windows then
        os.execute('if not exist "' .. path:gsub("/", "\\") .. '" mkdir "' .. path:gsub("/", "\\") .. '"')
    else
        os.execute(string.format("mkdir -p %q", path))
    end
end

local function report_filename(index, markdown)
    local key = markdown:match("^#%s+Молекула%s+([^\n]+)")
    if key then
        key = key:gsub("[^%w]+", "_"):gsub("^_+", ""):gsub("_+$", "")
    end
    if not key or key == "" then
        return string.format("%02d_report.md", index)
    end
    return string.format("%02d_molecule_%s.md", index, key)
end

local function template_candidates(name)
    local base = dirname(arg[0])
    local filename = name:match("%.lua$") and name or (name .. ".lua")
    return {
        name,
        filename,
        base .. "/templates/" .. filename,
        base .. "/../templates/" .. filename,
    }
end

local function load_template(name)
    for _, path in ipairs(template_candidates(name)) do
        local code = read_all(path)
        if code then
            return {
                name = name,
                code_blocks = {code},
                mode = "auto",
                source = path,
            }
        end
    end

    error("Cannot load template: " .. tostring(name))
end

if mode == "--help" or mode == "-h" then
    print("Quantum Analyzer - batch processing tool")
    print()
    print("Usage:")
    print("  quantum_analyzer --batch --template <name> --files <log> [<log> ...]")
    print("  quantum_analyzer --filter [options]")
    print("  quantum_analyzer <script.lua> [args...]")
    print()
    print("Options:")
    print("  --batch, -b     Batch processing mode")
    print("  --template, -t  Analysis template (default: standard_analysis)")
    print("  --files, -f     Gaussian log/out files")
    print("  --out-dir, -o   Save each rendered report as a separate Markdown file")
    print("  --filter        Filter mode (not implemented)")
    print("  --help, -h      Show this help")
    return
end

if mode == "--batch" then
    table.remove(arg, 1)

    local template = "standard_analysis"
    local files_list = {}
    local out_dir = nil

    local i = 1
    while i <= #arg do
        local val = arg[i]:gsub('"', '')
        if val == "--template" or val == "-t" then
            if not arg[i + 1] then
                error("--template requires a value")
            end
            template = arg[i + 1]:gsub('"', '')
            i = i + 2
        elseif val == "--out-dir" or val == "-o" then
            if not arg[i + 1] then
                error("--out-dir requires a value")
            end
            out_dir = arg[i + 1]:gsub('"', '')
            i = i + 2
        elseif val == "--files" or val == "-f" then
            i = i + 1
            while i <= #arg do
                local f = arg[i]:gsub('"', '')
                if f:sub(1, 1) == "-" then
                    break
                end
                table.insert(files_list, f)
                i = i + 1
            end
        else
            i = i + 1
        end
    end

    if #files_list == 0 then
        error("--batch requires at least one file after --files")
    end

    local files_str = table.concat(files_list, ", ")
    print("Processing files: " .. files_str)
    print("Template: " .. template)

    local sandbox = require("sandbox")
    local tpl = load_template(template)

    local code_blocks = tpl.code_blocks
    if not code_blocks then
        error("Template is empty: " .. tostring(template))
    end

    local run = sandbox.run(code_blocks, files_list, {
        mode = tpl.mode or "auto",
        format = "markdown",
    })

    print()
    print(string.format("Completed: %d successful, %d failed", #run.records - #run.errors, #run.errors))
    if out_dir then
        ensure_dir(out_dir)
    end
    for i, out in ipairs(run.outputs or {}) do
        if out_dir then
            local filename = report_filename(i, out)
            local out_path = out_dir .. "/" .. filename
            write_all(out_path, out .. "\n")
            print("Saved report: " .. out_path)
        end
        print()
        print(string.rep("=", 80))
        print(out)
        print(string.rep("=", 80))
    end
    return
elseif mode == "--filter" then
    print("--filter mode not yet implemented")
    return
else
    print("Quantum Analyzer - standalone interpreter")
    print("Usage:")
    print("  quantum_analyzer.exe <script.lua> [args...]")
    print("  quantum_analyzer.exe -e 'code'")
    print("  quantum_analyzer.exe --batch --template <name> --files <log> [<log> ...]")
    print("  quantum_analyzer.exe --filter [options]")
end
