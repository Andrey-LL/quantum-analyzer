set_project("quantum_analyzer")
add_rules("mode.debug", "mode.release")
set_languages("c11", "cxx17")

option("native")
    set_default(true)
    set_description("Enable -march=native for local builds")
option_end()

option("fast_math")
    set_default(false)
    set_description("Enable -ffast-math (may affect strict IEEE behavior)")
option_end()

option("use_vcpkg")
    set_default(os.host() == "windows")
    set_showmenu(true)
    set_description("Use vcpkg dependencies on Windows/MSVC")
option_end()

if os.host() == "windows" and (not get_config("plat") or get_config("plat") == "windows") then
    set_defaultplat("windows")
    set_defaultarchs("windows|x64")
    set_toolchains("msvc")
    set_runtimes("MT")
end

-- Общие каталоги вывода для всех целей сборки.
set_targetdir("$(projectdir)/build/lib")
set_objectdir("$(projectdir)/build/obj")

-- Системные зависимости, которые реально линкуются.
if is_plat("windows") and has_config("use_vcpkg") then
    add_requires("vcpkg::openblas", {alias = "openblas", system = true})
    add_requires("vcpkg::lapack", {alias = "lapack", system = true})
    add_requires("vcpkg::eigen3", {alias = "eigen3", system = true})
    add_requires("vcpkg::boost-dynamic-bitset", {alias = "boost_dynamic_bitset", system = true})
    add_requires("vcpkg::luajit", {alias = "luajit", system = true})
elseif is_plat("mingw") then
    -- MSYS2/UCRT64 packages are wired directly in add_mingw_deps().
else
    add_requires("openblas", {system = true})
    add_requires("luajit", {system = true})
end

-- Eigen и boost::dynamic_bitset используются только как заголовочные библиотеки. Не добавляем их
-- в add_requires/add_packages, чтобы xmake не искал и не линковал boost_*.
function add_header_only_deps()
    local function add_if_dir(dir)
        if os.isdir(dir) then
            add_includedirs(dir)
        end
    end

    if is_plat("windows") then
        local mingw_prefix = os.getenv("MINGW_PREFIX")
        if mingw_prefix then
            add_if_dir(path.join(mingw_prefix, "include"))
            add_if_dir(path.join(mingw_prefix, "include", "eigen3"))
        end
        add_defines("BOOST_ALL_NO_LIB")
    else
        add_if_dir("/usr/include")
        add_if_dir("/usr/include/eigen3")
        add_if_dir("/usr/local/include")
        add_if_dir("/usr/local/include/eigen3")
    end
end

function add_vcpkg_luajit_includedirs()
    if not is_plat("windows") then
        return
    end

    local root = get_config("vcpkg") or os.getenv("VCPKG_ROOT") or os.getenv("VCPKG_INSTALLATION_ROOT") or "C:/vcpkg"
    local triplet = get_config("vcpkg_triplet") or "x64-windows-static"
    local include_root = path.join(root, "installed", triplet, "include")

    local candidates = {
        path.join(include_root, "luajit"),
        path.join(include_root, "luajit-2.1"),
        include_root,
    }

    for _, dir in ipairs(candidates) do
        if os.isdir(dir) then
            add_includedirs(dir)
        end
    end
end

function add_mingw_deps()
    local prefix = os.getenv("MINGW_PREFIX") or "/ucrt64"
    add_includedirs(
        path.join(prefix, "include"),
        path.join(prefix, "include", "eigen3"),
        path.join(prefix, "include", "luajit-2.1"),
        path.join(prefix, "include", "luajit")
    )
    add_linkdirs(path.join(prefix, "lib"))
    add_links("openblas", "lapack", "luajit-5.1")
    add_defines("BOOST_ALL_NO_LIB")
end

includes("src/lib")
includes("src/app")

function qa_project_root()
    return os.projectdir()
end

function qa_analyzer_exe()
    local is_windows_target = is_plat("windows") or is_plat("mingw")
    return path.join(qa_project_root(), "bin", is_windows_target and "quantum_analyzer.exe" or "quantum_analyzer")
end

function qa_pandoc_filter_path()
    if os.host() == "windows" then
        return path.join(os.getenv("APPDATA"), "pandoc", "filters", "gaussian_filter.lua")
    end
    return path.join(os.getenv("HOME"), ".local", "share", "pandoc", "filters", "gaussian_filter.lua")
end

function qa_install_pandoc_filter()
    local project_root = qa_project_root()
    local out = qa_pandoc_filter_path()
    local s = io.readfile(path.join(project_root, "src", "app", "share", "pandoc", "gaussian_filter.lua.in"))
    s = s:gsub("@PROJECT_ROOT@", project_root)
    s = s:gsub("@QUANTUM_ANALYZER_EXE@", qa_analyzer_exe())
    os.mkdir(path.directory(out))
    io.writefile(out, s)
    print("Installed pandoc filter: " .. out)
    print("Quantum Analyzer path: " .. qa_analyzer_exe())
end

function qa_run_lua_test(mode, file)
    local envs = {
        PROJECT_ROOT = qa_project_root(),
        QA_TEST_MODE = mode,
    }
    if mode == "library" then
        os.execv("luajit", {file}, {envs = envs})
    else
        os.execv(qa_analyzer_exe(), {file}, {envs = envs})
    end
end

-- Диагностика зависимостей без автоматической установки.
-- Запуск: xmake run deps
target("deps")
    set_kind("phony")

    on_run(function ()
        local host_windows = os.host() == "windows"

        local function have_tool(name)
            local pathenv = os.getenv("PATH") or ""
            local sep = host_windows and ";" or ":"
            local exts = host_windows and {"", ".exe", ".bat", ".cmd"} or {""}
            for dir in pathenv:gmatch("([^" .. sep .. "]+)") do
                for _, ext in ipairs(exts) do
                    if os.isfile(path.join(dir, name .. ext)) then
                        return true
                    end
                end
            end
            return false
        end

        local function have_pkg(name)
            if not have_tool("pkg-config") then
                return false
            end
            return os.execv("pkg-config", {"--exists", name}, {try = true}) == 0
        end

        local checks = {}
        if host_windows then
            checks = {
                {name = "MSVC compiler (cl)", ok = have_tool("cl")},
                {name = "vcpkg command", ok = have_tool("vcpkg")},
                {name = "luarocks", ok = have_tool("luarocks"), optional = true},
                {name = "gnuplot", ok = have_tool("gnuplot"), optional = true},
                {name = "pandoc", ok = have_tool("pandoc"), optional = true},
                {name = "ctags", ok = have_tool("ctags"), optional = true},
            }
        else
            checks = {
                {name = "pkg-config", ok = have_tool("pkg-config")},
                {name = "luajit", ok = have_tool("luajit")},
                {name = "luarocks", ok = have_tool("luarocks")},
                {name = "openblas", ok = have_pkg("openblas")},
                {name = "lapack", ok = have_pkg("lapack")},
                {name = "eigen headers", ok = os.isdir("/usr/include/eigen3/Eigen") or os.isdir("/usr/local/include/eigen3/Eigen")},
                {name = "boost headers", ok = os.isfile("/usr/include/boost/dynamic_bitset.hpp") or os.isfile("/usr/local/include/boost/dynamic_bitset.hpp")},
                {name = "gnuplot", ok = have_tool("gnuplot"), optional = true},
                {name = "pandoc", ok = have_tool("pandoc"), optional = true},
                {name = "ctags", ok = have_tool("ctags"), optional = true},
            }
        end

        local missing_required = {}
        local missing_optional = {}
        for _, check in ipairs(checks) do
            if check.ok then
                print("[deps] ok: " .. check.name)
            elseif check.optional then
                missing_optional[#missing_optional + 1] = check.name
            else
                missing_required[#missing_required + 1] = check.name
            end
        end

        local function join(list)
            return table.concat(list, ", ")
        end

        if #missing_required > 0 then
            print("[deps] missing required: " .. join(missing_required))
        end
        if #missing_optional > 0 then
            print("[deps] missing optional: " .. join(missing_optional))
        end

        if #missing_required == 0 and #missing_optional == 0 then
            print("[deps] dependencies look installed")
            return
        end

        print("")
        if host_windows then
            print("Install required vcpkg packages manually:")
            print("    vcpkg install openblas:x64-windows-static lapack:x64-windows-static eigen3:x64-windows-static boost-dynamic-bitset:x64-windows-static luajit:x64-windows-static")
            print("")
            print("Then build with:")
            print("    xmake")
        else
            print("Install manually, system packages:")
            print("  Ubuntu/Debian:")
            print("    sudo apt update")
            print("    sudo apt install -y build-essential pkg-config libopenblas-dev liblapack-dev libeigen3-dev libboost-dev luajit luarocks gnuplot pandoc universal-ctags")
            print("")
            print("  Fedora:")
            print("    sudo dnf install -y gcc gcc-c++ make pkgconf-pkg-config openblas-devel lapack-devel eigen3-devel boost-devel luajit luarocks gnuplot pandoc ctags")
            print("")
            print("  Arch Linux:")
            print("    sudo pacman -S --needed base-devel pkgconf openblas lapack eigen boost luajit luarocks gnuplot pandoc ctags")
        end

        if #missing_required > 0 then
            raise("[deps] install required dependencies and run xmake run deps again")
        end
    end)

-- Установка фильтра pandoc для локальной интеграции отчётов.
target("pandoc_integration")
    set_kind("phony")
    on_run(function ()
        local root = os.projectdir()
        local is_windows_target = is_plat("windows") or is_plat("mingw")
        local exe = path.join(root, "bin", is_windows_target and "quantum_analyzer.exe" or "quantum_analyzer")
        local out = os.host() == "windows"
            and path.join(os.getenv("APPDATA"), "pandoc", "filters", "gaussian_filter.lua")
            or path.join(os.getenv("HOME"), ".local", "share", "pandoc", "filters", "gaussian_filter.lua")
        local s = io.readfile(path.join(root, "src", "app", "share", "pandoc", "gaussian_filter.lua.in"))
        s = s:gsub("@PROJECT_ROOT@", root)
        s = s:gsub("@QUANTUM_ANALYZER_EXE@", exe)
        os.mkdir(path.directory(out))
        io.writefile(out, s)
        print("Installed pandoc filter: " .. out)
        print("Quantum Analyzer path: " .. exe)
    end)

-- Генерация публичных примеров без shell/perl-склейки.
-- Запуск: xmake run examples
target("examples")
    set_kind("phony")
    add_deps("quantum_analyzer")

    on_run(function ()
        local root = os.projectdir()
        local exe = path.join(root, "bin", (is_plat("windows") or is_plat("mingw")) and "quantum_analyzer.exe" or "quantum_analyzer")
        local generated = path.join(root, "examples", "generated")
        local reports_dir = path.join(generated, "analysis_reports")
        os.mkdir(generated)
        os.rm(reports_dir)
        os.mkdir(reports_dir)

        os.execv(exe, {
            "--batch",
            "--template", "examples/analysis_demo.lua",
            "--out-dir", "examples/generated/analysis_reports",
            "--files",
            "examples/fixtures/methane_6-31g.log",
            "examples/fixtures/methane_sto-3g.log",
            "examples/fixtures/water_sto-3g.log",
        })
        local reports = os.files(path.join(reports_dir, "*.md"))
        table.sort(reports)
        assert(#reports == 2, "analysis_demo must generate exactly two grouped report files")
        print("Generated grouped analysis reports:")
        for _, report in ipairs(reports) do
            print("  " .. report)
        end

        local filter = os.host() == "windows"
            and path.join(os.getenv("APPDATA"), "pandoc", "filters", "gaussian_filter.lua")
            or path.join(os.getenv("HOME"), ".local", "share", "pandoc", "filters", "gaussian_filter.lua")
        local filter_tpl = io.readfile(path.join(root, "src", "app", "share", "pandoc", "gaussian_filter.lua.in"))
        filter_tpl = filter_tpl:gsub("@PROJECT_ROOT@", root)
        filter_tpl = filter_tpl:gsub("@QUANTUM_ANALYZER_EXE@", exe)
        os.mkdir(path.directory(filter))
        io.writefile(filter, filter_tpl)
        print("Installed pandoc filter: " .. filter)

        os.execv("pandoc", {
            "examples/pandoc_demo.md",
            "--lua-filter", filter,
            "-s",
            "-o", "examples/generated/pandoc_demo.html",
        })
        print("Generated examples/generated/pandoc_demo.html")
    end)

target("test-library")
    set_kind("phony")
    add_deps("lib_shared")

    on_run(function ()
        local root = os.projectdir()
        local envs = {PROJECT_ROOT = root, QA_TEST_MODE = "library"}
        for _, f in ipairs(os.files(path.join(root, "tests", "lib", "test_*.lua"))) do
            os.execv("luajit", {f}, {envs = envs})
        end
        os.execv("luajit", {path.join(root, "tests", "app", "test_sandbox.lua")}, {envs = envs})
        os.execv("luajit", {path.join(root, "tests", "app", "test_viz_gnuplot.lua")}, {envs = envs})
    end)

target("test-app")
    set_kind("phony")
    add_deps("quantum_analyzer")

    on_run(function ()
        local root = os.projectdir()
        local exe = path.join(root, "bin", (is_plat("windows") or is_plat("mingw")) and "quantum_analyzer.exe" or "quantum_analyzer")
        local envs = {
            PROJECT_ROOT = root,
            QA_TEST_MODE = "app",
            QA_TEST_WINDOWS = (is_plat("windows") or is_plat("mingw")) and "1" or "0",
        }
        for _, f in ipairs(os.files(path.join(root, "tests", "lib", "test_*.lua"))) do
            os.execv(exe, {f}, {envs = envs})
        end
        os.execv(exe, {path.join(root, "tests", "app", "test_sandbox.lua")}, {envs = envs})
        os.execv(exe, {path.join(root, "tests", "app", "test_viz_gnuplot.lua")}, {envs = envs})
        os.execv(exe, {path.join(root, "tests", "app", "test_router.lua")}, {envs = envs})
    end)

target("test-release-app")
    set_kind("phony")
    add_deps("quantum_analyzer")

    on_run(function ()
        local root = os.projectdir()
        local exe = path.join(root, "bin", (is_plat("windows") or is_plat("mingw")) and "quantum_analyzer.exe" or "quantum_analyzer")
        local envs = {
            PROJECT_ROOT = root,
            QA_TEST_MODE = "app",
            QA_TEST_WINDOWS = (is_plat("windows") or is_plat("mingw")) and "1" or "0",
        }
        for _, f in ipairs(os.files(path.join(root, "tests", "lib", "test_*.lua"))) do
            os.execv(exe, {f}, {envs = envs})
        end
        os.execv(exe, {path.join(root, "tests", "app", "test_router.lua")}, {envs = envs})
    end)

target("test-luarocks")
    set_kind("phony")
    add_deps("lib_shared")

    on_run(function ()
        local root = os.projectdir()
        if is_plat("windows") then
            os.execv("cmd", {"/c", path.join(root, "tests", "run_luarocks.cmd")})
        else
            os.execv("bash", {path.join(root, "tests", "run_luarocks.sh")})
        end
    end)

target("test-pandoc")
    set_kind("phony")
    add_deps("quantum_analyzer")

    on_run(function ()
        local root = os.projectdir()
        local exe = path.join(root, "bin", (is_plat("windows") or is_plat("mingw")) and "quantum_analyzer.exe" or "quantum_analyzer")
        local filter = os.host() == "windows"
            and path.join(os.getenv("APPDATA"), "pandoc", "filters", "gaussian_filter.lua")
            or path.join(os.getenv("HOME"), ".local", "share", "pandoc", "filters", "gaussian_filter.lua")
        local s = io.readfile(path.join(root, "src", "app", "share", "pandoc", "gaussian_filter.lua.in"))
        s = s:gsub("@PROJECT_ROOT@", root)
        s = s:gsub("@QUANTUM_ANALYZER_EXE@", exe)
        os.mkdir(path.directory(filter))
        io.writefile(filter, s)
        print("Installed pandoc filter: " .. filter)

        os.execv("pandoc", {
            "examples/pandoc_demo.md",
            "--lua-filter", filter,
            "-s",
            "-o", "examples/generated/pandoc_demo.html",
        })

        local html = io.readfile(path.join(root, "examples", "generated", "pandoc_demo.html")) or ""
        assert(html:find("pandoc_overlap_atoms_heatmap.svg", 1, true), "Pandoc demo missing overlap heatmap")
        assert(html:find("pandoc_mulliken_charges.svg", 1, true), "Pandoc demo missing charge bar chart")
        assert(not html:find("Ошибка", 1, true), "Pandoc demo contains execution error")
        print("Pandoc integration ok")
    end)

target("test-all")
    set_kind("phony")

    on_run(function ()
        os.execv("xmake", {"run", "test-luarocks"})
        os.execv("xmake", {"run", "test-library"})
        os.execv("xmake", {"run", "test-app"})
        os.execv("xmake", {"run", "test-pandoc"})
    end)

target("package-release")
    set_kind("phony")
    add_deps("quantum_analyzer")

    on_run(function ()
        local root = os.projectdir()
        local version = "1.0"
        local platform
        if is_plat("windows") or is_plat("mingw") then
            platform = "windows-x86_64"
        elseif is_plat("linux") then
            platform = "linux-x86_64"
        else
            platform = (os.host() or "unknown") .. "-" .. (os.arch() or "unknown")
        end

        local package_name = "quantum-analyzer-" .. version .. "-" .. platform
        local dist_dir = path.join(root, "dist")
        local stage = path.join(dist_dir, package_name)
        local is_windows_target = is_plat("windows") or is_plat("mingw")
        local exe = path.join(root, "bin", is_windows_target and "quantum_analyzer.exe" or "quantum_analyzer")

        os.mkdir(dist_dir)
        os.rm(stage)
        os.mkdir(stage)

        os.cp(path.join(root, "README.md"), stage)
        os.cp(path.join(root, "LICENSE"), stage)
        os.cp(exe, stage)
        if is_plat("mingw") then
            local copied = {}
            local function copy_runtime_dll(dll)
                if dll
                    and os.isfile(dll)
                    and not dll:match("[/\\]Windows[/\\]")
                    and not dll:match("libquantum_analyzer_core%.dll$")
                    and not copied[dll]
                then
                    os.cp(dll, stage)
                    copied[dll] = true
                end
            end

            local ldd_output = os.iorunv("ldd", {exe}) or ""
            for line in ldd_output:gmatch("[^\r\n]+") do
                copy_runtime_dll(line:match("=>%s+([^%s]+%.dll)") or line:match("^%s*([^%s]+%.dll)"))
            end

            local prefix = os.getenv("MINGW_PREFIX") or "/ucrt64"
            local bin_dir = path.join(prefix, "bin")
            local patterns = {
                "libgcc_s_*.dll",
                "libstdc++-6.dll",
                "libwinpthread-1.dll",
                "libgfortran-*.dll",
                "libquadmath-*.dll",
                "libgomp-1.dll",
                "libopenblas*.dll",
                "libblas*.dll",
                "liblapack*.dll",
                "libluajit-*.dll",
                "lua51.dll",
            }
            for _, pattern in ipairs(patterns) do
                for _, dll in ipairs(os.files(path.join(bin_dir, pattern))) do
                    copy_runtime_dll(dll)
                end
            end
        end
        os.cp(path.join(root, "bin", "plugins"), stage)
        os.cp(path.join(root, "bin", "templates"), stage)

        local examples_dst = path.join(stage, "examples")
        os.mkdir(examples_dst)
        os.cp(path.join(root, "examples", "README.md"), examples_dst)
        os.cp(path.join(root, "examples", "fixtures"), examples_dst)

        local archive
        if is_plat("windows") or is_plat("mingw") then
            archive = path.join(dist_dir, package_name .. ".zip")
            os.rm(archive)
            local function ps_quote(s)
                return "'" .. tostring(s):gsub("'", "''") .. "'"
            end
            os.execv("powershell", {
                "-NoProfile",
                "-ExecutionPolicy", "Bypass",
                "-Command",
                "Compress-Archive -Path " .. ps_quote(stage) ..
                    " -DestinationPath " .. ps_quote(archive) .. " -Force",
            })
        else
            archive = path.join(dist_dir, package_name .. ".tar.gz")
            os.rm(archive)
            os.execv("tar", {"-czf", archive, "-C", dist_dir, package_name})
        end

        print("Release package: " .. archive)
    end)
