target("quantum_analyzer")
    set_kind("binary")
    set_targetdir("$(projectdir)/bin")
    set_basename("quantum_analyzer")
    set_strip("none")

    -- Linux/MinGW need whole-archive linking below so embedded LuaJIT can
    -- resolve C API symbols through ffi.C. Keep lib as a build dependency
    -- there, but do not let Xmake link the static archive a second time.
    if is_plat("linux") or is_plat("mingw") then
        add_deps("lib", {links = false})
    else
        add_deps("lib")
    end
    add_linkdirs("$(projectdir)/build/lib")

    if is_plat("windows") and has_config("use_vcpkg") then
        add_packages("openblas", "lapack", "eigen3", "boost_dynamic_bitset", "luajit")
        add_vcpkg_luajit_includedirs()
    elseif is_plat("mingw") then
        add_mingw_deps()
    else
        add_packages("luajit", "openblas")
        add_header_only_deps()
    end
    set_policy("check.auto_ignore_flags", false)

    add_files("loader.c")

    -- add_deps обеспечивает линковку с target("lib").

    if is_plat("windows") then
        set_runtimes("MT")
        add_cxflags("/O2", "/fp:precise", "/EHsc", "/W3")
        add_defines("_CRT_SECURE_NO_WARNINGS")
        add_defines("GAUSSIAN_API=__declspec(dllexport)")
        add_defines("QUANTUM_ANALYZER_EMBEDDED")
    elseif is_plat("mingw") then
        add_cxflags("-O3", "-fPIC", "-Wall")
        add_defines("GAUSSIAN_API=__declspec(dllexport)")
        add_defines("QUANTUM_ANALYZER_EMBEDDED")
    else
        add_cxflags("-O3", "-fPIC", "-Wall")
        if has_config("native") then
            add_cxflags("-march=native")
        end
        if has_config("fast_math") then
            add_cxflags("-ffast-math")
        end
        add_defines("GAUSSIAN_API=__attribute__((visibility(\"default\")))")
        add_defines("QUANTUM_ANALYZER_EMBEDDED")
    end

    on_load(function (target)
        local build_dir = path.join(os.projectdir(), "build")
        local app_dir   = path.join(os.projectdir(), "src", "app")
        local lua_dir   = path.join(os.projectdir(), "src", "lib", "lua_core")
        local luajit_prog = "luajit"
        os.mkdir(build_dir)

        -- Lua-модули компилируются LuaJIT в bytecode и встраиваются как C-массивы.
        -- Анализ выполняется через sandbox: шаблон передаётся как кодовый блок.
        local modules = {
            {name = "chemistry_ffi",       src = path.join(lua_dir, "chemistry_ffi.lua")},
            {name = "quantum_analysis",    src = path.join(lua_dir, "quantum_analysis.lua")},
            {name = "sandbox",             src = path.join(app_dir, "sandbox.lua")},
            {name = "router",              src = path.join(app_dir, "router.lua")},
        }

        for _, mod in ipairs(modules) do
            local dst_c = path.join(build_dir, mod.name .. ".c")

            -- Компилируем исходный Lua-модуль напрямую, без промежуточных patched-файлов.
            os.execv(luajit_prog, {"-b", "-n", mod.name, mod.src, dst_c})

            -- LuaJIT bytecode C-файл дополняется stddef.h для size_t.
            local c_content = io.readfile(dst_c)
            if not c_content:match("#include%s+<stddef.h>") then
                c_content = c_content:gsub("(#ifdef __cplusplus\nextern \"C\"\n#endif)", "#include <stddef.h>\n%1")
            end

            -- Размер bytecode нужен loader.c для luaL_loadbuffer.
            c_content = c_content:gsub(
                "(const unsigned char luaJIT_BC_" .. mod.name .. "%[%] = %{[^}]+%});",
                "%1;\nconst size_t luaJIT_BC_" .. mod.name .. "_size = sizeof(luaJIT_BC_" .. mod.name .. ");"
            )
            io.writefile(dst_c, c_content)
            target:add("files", dst_c)
        end

        local api_h_path = path.join(os.projectdir(), "src", "lib", "core", "api.h")
        local api_h = io.readfile(api_h_path)
        local api_exports = {}
        for line in (api_h or ""):gmatch("[^\r\n]+") do
            line = line:gsub("//.*$", "")
            local func = line:match("QUANTUM_ANALYZER_API%s+.-%s+([%w_]+)%s*%(")
            if func then
                api_exports[#api_exports + 1] = func
            end
        end
        table.sort(api_exports)

        -- Экспорт C API для LuaJIT FFI: .def на Windows, динамическая таблица символов на Linux.
        if is_plat("windows") or is_plat("mingw") then
            local def_file = path.join(build_dir, "logic.def")
            local def_content = {"EXPORTS"}
            for _, func in ipairs(api_exports) do
                def_content[#def_content + 1] = "    " .. func
            end
            io.writefile(def_file, table.concat(def_content, "\n"))
            target:add("files", def_file)
            if is_plat("mingw") then
                for _, func in ipairs(api_exports) do
                    target:add("ldflags", "-Wl,-u," .. func)
                end
                target:add("ldflags",
                    "-Wl,--export-all-symbols",
                    "-Wl,--whole-archive",
                    "build/lib/libquantum_analyzer_core.a",
                    "-Wl,--no-whole-archive",
                    "-lopenblas",
                    "-llapack",
                    "-lluajit-5.1",
                    {force = true}
                )
            end
        elseif is_plat("linux") then
            for _, func in ipairs(api_exports) do
                target:add("ldflags", "-Wl,-u," .. func)
                target:add("ldflags", "-Wl,--export-dynamic-symbol=" .. func)
            end
            target:add("ldflags",
                "-static-libgcc",
                "-static-libstdc++",
                "-Wl,--whole-archive",
                "build/lib/libquantum_analyzer_core.a",
                "-Wl,--no-whole-archive",
                "-lopenblas",
                {force = true}
            )
        end
    end)

    -- Ресурсы времени выполнения копируются рядом со встроенным бинарником.
    after_build(function (target)
        local targetdir = target:targetdir()
        local plugins_src = path.join(os.projectdir(), "src", "app", "share")
        local plugins_dst = path.join(targetdir, "plugins")
        local templates_src = path.join(os.projectdir(), "src", "app", "templates")
        local templates_dst = path.join(targetdir, "templates")
        os.rm(plugins_dst)
        os.rm(templates_dst)
        os.mkdir(plugins_dst)
        os.mkdir(templates_dst)
        os.cp(path.join(plugins_src, "*"), plugins_dst)
        os.cp(path.join(templates_src, "*"), templates_dst)
        print("Copied plugins/ to " .. plugins_dst)
        print("Copied templates/ to " .. templates_dst)
    end)
