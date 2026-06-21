-- Общая конфигурация статического и динамического вариантов quantum_analyzer_core.
function configure_library(target_name, kind)
    target(target_name)
    set_kind(kind)
    set_basename("quantum_analyzer_core")
    -- Целевой каталог задан в корневом xmake.lua.

    add_files("core/*.cpp")
    if is_plat("windows") and has_config("use_vcpkg") then
        add_packages("openblas", "lapack", "eigen3", "boost_dynamic_bitset", "luajit")
    elseif is_plat("mingw") then
        add_mingw_deps()
    else
        add_packages("luajit", "openblas")
        add_header_only_deps()
    end
    add_defines("QUANTUM_ANALYZER_API_EXPORT")

    add_headerfiles("core/api.h", {install = true})
    add_headerfiles("core/internal.h", {install = false})

    if is_plat("windows") then
        set_runtimes("MT")
        add_cxflags("/O2", "/fp:precise", "/EHsc", "/W3")
        add_defines("_CRT_SECURE_NO_WARNINGS")
    elseif is_plat("mingw") then
        add_cxflags("-O3", "-fPIC", "-fvisibility=default")
    else
        add_cxflags("-O3", "-fPIC", "-fvisibility=default")
        if has_config("native") then
            add_cxflags("-march=native")
        end
        if has_config("fast_math") then
            add_cxflags("-ffast-math")
        end
    end
end

-- Статическое ядро для встроенного приложения.
configure_library("lib", "static")

-- Динамическое ядро для LuaJIT/LuaRocks режима.
configure_library("lib_shared", "shared")

-- Настройки динамической библиотеки.
target("lib_shared")
    if is_plat("linux") then
        add_ldflags("-Wl,-rpath,$ORIGIN", "-lstdc++")
    end
    if is_mode("release") then
        set_strip("all")
    end
