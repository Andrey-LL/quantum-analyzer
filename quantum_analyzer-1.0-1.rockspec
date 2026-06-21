package = "quantum_analyzer"
version = "1.0-1"
rockspec_format = "3.0"

source = {
    url = "file://."
}

description = {
    summary = "Quantum analysis library for Gaussian file processing",
    homepage = "https://github.com/Andrey-LL/quantum-analyzer",
    license = "MIT",
    labels = {"quantum", "chemistry", "gaussian", "analysis"}
}

dependencies = {
    "lua >= 5.1",
}

build = {
    type = "command",

    build_command = [[
        xmake f -m release --confirm=n
        xmake build lib_shared
    ]],

    install_command = [[
        cp -f build/lib/libquantum_analyzer_core.so "$(LIBDIR)/"
        cp -f src/lib/lua_core/chemistry_ffi.lua "$(LUADIR)/"
        cp -f src/lib/lua_core/quantum_analysis.lua "$(LUADIR)/"
    ]],

    platforms = {
        windows = {
            build_command = [[
                xmake f -m release --confirm=n -p windows
                xmake build lib_shared
            ]],
            install_command = [[
                copy /Y build\lib\quantum_analyzer_core.dll "$(LIBDIR)\"
                copy /Y src\lib\lua_core\chemistry_ffi.lua "$(LUADIR)\"
                copy /Y src\lib\lua_core\quantum_analysis.lua "$(LUADIR)\"
            ]]
        }
    }
}
