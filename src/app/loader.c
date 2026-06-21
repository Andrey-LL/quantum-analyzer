#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
#include "luajit.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stddef.h>

#ifdef _WIN32
#include <windows.h>
#endif

extern const unsigned char luaJIT_BC_chemistry_ffi[];
extern const size_t luaJIT_BC_chemistry_ffi_size;
extern const unsigned char luaJIT_BC_quantum_analysis[];
extern const size_t luaJIT_BC_quantum_analysis_size;
extern const unsigned char luaJIT_BC_router[];
extern const size_t luaJIT_BC_router_size;
extern const unsigned char luaJIT_BC_sandbox[];
extern const size_t luaJIT_BC_sandbox_size;

static int load_embedded(lua_State *L) {
    const char* name = lua_tostring(L, lua_upvalueindex(1));
    const unsigned char* bc = (const unsigned char*)lua_touserdata(L, lua_upvalueindex(2));
    size_t size = (size_t)lua_tointeger(L, lua_upvalueindex(3));
    if (luaL_loadbuffer(L, (const char*)bc, size, name) != LUA_OK) lua_error(L);
    lua_call(L, 0, 1);
    return 1;
}

static void reg_module(lua_State *L, const char* name, const unsigned char* bc, size_t size) {
    lua_getglobal(L, "package");
    lua_getfield(L, -1, "preload");
    lua_pushstring(L, name);
    lua_pushlightuserdata(L, (void*)bc);
    lua_pushinteger(L, size);
    lua_pushcclosure(L, load_embedded, 3);
    lua_setfield(L, -2, name);
    lua_pop(L, 2);
}

int main(int argc, char **argv) {
    lua_State *L = luaL_newstate();
    luaL_openlibs(L);

    lua_newtable(L);
    for(int i=0; i<argc; i++) {
        lua_pushstring(L, argv[i]);
        lua_rawseti(L, -2, i);
    }
    lua_setglobal(L, "arg");

    /* Включает embedded-режим LuaJIT FFI: chemistry_ffi.lua использует ffi.C. */
    lua_pushboolean(L, 1);
    lua_setglobal(L, "QUANTUM_ANALYZER_EMBEDDED");

    lua_getglobal(L, "package");
    lua_getfield(L, -1, "path");
    const char *old_path = lua_tostring(L, -1);

    char exe_path[4096];
#ifdef _WIN32
    GetModuleFileNameA(NULL, exe_path, sizeof(exe_path));
#else
    if (argv[0]) {
        if (!realpath(argv[0], exe_path)) {
            exe_path[0] = '\0';
        }
    } else {
        exe_path[0] = '\0';
    }
#endif

    if (exe_path[0] && old_path) {
        char project_root[4096];
        strncpy(project_root, exe_path, sizeof(project_root) - 1);
        project_root[sizeof(project_root) - 1] = '\0';
        for (char *p = project_root; *p; p++) {
            if (*p == '\\') *p = '/';
        }
        char *bin_pos = strstr(project_root, "/bin/");
        if (bin_pos) {
            *bin_pos = '\0';
            lua_pushfstring(L, "%s/lua/?.lua;%s/lua/?/init.lua;%s",
                           project_root, project_root, old_path);
            lua_setfield(L, -3, "path");
        }
    }
    lua_pop(L, 2);

    reg_module(L, "chemistry_ffi", luaJIT_BC_chemistry_ffi, luaJIT_BC_chemistry_ffi_size);
    reg_module(L, "quantum_analysis", luaJIT_BC_quantum_analysis, luaJIT_BC_quantum_analysis_size);
    reg_module(L, "sandbox", luaJIT_BC_sandbox, luaJIT_BC_sandbox_size);
    reg_module(L, "router", luaJIT_BC_router, luaJIT_BC_router_size);

    if (argc > 1 && (strcmp(argv[1], "--batch") == 0 ||
                     strcmp(argv[1], "--filter") == 0 ||
                     strcmp(argv[1], "--help") == 0 ||
                     strcmp(argv[1], "-h") == 0)) {
        int rc = 0;
        lua_getglobal(L, "require");
        lua_pushstring(L, "router");
        if (lua_pcall(L, 1, 0, 0) != LUA_OK) {
            fprintf(stderr, "Error: %s\n", lua_tostring(L, -1));
            rc = 1;
        }
        lua_close(L);
        return rc;
    }

    if (argc > 1) {
        if (strcmp(argv[1], "-e") == 0 && argc > 2) {
            char code[4096];
            strncpy(code, argv[2], sizeof(code) - 1);
            code[sizeof(code) - 1] = '\0';
            size_t len = strlen(code);
            if (len > 0 && code[0] == '"') {
                memmove(code, code + 1, len);
                len--;
            }
            if (len > 0 && code[len - 1] == '"') {
                code[len - 1] = '\0';
            }
            if (luaL_dostring(L, code) != LUA_OK) {
                fprintf(stderr, "Error: %s\n", lua_tostring(L, -1));
                return 1;
            }
        } else if (luaL_dofile(L, argv[1]) != LUA_OK) {
            fprintf(stderr, "Error: %s\n", lua_tostring(L, -1));
            return 1;
        }
    } else {
        printf("Quantum Analyzer - Standalone Interpreter\n");
    }

    lua_close(L);
    return 0;
}
