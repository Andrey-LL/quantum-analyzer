# Quantum Analyzer

Quantum Analyzer is an educational C++ / LuaJIT project for reading Gaussian
output files and running small quantum-chemistry analysis workflows. The project
is prepared as a public demo for a junior C++ developer portfolio: it shows a
C++ computational core, a C API boundary, LuaJIT FFI bindings, and an Xmake
build setup.

## What it does

The program reads Gaussian `.log` / `.out` files, extracts molecular metadata
and matrices, and exposes numerical primitives to Lua analysis scripts.

Current capabilities include:

- opening Gaussian output files through a C API;
- reading basis size, atoms, electron counts, method and basis metadata;
- accessing matrices such as overlap and density matrices;
- running matrix operations and quantum-analysis helpers;
- generating grouped Markdown tables from analysis templates;
- generating compact one-molecule summaries and method-comparison reports;
- using the same C++ core from a standalone embedded app or from LuaJIT through
  a shared library.

## Architecture

The project has three main layers:

- `src/lib/core/` - C++17 core with parsing, matrix operations and numerical
  analysis primitives.
- `src/lib/lua_core/` - LuaJIT FFI wrappers around the exported C API.
- `src/app/` - standalone LuaJIT-based application that embeds Lua modules and
  runs analysis templates against Gaussian files.

Architecture overview is described directly in this README and source tree layout below.

## Tech stack

- C++17
- C API for ABI-stable LuaJIT FFI access
- LuaJIT and LuaJIT FFI
- Xmake
- OpenBLAS / LAPACK
- Eigen headers
- Boost headers, mainly `boost::dynamic_bitset`
- Cross-platform test scripts (`.sh` and `.cmd`)

## Build

Install required system dependencies first.

Linux (Debian/Ubuntu):

```bash
sudo apt update
sudo apt install -y build-essential pkg-config libopenblas-dev liblapack-dev libeigen3-dev libboost-dev luajit luarocks gnuplot pandoc universal-ctags
```

Check project dependencies:

```bash
xmake run deps
```

Windows (MSVC + vcpkg):

```powershell
vcpkg install openblas:x64-windows-static lapack:x64-windows-static eigen3:x64-windows-static boost-dynamic-bitset:x64-windows-static luajit:x64-windows-static
xmake
```

Windows prerequisites:

- Visual Studio Build Tools (MSVC x64 toolchain)
- `xmake` in `PATH`
- `vcpkg` in `PATH` and configured (`VCPKG_ROOT`)

Build the static C++ core:

```bash
xmake --build lib
```

Build the shared C++ core for LuaJIT dynamic loading:

```bash
xmake --build lib_shared
```

Build the standalone application:

```bash
xmake --build quantum_analyzer
```

Or build the default target:

```bash
xmake
```

## Run

Run a simple embedded LuaJIT command:

```bash
bin/quantum_analyzer -e 'print("Quantum Analyzer is running")'
```

Run a minimal embedded FFI smoke check:

```bash
bin/quantum_analyzer -e 'local ffi=require("ffi"); ffi.cdef[[typedef void* GaussianFileHandle; GaussianFileHandle gaussian_open(const char* filename); void gaussian_close(GaussianFileHandle file);]]; local h=ffi.C.gaussian_open("examples/fixtures/methane_6-31g.log"); assert(h ~= nil); ffi.C.gaussian_close(h); print("embedded ffi.C api ok")'
```

Run a shared-library FFI smoke check:

```bash
luajit -e 'local ffi=require("ffi"); ffi.cdef[[typedef void* GaussianFileHandle; GaussianFileHandle gaussian_open(const char* filename); void gaussian_close(GaussianFileHandle file);]]; local lib=ffi.load("build/lib/libquantum_analyzer_core.so"); local h=lib.gaussian_open("examples/fixtures/methane_6-31g.log"); assert(h ~= nil); lib.gaussian_close(h); print("shared ffi.load api ok")'
```

Run integration tests (platform auto-detected):

```bash
xmake run test-all
```

Run LuaRocks install smoke only:

```bash
xmake run test-luarocks
```

## Analysis templates

The standalone application can run built-in Lua analysis templates:

```bash
bin/quantum_analyzer --batch --template standard_analysis --files examples/fixtures/methane_6-31g.log
bin/quantum_analyzer --batch --template summary_analysis --files examples/fixtures/methane_6-31g.log examples/fixtures/methane_sto-3g.log
bin/quantum_analyzer --batch --template single_molecule_overview --files examples/fixtures/methane_6-31g.log
bin/quantum_analyzer --batch --template method_comparison_analysis --files examples/fixtures/methane_6-31g.log examples/fixtures/methane_sto-3g.log
```

Available templates:

| Template | Purpose |
| --- | --- |
| `standard_analysis` | Detailed block report with matrix preview, atomic metrics, bonds and lone-pair tables. |
| `summary_analysis` | Grouped basis-comparison report for multiple calculations of the same molecule. |
| `single_molecule_overview` | Compact one-molecule report with charges, charge deltas, valence and bond indices. |
| `method_comparison_analysis` | Mulliken/Lowdin charge comparison and per-basis method-difference summary. |
| `lone_pair_analysis` | Lone-pair focused report with valence, lone-pair populations and bond indices. |

Use `--out-dir` to save grouped reports as Markdown files:

```bash
bin/quantum_analyzer --batch --template summary_analysis --out-dir examples/generated/analysis_reports --files examples/fixtures/methane_6-31g.log examples/fixtures/methane_sto-3g.log examples/fixtures/water_sto-3g.log
```

## Release packages

Recommended release path: use GitHub Actions. The release workflow builds both
platforms on GitHub runners, tests the packages, and publishes the release
assets:

1. Open `Actions` on GitHub.
2. Run the `Release` workflow.
3. Keep the default tag `v1.0` for the initial release.

The workflow produces and uploads:

- `quantum-analyzer-1.0-linux-x86_64.tar.gz`
- `quantum-analyzer-1.0-windows-x86_64.zip`

For local packaging of the current platform only:

```bash
xmake
xmake run package-release
```

This creates one archive under `dist/` for the current platform:

- Linux: `quantum-analyzer-1.0-linux-x86_64.tar.gz`
- Windows: `quantum-analyzer-1.0-windows-x86_64.zip`

## Repository layout

```text
src/lib/core/       C++ core and public C API
src/lib/lua_core/   LuaJIT FFI wrapper modules
src/app/            Standalone embedded LuaJIT application
src/app/templates/  Lua analysis templates
src/app/share/      External integration resources
tests/              Lua integration tests
examples/           Example Gaussian output files
```

## Notes

This repository is a learning and portfolio project. The mathematical and
domain-specific routines are kept as part of the original implementation; public
cleanup should focus on packaging, documentation, tests and build reliability.

Grouped mode assumes the same report block structure for all files inside one
molecular group and merges values within that shared structure.

Safe mode limits user template Lua surface (`os`/`io`/`debug`/`require` are not
injected), while trusted system plugins may still invoke external tools such as
`gnuplot`.
