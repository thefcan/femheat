# femheat

[![CI](https://github.com/thefcan/femheat/actions/workflows/ci.yml/badge.svg)](https://github.com/thefcan/femheat/actions/workflows/ci.yml)

A from-scratch **2D finite element method (FEM) solver** for steady-state heat
conduction, written in modern C++17.

> **Status:** early development — building milestone by milestone. See the
> roadmap below.

## Overview

femheat solves the steady-state heat equation

```
-∇·(k ∇T) = f      on a domain Ω,
```

with Dirichlet and Neumann boundary conditions, using linear finite elements
(1D `Line2`, 2D linear triangles `Tri3`). The implementation is validated with
the **Method of Manufactured Solutions (MMS)**: we verify that the observed
spatial convergence order approaches 2, as theory predicts for linear elements.
Solutions are exported to HDF5 + XDMF for visualization in ParaView.

## Tech stack

| Concern            | Tool                                                   |
| ------------------ | ------------------------------------------------------ |
| Language / build   | C++17, CMake (≥ 3.16)                                   |
| Linear algebra     | [Eigen 3](https://eigen.tuxfamily.org) — sparse + CG   |
| File I/O           | [HDF5](https://www.hdfgroup.org/) via [HighFive](https://github.com/BlueBrain/HighFive) |
| CLI                | [Boost.Program_options](https://www.boost.org/)        |
| Testing            | GoogleTest + CTest                                     |
| CI                 | GitHub Actions (ubuntu-latest)                         |

Eigen, HighFive and GoogleTest are pulled in automatically via CMake
`FetchContent`; HDF5 and Boost are located with `find_package`.

## Build & test

### Linux (Ubuntu)

```bash
sudo apt-get update
sudo apt-get install -y build-essential cmake libhdf5-dev libboost-program-options-dev
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

### macOS

```bash
brew install cmake hdf5 boost
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

A fresh clone builds with just `cmake -B build && cmake --build build`.

## Roadmap

- [x] **M0** — Skeleton & toolchain (CMake, dependencies, CI, smoke tests)
- [ ] **M1** — 1D FEM core (`Line2`, assembly, Dirichlet BC, Eigen solve)
- [ ] **M2** — 1D MMS validation + convergence study; Neumann BC
- [ ] **M3** — 2D linear triangles (`Tri3`) + 2D MMS validation
- [ ] **M4** — HDF5/XDMF output + Boost.Program_options CLI
- [ ] **M5** — Docs, real convergence table, ParaView screenshot, `v0.1` tag

## License

[MIT](LICENSE) © 2026 Furkan Karafil
