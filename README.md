# femheat

*A from-scratch 2D finite-element solver for steady-state heat conduction, written
in modern C++17 and validated by the Method of Manufactured Solutions.*

[![CI](https://github.com/thefcan/femheat/actions/workflows/ci.yml/badge.svg)](https://github.com/thefcan/femheat/actions/workflows/ci.yml)
[![C++17](https://img.shields.io/badge/C%2B%2B-17-00599C?logo=cplusplus&logoColor=white)](https://en.cppreference.com/w/cpp/17)
[![CMake](https://img.shields.io/badge/CMake-%E2%89%A5%203.16-064F8C?logo=cmake&logoColor=white)](https://cmake.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

<p align="center">
  <img src="docs/temperature_field.png" alt="2D temperature field T = sin(pi x) sin(pi y)" width="420"/>
  <br/>
  <em>Temperature field from <code>femheat --case mms2d</code> on the unit square
  (real solver output, rendered from the exported node values).</em>
</p>

## Highlights

- **From-scratch FEM** — Galerkin weak form, element matrices and sparse assembly
  written by hand, with no external FEM library.
- **1D and 2D** — linear segment (`Line2`) and triangle (`Tri3`) elements behind a
  single polymorphic interface, so one assembly/solve path serves both.
- **Validated, not just running** — an MMS convergence study confirms the
  theoretical **second-order** accuracy: observed L2 order **2.00 (1D)** and
  **1.996 (2D)**.
- **Sparse SPD solve** — Jacobi-preconditioned conjugate gradient (Eigen).
- **Real I/O** — HDF5 + XDMF output for ParaView and a Boost.Program_options CLI.
- **Engineered** — 14 GoogleTest cases, a warning-free build (`-Wall -Wextra`), and
  GitHub Actions CI on every push.

## Quick start

```bash
git clone https://github.com/thefcan/femheat.git
cd femheat
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build

# solve the 2D manufactured problem and write ParaView output
./build/femheat --case mms2d --nx 41 --ny 41 --out result.h5
```

Eigen, HighFive and GoogleTest are fetched automatically; you only need a C++17
compiler, CMake, HDF5 and Boost on the system. See **Build & test** below for the
per-platform package list.

## What it solves

femheat solves the steady-state heat equation

$$-\nabla \cdot (k\,\nabla T) = f \quad \text{in } \Omega,$$

with Dirichlet ($T = g$) and Neumann (prescribed flux) boundary conditions, in
both 1D and 2D.

## How it works

By the Galerkin finite element method: multiplying the PDE by a test function $v$
and integrating by parts gives the **weak form**

$$\int_\Omega k\,\nabla T \cdot \nabla v \, d\Omega \;=\; \int_\Omega f\,v \, d\Omega \;+\; \text{(boundary flux terms)}.$$

The temperature is approximated by **piecewise-linear shape functions** on a mesh
of elements — 2-node segments (`Line2`) in 1D, 3-node triangles (`Tri3`) in 2D —
which turns the weak form into a sparse **symmetric positive-definite** system
$K\,T = f$. Dirichlet conditions are imposed with a penalty term; Neumann
conditions enter naturally through the boundary integral. The system is solved
with a Jacobi-preconditioned **conjugate gradient** iteration (Eigen).

Because the elements are linear, the discretization is **second-order accurate**:
halving the mesh size cuts the L2 error by ≈ 4×. We confirm this with the Method
of Manufactured Solutions (MMS) — pick an exact $u$, derive the source $f$ it
implies, solve, and measure the observed convergence order (see the convergence
study below).

## Architecture

```
                          ┌──────────────┐
        mesh + material   │  HeatProblem │   orchestration  (solve())
        + source + BCs    └──────┬───────┘
                                 │
   Mesh1D / TriMesh  ──build──►  IElement  ◄── Material (k)
     (geometry)                 Line2 / Tri3
                                 │  elementStiffness() / elementLoad()
                                 ▼
                           ┌───────────┐        IBoundaryCondition
                           │ Assembler │ ◄────  DirichletBC (penalty)
                           └─────┬─────┘        NeumannBC   (flux)
                    sparse K, f  │
                                 ▼
                          ┌──────────────┐
                          │ LinearSolver │   Eigen ConjugateGradient (SPD)
                          └──────┬───────┘
                       nodal T   │
                                 ▼
                          ┌──────────────┐
                          │  Hdf5Writer  │ → result.h5 + result.xdmf (ParaView)
                          └──────────────┘
```

- `IElement` is an abstract element interface (runtime **polymorphism**), so the
  same `Assembler` serves both 1D and 2D.
- `IBoundaryCondition` is a **strategy** interface (`DirichletBC`, `NeumannBC`).

## Tech stack

| Concern         | Tool                                                                  |
| --------------- | --------------------------------------------------------------------- |
| Language / build | C++17, CMake (≥ 3.16)                                                 |
| Linear algebra  | [Eigen 3](https://eigen.tuxfamily.org) — sparse matrices + CG          |
| File I/O        | [HDF5](https://www.hdfgroup.org/) via [HighFive](https://github.com/BlueBrain/HighFive), with a companion XDMF |
| CLI             | [Boost.Program_options](https://www.boost.org/)                       |
| Testing         | GoogleTest + CTest                                                     |
| CI              | GitHub Actions (ubuntu-latest)                                         |

Eigen, HighFive and GoogleTest are fetched automatically via CMake
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

## Running the solver

```bash
# 2D manufactured solution on the unit square -> HDF5 + XDMF for ParaView
./build/femheat --case mms2d --nx 40 --ny 40 --k 1.0 --out result.h5

# 1D manufactured solution, reports its L2 error
./build/femheat --case mms1d --nx 80 --k 1.0

# options
./build/femheat --help
```

Open `result.xdmf` in [ParaView](https://www.paraview.org/) to visualize the
temperature field. The optional `--csv field.csv` flag also dumps the node
values, which `scripts/plot_field.py` renders to a PNG with **no third-party
Python packages** required:

```bash
./build/femheat --case mms2d --nx 81 --ny 81 --csv field.csv --out result.h5
python3 scripts/plot_field.py field.csv 81 81 docs/temperature_field.png
```

## Using the library

The CLI is a thin front-end over `femheat_lib`. The library API mirrors the
math — build a mesh, define the problem, add boundary conditions, solve:

```cpp
#include <iostream>

#include "femheat/hdf5_writer.hpp"
#include "femheat/heat_problem.hpp"
#include "femheat/material.hpp"
#include "femheat/tri_mesh.hpp"

using namespace femheat;

int main() {
  // 1. Triangulate the unit square with a 41 x 41 node grid.
  TriMesh mesh(/*Lx=*/1.0, /*Ly=*/1.0, /*nx=*/41, /*ny=*/41);

  // 2. Conductivity k = 1, with a constant volumetric heat source f = 10.
  HeatProblem problem(mesh, Material(/*k=*/1.0));
  problem.setConstantSource(10.0);

  // 3. Hold T = 0 on the whole boundary (penalty Dirichlet).
  for (int node : mesh.boundaryNodes()) {
    problem.addDirichlet(node, 0.0);
  }

  // 4. Assemble and solve the sparse SPD system  K T = f  with CG.
  const auto result = problem.solveDetailed();
  std::cout << "CG converged in " << result.iterations << " iterations\n";

  // 5. Write result.h5 + result.xdmf for ParaView.
  Hdf5Writer::writeTriMesh("result.h5", mesh, result.x);
}
```

## Validation — convergence study

Both solvers are validated with MMS. The observed order is

$$p = \frac{\log(e_{\text{coarse}} / e_{\text{fine}})}{\log(h_{\text{coarse}} / h_{\text{fine}})},$$

and the regression tests assert $p \in [1.7, 2.3]$. The numbers below are the
actual output of the test suite.

**1D**, $u(x) = \sin(\pi x)$:

| elements `n` | L2 error    | observed order |
| -----------: | ----------- | -------------: |
|           10 | 9.0855e-03  |        —       |
|           20 | 2.2747e-03  |      1.998     |
|           40 | 5.6887e-04  |      1.999     |
|           80 | 1.4223e-04  |      2.000     |

**2D**, $u(x, y) = \sin(\pi x)\,\sin(\pi y)$:

| nodes / side | L2 error    | observed order |
| -----------: | ----------- | -------------: |
|           11 | 1.5778e-02  |        —       |
|           21 | 3.9913e-03  |      1.983     |
|           41 | 1.0008e-03  |      1.996     |

The error drops by ≈ 4× per refinement and the order approaches 2, exactly as
theory predicts for linear elements. ✅

## Project layout

```
femheat/
├─ include/femheat/   public headers (Mesh, IElement, Assembler, ... )
├─ src/               library implementation  -> femheat_lib
├─ apps/femheat.cpp   CLI executable          -> femheat
├─ tests/             GoogleTest unit + regression tests -> femheat_tests
├─ scripts/           stdlib-only field plotter
└─ docs/              figures
```

## Scope & limitations

femheat is a focused, correctness-first solver with a deliberately small and
honestly bounded scope:

- **Steady state only** — no time-dependent term.
- **Linear elements** — `Line2` in 1D, `Tri3` in 2D (second-order accurate).
- **Structured meshes** — uniform 1D intervals and structured triangulations of a
  rectangle; no unstructured-mesh import yet.
- **Isotropic, single material** — `Material` stores one global conductivity `k`.
- **Serial** — no threading or MPI.

The architecture was built to grow past these. Because every element lives behind
the `IElement` interface and the assembly/solve path is element-agnostic, the
natural next steps — **3D tetrahedra**, a **transient term**, **per-element
materials**, or **higher-order elements** — slot in without touching the core
pipeline.

## References

- T. J. R. Hughes, *The Finite Element Method: Linear Static and Dynamic Finite
  Element Analysis*, Dover, 2000.
- P. J. Roache, "Code Verification by the Method of Manufactured Solutions,"
  *Journal of Fluids Engineering*, 124(1), 4–10, 2002.
- G. Guennebaud, B. Jacob, et al., *Eigen v3*, <https://eigen.tuxfamily.org>.

## License

[MIT](LICENSE) © 2026 Furkan Karafil
