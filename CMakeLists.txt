cmake_minimum_required(VERSION 3.16)

project(femheat
  VERSION 0.1.0
  DESCRIPTION "From-scratch 2D finite element heat-conduction solver"
  LANGUAGES CXX)

# ---------------------------------------------------------------------------
# Global settings
# ---------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default to a Release build (optimized) when the user did not pick one.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Emit compile_commands.json for clang-tidy / editor integration.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Prefer modern Boost config-mode discovery (BoostConfig.cmake) on CMake >= 3.30.
if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()

option(FEMHEAT_BUILD_TESTS "Build the femheat unit/regression tests" ON)

# Warning flags applied to *our* targets only (never to fetched dependencies).
add_library(femheat_warnings INTERFACE)
target_compile_options(femheat_warnings INTERFACE
  $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall;-Wextra>)

# ---------------------------------------------------------------------------
# Third-party dependencies
# ---------------------------------------------------------------------------
include(FetchContent)

# Some third-party FetchContent subprojects (e.g. HighFive 2.9) still declare a
# cmake_minimum_required below 3.5, which CMake >= 4 rejects outright. Provide a
# policy floor so they configure cleanly. Older CMake (the CI image) accepts
# those projects natively and simply ignores this variable, so the build works
# on both toolchains.
set(CMAKE_POLICY_VERSION_MINIMUM 3.5 CACHE STRING "Policy floor for legacy subprojects")

# Keep fetched projects from registering their own ctest tests in our build.
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)

# --- Eigen 3 (header-only sparse linear algebra) ---------------------------
set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_PKGCONFIG OFF CACHE BOOL "" FORCE)
FetchContent_Declare(Eigen3
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG 3.4.0
  GIT_SHALLOW TRUE
  SYSTEM)
FetchContent_MakeAvailable(Eigen3)

# --- HDF5 + HighFive (file I/O) --------------------------------------------
find_package(HDF5 REQUIRED COMPONENTS C)
set(HIGHFIVE_USE_BOOST OFF CACHE BOOL "" FORCE)
set(HIGHFIVE_USE_EIGEN OFF CACHE BOOL "" FORCE)
set(HIGHFIVE_UNIT_TESTS OFF CACHE BOOL "" FORCE)
set(HIGHFIVE_EXAMPLES OFF CACHE BOOL "" FORCE)
set(HIGHFIVE_BUILD_DOCS OFF CACHE BOOL "" FORCE)
FetchContent_Declare(HighFive
  GIT_REPOSITORY https://github.com/BlueBrain/HighFive.git
  GIT_TAG v2.9.0
  GIT_SHALLOW TRUE
  SYSTEM)
FetchContent_MakeAvailable(HighFive)

# --- Boost.Program_options (CLI parsing) -----------------------------------
find_package(Boost REQUIRED COMPONENTS program_options)

# --- GoogleTest (unit / regression tests) ----------------------------------
if(FEMHEAT_BUILD_TESTS)
  set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  FetchContent_Declare(googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.14.0
    GIT_SHALLOW TRUE
    SYSTEM)
  FetchContent_MakeAvailable(googletest)
endif()

# ---------------------------------------------------------------------------
# Library target: femheat_lib
# ---------------------------------------------------------------------------
add_library(femheat_lib
  src/build_info.cpp
  src/line2.cpp
  src/mesh1d.cpp
  src/boundary_condition.cpp
  src/assembler.cpp
  src/linear_solver.cpp
  src/heat_problem.cpp
  src/mms.cpp
  src/tri3.cpp
  src/tri_mesh.cpp)
target_include_directories(femheat_lib PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(femheat_lib PRIVATE
  FEMHEAT_VERSION="${PROJECT_VERSION}")
# HighFive 2.x exposes a plain INTERFACE target named "HighFive" (it carries the
# HDF5 include/link dependencies transitively) when consumed via FetchContent.
target_link_libraries(femheat_lib
  PUBLIC  Eigen3::Eigen HighFive
  PRIVATE femheat_warnings)

# ---------------------------------------------------------------------------
# Executable target: femheat (CLI)
# ---------------------------------------------------------------------------
add_executable(femheat apps/femheat.cpp)
target_link_libraries(femheat
  PRIVATE femheat_lib Boost::program_options femheat_warnings)

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
if(FEMHEAT_BUILD_TESTS)
  enable_testing()
  add_executable(femheat_tests
    tests/test_smoke.cpp
    tests/test_line2.cpp
    tests/test_assembler.cpp
    tests/test_solve1d.cpp
    tests/test_mms1d.cpp
    tests/test_neumann.cpp
    tests/test_tri3.cpp
    tests/test_mms2d.cpp)
  target_link_libraries(femheat_tests
    PRIVATE femheat_lib GTest::gtest_main femheat_warnings)
  include(GoogleTest)
  gtest_discover_tests(femheat_tests)
endif()
