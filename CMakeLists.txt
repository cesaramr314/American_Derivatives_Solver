cmake_minimum_required(VERSION 3.16)
project(american_derivatives_solver LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Tipo de build" FORCE)
endif()

# Warnings. Each compilator uses its own syntax.
if(MSVC)
  add_compile_options(/W4 /utf-8)   # Visual Studio
else()
  add_compile_options(-Wall -Wextra -Wpedantic)   # g++ y clang
endif()

# ---- Libraries ----
add_library(amopt
  src/tridiag.cpp
  src/lcp.cpp
  src/pde_solver.cpp
  src/analytic.cpp
  src/binomial.cpp
)
target_include_directories(amopt PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# ---- Executables ----
add_executable(amopt_price apps/price.cpp)
target_link_libraries(amopt_price PRIVATE amopt)

add_executable(amopt_benchmark apps/benchmark.cpp)
target_link_libraries(amopt_benchmark PRIVATE amopt)

# ---- Tests ----
enable_testing()
add_executable(amopt_tests tests/test_solver.cpp)
target_link_libraries(amopt_tests PRIVATE amopt)
add_test(NAME amopt_tests COMMAND amopt_tests)