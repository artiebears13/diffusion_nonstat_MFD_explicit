cmake_minimum_required(VERSION 3.20)
project(diffusion_mfd_solver)

set(CMAKE_CXX_STANDARD 17)

add_executable(diffusion_mfd_solver solver_linear_algo.cpp mpi_function_xy.hpp linear_function.hpp mpi_function.hpp main_mpi.cpp)