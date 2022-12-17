cmake_minimum_required(VERSION 3.20)
project(diffusion_mfd_solver)

set(CMAKE_CXX_STANDARD 17)

add_executable(diffusion_mfd_solver main3.cpp to_vtk.hpp)