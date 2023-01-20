# There is 3 ways to run algo.

# 1. With classes - the slowest but easiest to understand. Dimentions: Nx, Ny, Nz - selects manually, dt selects manually too
   ## 2 ways to run code:
 - ```g++ main.cpp && ./a.out``` 
 - copy in ```CMakeLists.txt``` :
```
cmake_minimum_required(VERSION 3.20)
project(diffusion_mfd_solver)

set(CMAKE_CXX_STANDARD 17)

add_executable(diffusion_mfd_solver main.cpp mfd_solver.h mfd_solver.cpp)
```

# 2. ```solver_linear_algo.cpp``` - is fast linear algo, also Nx,Ny,Nz - selects manually, df - selects automatically
## 2 ways to run code:
- ```g++ main.cpp && ./a.out``` 
- - copy in ```CMakeLists.txt``` :
```
cmake_minimum_required(VERSION 3.20)
project(diffusion_mfd_solver)

set(CMAKE_CXX_STANDARD 17)

add_executable(diffusion_mfd_solver solver_linear_algo.cpp)
```
Also with this program you can build ```*.vtk``` file: you need to uncomment lines:
```
//#include "to_vtk.hpp"   //<----uncomment to build .vtk
//    write_to_file(u_actual,n,m,kk);      //<----uncomment to build .vtk
```
and then run py script: ```python3 to_vtk.py```

results saved to ```numeric.vtk``` and then you can open it for example with [ParaView](https://www.paraview.org/)

# MPI RUN

instructions to parallel algo:
running:
```mpic++ main_mpi.cpp -O2```

```mpi_function.hpp``` - paralleling in one directions - x axis
```mpi_function_xy.hpp```  - paralleling in two directions - x and y