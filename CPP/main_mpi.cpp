#include <iostream>
#include <cmath>
#include <chrono>
#include <mpi.h>
#include "linear_function.hpp"



int main(int argc, char **argv) {
    MPI_Status status;
    int size, id;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    const unsigned n =(Nx-2)/size+2;
    const unsigned m =Ny;
    const unsigned kk= Nz;
    double u_act[len];
    double u_n[len];
    double *u_actual = &u_act[0];
    double *u_next = &u_n[0];
    fill_boundary(u_actual,u_next);

    clock_t tStart = clock();
    solver(u_actual,u_next);
    double time = (double) (clock() - tStart) / CLOCKS_PER_SEC;

    std::cout << "time: " << time << std::endl;
    get_mist(u_actual);

    return 0;
}

