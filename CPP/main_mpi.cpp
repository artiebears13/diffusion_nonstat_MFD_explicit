#include <iostream>
#include <cmath>
#include <chrono>
#include <mpi.h>
#include "mpi_function_xy.hpp"



int main(int argc, char **argv) {
    MPI_Status status;
    int size, id;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    double u_act[n*m*kk];
    double u_n[n*m*kk];
    double *u_actual = &u_act[0];
    double *u_next = &u_n[0];
    fill_boundary(u_actual,u_next);


    solver(u_actual,u_next,id,size);
//    get_mist(u_actual,id);

    MPI_Finalize();
    return 0;
}

