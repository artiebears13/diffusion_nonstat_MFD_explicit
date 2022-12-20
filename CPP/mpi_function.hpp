#include <iostream>
#include <cmath>
#include <chrono>
//#include <mpi.h>


#define get_analytic_u(x, y, z) (sin(M_PI * x) * sin(M_PI * y) * sin(M_PI * z) *(1 - exp(-((D0 + D1 + D2) * M_PI*M_PI * 1))))
//#define f(x, y, z) ((D0+D1+D2)*(M_PI*M_PI*sin(M_PI*x))* sin(M_PI * y) * sin(M_PI * z))

const unsigned Nx = 14;
const unsigned Ny = 10;
const unsigned Nz = 10;
const unsigned npx = 4;
const unsigned npy = 1;
const unsigned n = (Nx - 2) / npx + 2;
//const unsigned m = (Ny - 2) / npy + 2;
const unsigned m = Ny;
const unsigned kk = Nz;
const unsigned len = n * m * kk;
const double hx = 1. / (Nx - 1) * 1.0;
const double hy = 1. / (Ny - 1) * 1.0;
const double hz = 1. / (Nz - 1) * 1.0;
const double D0 = 0.25;
const double D1 = 0.15;
const double D2 = 0.1;
const double dt = 0.9 * (0.5 / (D0 / (hx * hx) + D1 / (hy * hy) + D2 / (hz * hz)));


inline double f(double x, double y, double z) {
    return ((D0 + D1 + D2) * (M_PI * M_PI * sin(M_PI * x)) * sin(M_PI * y) * sin(M_PI * z));
}


inline double get_u(double *const &u, int i, int j, int k) {
    return *(u + k * n * m + j * n + i);
}

inline void set_u(double *const &u, int i, int j, int k, double const &value) {
    *(u + k * n * m + j * n + i) = value;
}


inline double Lx(double *const &u_actual, int i, int j, int k) {
    return (get_u(u_actual, i - 1, j, k) - 2 * get_u(u_actual, i, j, k) + get_u(u_actual, i + 1, j, k))
           / (hx * hx);
}

inline double Ly(double *const &u_actual, int i, int j, int k) {
    return (get_u(u_actual, i, j - 1, k) - 2 * get_u(u_actual, i, j, k) + get_u(u_actual, i, j + 1, k))
           / (hy * hy);
}

inline double Lz(double *const &u_actual, int i, int j, int k) {
    return (get_u(u_actual, i, j, k - 1) - 2 * get_u(u_actual, i, j, k) + get_u(u_actual, i, j, k + 1))
           / ((hz) * (hz));
}

inline void swap(double *&u_actual, double *&u_next) {

    double *temp = u_actual;
    u_actual = u_next;
    u_next = temp;
}


inline void get_mist(double *&u_actual, int id) {
    double mist = 0.;
    for (int k = 0; k < kk; ++k) {
        for (int j = 0; j < m; ++j) {
            for (int i = 1; i < n - 1; ++i) {
                if (id==1){
                    std::cout << "id: " << id << " global i: " << (id * (n - 2) + i) <<" local i: "<< i <<" j: "<<j*hy<<" k: "<<k*hz<<" VALUE: " << get_u(u_actual, i, j, k)
                              << " and " << get_analytic_u((id * (n - 2) + i) * hx, j * hy, k * hz) << std::endl;
                }
                if (std::fabs(get_u(u_actual, i, j, k) - get_analytic_u((id * (n-2) + i) * hx, j * hy, k * hz)) >
                    mist) {
                    mist = std::fabs(
                            get_u(u_actual, i, j, k) - get_analytic_u((id * (n-2) + i) * hx, j * hy, k * hz));


                }
            }
        }
    }

    std::cout << "id: " << id << " MIST: " << mist << std::endl;
    std::cout << "dt: " << dt << std::endl;
}


inline void fill_boundary(double *&u_actual, double *&u_next) {
    for (int k = 0; k < kk; ++k) {
        for (int j = 0; j < m; ++j) {
            for (int i = 0; i < n; ++i) {
                set_u(u_actual, i, j, k, 0.);
                set_u(u_next, i, j, k, 0.);
            }
        }
    }
}


inline void pack_left(double *&buffer, double *&u_actual) {  //i=1
    for (int k = 0; k < kk; ++k) {
        for (int j = 0; j < m; ++j) {
            *(buffer + k * Ny + j) = get_u(u_actual, 1, j, k);
            if (j == 5 and k == 5) {
//                std::cout << "BUFFER " << *(buffer + k * Ny + j) << "GET_U" << get_u(u_actual, 1, j, k) << std::endl;
            }
        }
    }
}


inline void pack_right(double *&buffer, double *&u_actual) {  //i=n-2
    for (int k = 0; k < kk; ++k) {
        for (int j = 0; j < m; ++j) {
            *(buffer + k * Ny + j) = get_u(u_actual, n - 2, j, k);
        }
    }
}


inline void unpack_left(double *&buffer, double *&u_actual) {  //i=n-2
    for (int k = 0; k < kk; ++k) {
        for (int j = 0; j < m; ++j) {
            set_u(u_actual, 0, j, k, *(buffer + k * Ny + j));
        }
    }
}

inline void unpack_right(double *&buffer, double *&u_actual) {  //i=n-2
    for (int k = 0; k < kk; ++k) {
        for (int j = 0; j < m; ++j) {
            set_u(u_actual, n - 1, j, k, *(buffer + k * Ny + j));
        }
    }
}

inline void solver(double *&u_actual, double *&u_next, int id, int size) {
    double u;
    MPI_Request *request;
    MPI_Status status;
    double buf_send_left[m * kk];
    double *buffer_send_left = &buf_send_left[0];

    double buf_send_right[m * kk];
    double *buffer_send_right = &buf_send_right[0];

    double buf_rec_left[m * kk];
    double *buffer_rec_left = &buf_rec_left[0];

    double buf_rec_right[m * kk];
    double *buffer_rec_right = &buf_rec_right[0];
    MPI_Barrier(MPI_COMM_WORLD);
    int counter = 0;
    for (int l = 1; l * dt <= 1; l++) {

//        std::cout << "==================================================================" << std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        if (id == 0) {

            pack_right(buffer_send_right, u_actual);
//            std::cout << "id "<<id<<" ===BUFFER " << *(buffer_send_right + 5 * Ny + 5) << " GET_U " << get_u(u_actual, n - 2, 5, 5) <<" TO "<<id+1<< std::endl;
//            MPI_Isend(buffer_send_right,m*kk,MPI_DOUBLE,id+1,0,MPI_COMM_WORLD,request);
            MPI_Send(buffer_send_right, m * kk, MPI_DOUBLE, id + 1, 0, MPI_COMM_WORLD);
//
//            MPI_Irecv(buffer_rec_right,m*kk,MPI_DOUBLE,id+1,0,MPI_COMM_WORLD,request);
            MPI_Recv(buffer_rec_right, m * kk, MPI_DOUBLE, id + 1, 0, MPI_COMM_WORLD, &status);
            unpack_right(buffer_rec_right, u_actual);
//            std::cout << "id "<<id<<" ===BUFFER " << *(buffer_rec_right + 5 * Ny + 5) << " GET_U " << get_u(u_actual, n-1, 5, 5) <<" FROM "<<id+1<< std::endl;


        } else if (id == size - 1) {
            if (id % 2 == 0) {
                pack_left(buffer_send_left, u_actual);
//                std::cout << "id "<<id<<" ===BUFFER " << *(buffer_send_left + 5 * Ny + 5) << " GET_U " << get_u(u_actual, 1, 5, 5) <<" TO "<<id-1<< std::endl;


//            MPI_Isend(buffer_send_left,m*kk,MPI_DOUBLE,id-1,0,MPI_COMM_WORLD,request);
                MPI_Send(buffer_send_left, m * kk, MPI_DOUBLE, id - 1, 0, MPI_COMM_WORLD);

//            MPI_Irecv(buffer_rec_left,m*kk,MPI_DOUBLE,id-1,0,MPI_COMM_WORLD,request);
                MPI_Recv(buffer_rec_left, m * kk, MPI_DOUBLE, id - 1, 0, MPI_COMM_WORLD, &status);

                unpack_left(buffer_rec_left, u_actual);
//                std::cout << "id "<<id<<" ===BUFFER " << *(buffer_rec_left + 5 * Ny + 5) << " GET_U " << get_u(u_actual, 0, 5, 5) <<" FROM "<<id-1<< std::endl;
            } else {
                MPI_Recv(buffer_rec_left, m * kk, MPI_DOUBLE, id - 1, 0, MPI_COMM_WORLD, &status);
                unpack_left(buffer_rec_left, u_actual);
//                std::cout << "id "<<id<<" ===BUFFER " << *(buffer_rec_left + 5 * Ny + 5) << " GET_U " << get_u(u_actual, 0, 5, 5) <<" FROM "<<id-1<< std::endl;
                pack_left(buffer_send_left, u_actual);
//                std::cout << "id "<<id<<" ===BUFFER " << *(buffer_send_left + 5 * Ny + 5) << " GET_U " << get_u(u_actual, 1, 5, 5) <<" TO "<<id-1<< std::endl;

//            MPI_Isend(buffer_send_left,m*kk,MPI_DOUBLE,id-1,0,MPI_COMM_WORLD,request);
                MPI_Send(buffer_send_left, m * kk, MPI_DOUBLE, id - 1, 0, MPI_COMM_WORLD);
            }

        } else {
            if (id % 2 == 0) {
                pack_left(buffer_send_left, u_actual);
//                std::cout << "id "<<id<<" ===BUFFER " << *(buffer_send_left + 5 * Ny + 5) << " GET_U " << get_u(u_actual, 1, 5, 5) <<" TO "<<id-1<< std::endl;
                MPI_Send(buffer_send_left, m * kk, MPI_DOUBLE, id - 1, 0, MPI_COMM_WORLD);

                pack_right(buffer_send_right, u_actual);
//                std::cout << "id "<<id<<" ===BUFFER " << *(buffer_send_right + 5 * Ny + 5) << " GET_U " << get_u(u_actual, n - 2, 5, 5) <<" TO "<<id+1<< std::endl;
                MPI_Send(buffer_send_right, m * kk, MPI_DOUBLE, id + 1, 0, MPI_COMM_WORLD);

                MPI_Recv(buffer_rec_left, m * kk, MPI_DOUBLE, id - 1, 0, MPI_COMM_WORLD, &status);
                unpack_left(buffer_rec_left, u_actual);
//                std::cout << "id "<<id<<" ===BUFFER " << *(buffer_rec_left + 5 * Ny + 5) << " GET_U " << get_u(u_actual, 0, 5, 5) <<" FROM "<<id-1<< std::endl;

                MPI_Recv(buffer_rec_right, m * kk, MPI_DOUBLE, id + 1, 0, MPI_COMM_WORLD, &status);
                unpack_right(buffer_rec_right, u_actual);
//                std::cout << "id "<<id<<" ===BUFFER " << *(buffer_rec_right + 5 * Ny + 5) << " GET_U " << get_u(u_actual, n-1, 5, 5) <<" FROM "<<id+1<< std::endl;
            } else {
                MPI_Recv(buffer_rec_left, m * kk, MPI_DOUBLE, id - 1, 0, MPI_COMM_WORLD, &status);
                unpack_left(buffer_rec_left, u_actual);
//                std::cout << "id "<<id<<" ===BUFFER " << *(buffer_rec_left + 5 * Ny + 5) << " GET_U " << get_u(u_actual, 0, 5, 5) <<" FROM "<<id-1<< std::endl;

                MPI_Recv(buffer_rec_right, m * kk, MPI_DOUBLE, id + 1, 0, MPI_COMM_WORLD, &status);
                unpack_right(buffer_rec_right, u_actual);
//                std::cout << "id "<<id<<" ===BUFFER " << *(buffer_rec_right + 5 * Ny + 5) << " GET_U " << get_u(u_actual, n-1, 5, 5) <<" FROM "<<id+1<< std::endl;

                pack_left(buffer_send_left, u_actual);
//                std::cout << "id "<<id<<" ===BUFFER " << *(buffer_send_left + 5 * Ny + 5) << " GET_U " << get_u(u_actual, 1, 5, 5) <<" TO "<<id-1<< std::endl;
                MPI_Send(buffer_send_left, m * kk, MPI_DOUBLE, id - 1, 0, MPI_COMM_WORLD);

                pack_right(buffer_send_right, u_actual);
//                std::cout << "id "<<id<<" ===BUFFER " << *(buffer_send_right + 5 * Ny + 5) << " GET_U " << get_u(u_actual, n - 2, 5, 5) <<" TO "<<id+1<< std::endl;
                MPI_Send(buffer_send_right, m * kk, MPI_DOUBLE, id + 1, 0, MPI_COMM_WORLD);
            }

        }
        MPI_Barrier(MPI_COMM_WORLD);

        for (unsigned k = 1; k < kk - 1; ++k) {
                for (unsigned j = 1; j < m - 1; ++j) {
                for (unsigned i = 1; i < n-1; ++i) {
                    u = get_u(u_actual, i, j, k) +
                        dt * (f((id * (n-2) + i) * hx, j * hy, k * hz)
                              + D0 * Lx(u_actual, i, j, k)
                              + D1 * Ly(u_actual, i, j, k)
                              + D2 * Lz(u_actual, i, j, k));
//                    std::cout << "dt: " << dt << " get_u: " << get_u(u_actual, i, j, k) << " " << " f(i, j, k): "
//                              << f(i * hx, j * hy, k * hz) <<
//                              " D0 * Lx(i, j, k): " << D0 * Lx(u_actual, i, j, k) << " D1 * Ly(i, j, k): "
//                              << D1 * Ly(u_actual, i, j, k)
//                              << " D2 * Lz(i, j, k): " << D2 * Lz(u_actual, i, j, k) << " u: "
//                              << u << std::endl;

//                        std::cout<<"i: "<< k * n * m + j * n + i <<"   "<< *(u_actual+k * n * m + j * n + i)<<"   " << get_u(u_actual, i,j,k)<<std::endl;

                    set_u(u_next, i, j, k, u);
//                    }

                }


            }

        }
        swap(u_actual, u_next);
        MPI_Barrier(MPI_COMM_WORLD);
        counter++;
    }
    std::cout << "id " << id << " counter: " << counter << std::endl;

}
