#include <iostream>
#include <cmath>
#include <chrono>
//#include <mpi.h>


#define get_analytic_u(x, y, z) (sin(M_PI * x) * sin(M_PI * y) * sin(M_PI * z) *(1 - exp(-((D0 + D1 + D2) * M_PI*M_PI * 1))))


const unsigned Nx = 122;
const unsigned Ny = 62;
const unsigned Nz = 22;
const unsigned npx = 2;
const unsigned npy = 2;
const unsigned n = (Nx - 2) / npx + 2;
const unsigned m = (Ny - 2) / npy + 2;
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


inline void get_mist(double *&u_actual, int id,int row, int col) {
    double mist = 0.;
    for (int k = 0; k < kk; ++k) {
        for (int j = 0; j < m; ++j) {
            for (int i = 1; i < n - 1; ++i) {

                if (std::fabs(get_u(u_actual, i, j, k) - get_analytic_u((col * (n-2) + i) * hx, (row * (m-2) + j) * hy, k * hz)) >
                    mist) {
                    mist = std::fabs(
                            get_u(u_actual, i, j, k) - get_analytic_u((col * (n-2) + i) * hx, (row * (m-2) + j) * hy, k * hz));
                    //                 if (id==3){
                    //        std::cout << "id: " << id <<" row: " <<row <<" col:" << col <<" global i: " << (col * (n - 2) + i) << " local i: " << i
                    //                  << " global j: " << (row * (m - 2) + j) << " local j: " << j  << " k: " << k
                    //                  << " VALUE: " << get_u(u_actual, i, j, k)
                    //                  << " and " << get_analytic_u((col * (n - 2) + i) * hx, (row * (m - 2) + j) * hy, k * hz)
                    //                  << std::endl;
                    //    }
                }
            }
        }
    }
    double mistake;
    MPI_Reduce(&mist,&mistake,1,MPI_DOUBLE, MPI_MAX,0,MPI_COMM_WORLD);

    if (id==0){
        std::cout << "id: " << id << " MIST: " << mistake << std::endl;
        std::cout << "dt: " << dt << std::endl;
    }
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

inline void pack_up(double *&buffer, double *&u_actual) {  //i=1
    for (int k = 0; k < kk; ++k) {
        for (int i = 0; i < n; ++i) {
            *(buffer + k * n + i) = get_u(u_actual, i, m-2, k);
        }
    }
}

inline void pack_down(double *&buffer, double *&u_actual) {  //i=1
    for (int k = 0; k < kk; ++k) {
        for (int i = 0; i < n; ++i) {
            *(buffer + k * n + i) = get_u(u_actual, i, 1, k);
        }
    }
}

inline void pack_left(double *&buffer, double *&u_actual) {  //i=1
    for (int k = 0; k < kk; ++k) {
        for (int j = 0; j < m; ++j) {
            *(buffer + k * m + j) = get_u(u_actual, 1, j, k);
        }
    }
}


inline void pack_right(double *&buffer, double *&u_actual) {  //i=n-2
    for (int k = 0; k < kk; ++k) {
        for (int j = 0; j < m; ++j) {
            *(buffer + k * m + j) = get_u(u_actual, n - 2, j, k);
        }
    }
}

inline void unpack_up(double *&buffer, double *&u_actual) {  //i=n-2
    for (int k = 0; k < kk; ++k) {
        for (int i = 0; i < n; ++i) {
            set_u(u_actual, i, m-1, k, *(buffer + k * n + i));
        }
    }
}

inline void unpack_down(double *&buffer, double *&u_actual) {  //i=n-2
    for (int k = 0; k < kk; ++k) {
        for (int i = 0; i < n; ++i) {
            set_u(u_actual, i, 0, k, *(buffer + k * n + i));
        }
    }
}

inline void unpack_left(double *&buffer, double *&u_actual) {  //i=n-2
    for (int k = 0; k < kk; ++k) {
        for (int j = 0; j < m; ++j) {
            set_u(u_actual, 0, j, k, *(buffer + k * m + j));
        }
    }
}

inline void unpack_right(double *&buffer, double *&u_actual) {  //i=n-2
    for (int k = 0; k < kk; ++k) {
        for (int j = 0; j < m; ++j) {
            set_u(u_actual, n - 1, j, k, *(buffer + k * m + j));
        }
    }
}

inline void chatting_horizontal(int &id, int &row, int &col,double* &u_actual, double* &buffer_send_right,double* &buffer_rec_right,
                                double* &buffer_send_left,double* &buffer_rec_left, MPI_Status status){
    if (col == 0) {

        pack_right(buffer_send_right, u_actual);
        MPI_Send(buffer_send_right, m * kk, MPI_DOUBLE, id + 1, 0, MPI_COMM_WORLD);

        MPI_Recv(buffer_rec_right, m * kk, MPI_DOUBLE, id + 1, 0, MPI_COMM_WORLD, &status);
        unpack_right(buffer_rec_right, u_actual);

    } else if (col == npx - 1) {
        if (col % 2 == 0) {
            pack_left(buffer_send_left, u_actual);
            MPI_Send(buffer_send_left, m * kk, MPI_DOUBLE, id - 1, 0, MPI_COMM_WORLD);

            MPI_Recv(buffer_rec_left, m * kk, MPI_DOUBLE, id - 1, 0, MPI_COMM_WORLD, &status);
            unpack_left(buffer_rec_left, u_actual);
        } else {
            MPI_Recv(buffer_rec_left, m * kk, MPI_DOUBLE, id - 1, 0, MPI_COMM_WORLD, &status);
            unpack_left(buffer_rec_left, u_actual);

            pack_left(buffer_send_left, u_actual);
            MPI_Send(buffer_send_left, m * kk, MPI_DOUBLE, id - 1, 0, MPI_COMM_WORLD);
        }

    } else {
        if (col % 2 == 0) {
            pack_left(buffer_send_left, u_actual);
            MPI_Send(buffer_send_left, m * kk, MPI_DOUBLE, id - 1, 0, MPI_COMM_WORLD);

            pack_right(buffer_send_right, u_actual);
            MPI_Send(buffer_send_right, m * kk, MPI_DOUBLE, id + 1, 0, MPI_COMM_WORLD);

            MPI_Recv(buffer_rec_left, m * kk, MPI_DOUBLE, id - 1, 0, MPI_COMM_WORLD, &status);
            unpack_left(buffer_rec_left, u_actual);

            MPI_Recv(buffer_rec_right, m * kk, MPI_DOUBLE, id + 1, 0, MPI_COMM_WORLD, &status);
            unpack_right(buffer_rec_right, u_actual);
        } else {
            MPI_Recv(buffer_rec_left, m * kk, MPI_DOUBLE, id - 1, 0, MPI_COMM_WORLD, &status);
            unpack_left(buffer_rec_left, u_actual);

            MPI_Recv(buffer_rec_right, m * kk, MPI_DOUBLE, id + 1, 0, MPI_COMM_WORLD, &status);
            unpack_right(buffer_rec_right, u_actual);

            pack_left(buffer_send_left, u_actual);
            MPI_Send(buffer_send_left, m * kk, MPI_DOUBLE, id - 1, 0, MPI_COMM_WORLD);

            pack_right(buffer_send_right, u_actual);
            MPI_Send(buffer_send_right, m * kk, MPI_DOUBLE, id + 1, 0, MPI_COMM_WORLD);
        }

    }


}

inline void chatting_vertical(int &id, int &row, int size,double* &u_actual, double* &buffer_send_up,double* &buffer_rec_up,
                              double* &buffer_send_down,double* &buffer_rec_down, MPI_Status status){

    if (row == 0) {

        pack_up(buffer_send_up, u_actual);
        MPI_Send(buffer_send_up, n * kk, MPI_DOUBLE, id + npx, 0, MPI_COMM_WORLD);

        MPI_Recv(buffer_rec_up, n * kk, MPI_DOUBLE, id + npx, 0, MPI_COMM_WORLD, &status);
        unpack_up(buffer_rec_up, u_actual);

    } else if (row == npy - 1) {
        if (row % 2 == 0) {
            pack_down(buffer_send_down, u_actual);
            MPI_Send(buffer_send_down, n * kk, MPI_DOUBLE, id - npx, 0, MPI_COMM_WORLD);

            MPI_Recv(buffer_rec_down, n * kk, MPI_DOUBLE, id - npx, 0, MPI_COMM_WORLD, &status);
            unpack_down(buffer_rec_down, u_actual);
        } else {
            MPI_Recv(buffer_rec_down, n * kk, MPI_DOUBLE, id - npx, 0, MPI_COMM_WORLD, &status);
            unpack_down(buffer_rec_down, u_actual);

            pack_down(buffer_send_down, u_actual);
            MPI_Send(buffer_send_down, n * kk, MPI_DOUBLE, id - npx, 0, MPI_COMM_WORLD);
        }

    } else {
        if (row % 2 == 0) {
            pack_down(buffer_send_down, u_actual);
            MPI_Send(buffer_send_down, n * kk, MPI_DOUBLE, id - npx, 0, MPI_COMM_WORLD);

            pack_up(buffer_send_up, u_actual);
            MPI_Send(buffer_send_up, n * kk, MPI_DOUBLE, id + npx, 0, MPI_COMM_WORLD);

            MPI_Recv(buffer_rec_down, n * kk, MPI_DOUBLE, id - npx, 0, MPI_COMM_WORLD, &status);
            unpack_down(buffer_rec_down, u_actual);

            MPI_Recv(buffer_rec_up, n * kk, MPI_DOUBLE, id + npx, 0, MPI_COMM_WORLD, &status);
            unpack_up(buffer_rec_up, u_actual);
        } else {
            MPI_Recv(buffer_rec_down, n * kk, MPI_DOUBLE, id - npx, 0, MPI_COMM_WORLD, &status);
            unpack_down(buffer_rec_down, u_actual);

            MPI_Recv(buffer_rec_up, n * kk, MPI_DOUBLE, id + npx, 0, MPI_COMM_WORLD, &status);
            unpack_up(buffer_rec_up, u_actual);

            pack_down(buffer_send_down, u_actual);
            MPI_Send(buffer_send_down, n * kk, MPI_DOUBLE, id - npx, 0, MPI_COMM_WORLD);

            pack_up(buffer_send_up, u_actual);
            MPI_Send(buffer_send_up, n * kk, MPI_DOUBLE, id + npx, 0, MPI_COMM_WORLD);
        }

    }
}


inline void solver(double *&u_actual, double *&u_next, int id, int size) {
    double u,t;
    int row,col;
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

    double buf_send_down[n * kk];
    double *buffer_send_down = &buf_send_down[0];

    double buf_send_up[n * kk];
    double *buffer_send_up = &buf_send_up[0];

    double buf_rec_down[n * kk];
    double *buffer_rec_down = &buf_rec_down[0];

    double buf_rec_up[n * kk];
    double *buffer_rec_up = &buf_rec_up[0];



    MPI_Barrier(MPI_COMM_WORLD);

//sending and receiving is chess law if id%2=0 it first send the receive, overvise, receive, then send
    t=MPI_Wtime();
    row=id/npx;
    col=id%npx;
//    std::cout<<"id: "<<id<<" col: "<<col<<" row: "<<row<<std::endl;
    for (int l = 1; l * dt <= 1; l++) {
        chatting_horizontal(id,row,col,u_actual,buffer_send_right,buffer_rec_right,buffer_send_left,buffer_rec_left,status);
        chatting_vertical(id,row,size,u_actual,buffer_send_up,buffer_rec_up,buffer_send_down,buffer_rec_down,status);

        for (unsigned k = 1; k < kk - 1; ++k) {
            for (unsigned j = 1; j < m - 1; ++j) {
                for (unsigned i = 1; i < n-1; ++i) {
                    u = get_u(u_actual, i, j, k) +
                        dt * (f((col * (n-2) + i) * hx, (row * (m-2) + j) * hy, k * hz)
                              + D0 * Lx(u_actual, i, j, k)
                              + D1 * Ly(u_actual, i, j, k)
                              + D2 * Lz(u_actual, i, j, k));
                    set_u(u_next, i, j, k, u);
                }
            }
        }
        swap(u_actual, u_next);
    }
    t=MPI_Wtime()-t;
    if (id==0) {
        std::cout << "time " << t  << std::endl;
    }
    get_mist(u_actual,id,row,col);
}