#include <iostream>
#include <cmath>
#include <chrono>
//#include "to_vtk.hpp"   //<----uncomment to build .vtk

#define get_analytic_u(x, y, z) (sin(M_PI * x) * sin(M_PI * y) * sin(M_PI * z) *(1 - exp(-((D0 + D1 + D2) * M_PI*M_PI * 1))))
//#define f(x, y, z) ((D0+D1+D2)*(M_PI*M_PI*sin(M_PI*x))* sin(M_PI * y) * sin(M_PI * z))

const unsigned Nx=12;  //<---your value here
const unsigned Ny=10;  //<---your value here
const unsigned Nz=10;  //<---your value here
const unsigned n =Nx;
const unsigned m =Ny;
const unsigned kk= Nz;
const unsigned len= n * m * kk;
const double hx= 1./(n-1)*1.0;
const double hy= 1. / (m - 1) * 1.0;
const double hz= 1. / (kk - 1) * 1.0;
const double D0= 0.25;
const double D1= 0.15;
const double D2= 0.1;
const double dt=0.9*(0.5 /(D0 / (hx*hx) + D1 / (hy*hy) + D2 / (hz*hz)));





inline double f(double x, double y, double z) {return((D0+D1+D2)*(M_PI*M_PI*sin(M_PI*x))* sin(M_PI * y) * sin(M_PI * z));}


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


inline void get_mist(double* &u_actual){
    double mist = 0.;
    for (int k = 0; k < kk; ++k) {
        for (int j = 0; j < m; ++j) {
            for (int i = 0; i < n; ++i) {
                if (std::abs(get_u(u_actual, i, j, k) - get_analytic_u(i * hx, j * hy, k * hz)) > mist) {
                    mist = std::abs(get_u(u_actual, i, j, k) - get_analytic_u(i * hx, j * hy, k * hz));
                }
            }
        }
    }

    std::cout << "mist: " << mist << std::endl;
    std::cout << "dt: " << dt << std::endl;
}

int main() {

    double u;
    double u_act[len];
    double u_n[len];

    double *u_actual = &u_act[0];
    double *u_next = &u_n[0];


    for (int k = 0; k < kk; ++k) {
        for (int j = 0; j < m; ++j) {
            for (int i = 0; i < n; ++i) {
                set_u(u_actual, i, j, k, 0.);
                set_u(u_next, i, j, k, 0.);
            }
        }
    }

    clock_t tStart = clock();
    for (int l = 1; l * dt <= 1; l++) {
        for (unsigned k = 1; k < kk-1; ++k) {
            for (unsigned j = 1; j < m-1; ++j) {
                for (unsigned i = 1; i < n-1; ++i) {
//                    if (i == 0) { set_u(u_next, i, j, k, 0.); }
//                    else if (i == n - 1) { set_u(u_next, i, j, k, 0.); }
//                    else if (j == 0) { set_u(u_next, i, j, k, 0.); }
//                    else if (j == m - 1) { set_u(u_next, i, j, k, 0.); }
//
//                    else if (k == 0) { set_u(u_next, i, j, k, 0.); }
//                    else if (k == kk - 1) { set_u(u_next, i, j, k, 0.); }
//                    else {
                    u = get_u(u_actual, i, j, k) +
                        dt * (f(i * hx, j * hy, k * hz)
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
    }
//    write_to_file(u_actual,n,m,kk);      //<----uncomment to build .vtk
    double time = (double) (clock() - tStart) / CLOCKS_PER_SEC;
    std::cout << "time: " << time << std::endl;
    get_mist(u_actual);
    return 0;
}

