#include <iostream>
#include <cmath>
#include <chrono>
#include <omp.h>

#define N 65
#define n N
#define m N
#define kk N
#define hx 1./(n-1)*1.0
#define hy 1. / (m - 1) * 1.0
#define hz 1. / (kk - 1) * 1.0
#define D0 0.25
#define D1 0.15
#define D2 0.1
//#define dt 0.9*(0.5 /(D0 / (hx*hx) + D1 / (hy*hy) + D2 / (hz*hz)))
#define dt (0.45 /(D0 / (hx*hx) + D1 / (hy*hy) + D2 / (hz*hz)))
//#define dt 0.001
#define get_analytic_u(x, y, z) (sin(M_PI * x) * sin(M_PI * y) * sin(M_PI * z) *(1 - exp(-((D0 + D1 + D2) * M_PI*M_PI * 1))))
#define f(x, y, z) ((D0+D1+D2)*(M_PI*M_PI*sin(M_PI*x))* sin(M_PI * y) * sin(M_PI * z))


inline double get_u(double *const u, int i, int j, int k) {
    return *(u + i * n * m + j * kk + k);
}

inline void set_u(double *u, int i, int j, int k, double value) {
    *(u + i * n * m + j * kk + k) = value;
}

inline void set_value_next(double *u_next, int i, int j, int k, const double value) {
    *(u_next + i * n * m + j * kk + k) = value;
}

inline void fill_boundary_x(double *u_actual) {
    for (int j = 0; j < m; ++j) {
        for (int k = 0; k < kk; ++k) {
            set_u(u_actual, 0, j, k, 0.);
            set_u(u_actual, n - 1, j, k, 0.);
        }
    }
}

inline void fill_boundary_y(double *const &u_actual) {
    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < kk; ++k) {
            set_u(u_actual, i, 0, k, 0.);
            set_u(u_actual, i, m - 1, k, 0.);
        }
    }
}


void fill_boundary_z(double *&u_actual) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            set_u(u_actual, i, j, 0, 0.);
            set_u(u_actual, i, j, kk - 1, 0.);
        }
    }
}

inline double Lx(double *u_actual, int i, int j, int k) {
    return (get_u(u_actual, i - 1, j, k) - 2 * get_u(u_actual, i, j, k) + get_u(u_actual, i + 1, j, k))
           / (hx * hx);
}

inline double Ly(double *u_actual, int i, int j, int k) {
    return (get_u(u_actual, i, j - 1, k) - 2 * get_u(u_actual, i, j, k) + get_u(u_actual, i, j + 1, k))
           / (hy * hy);
}

inline double Lz(double *u_actual, int i, int j, int k) {
    return (get_u(u_actual, i, j, k - 1) - 2 * get_u(u_actual, i, j, k) + get_u(u_actual, i, j, k + 1))
           / ((hz) * (hz));
}

void swap(double *&u_actual, double *&u_next) {
    auto tmp = u_actual;
    u_actual = u_next;
    u_next = tmp;
}

int main() {
    int len = n * m * kk;
    double *u_actual = new double[len];
    double *u_next = new double[len];
    for (int i = 1; i < n; ++i) {
        for (int j = 1; j < m; ++j) {
            for (int k = 1; k < kk; ++k) {
                set_u(u_actual, i, j, k, 0.);
//                set_u(u_next,i,j,k,0.);
            }
        }
    }

//    set_u(u_actual, 1, 0, 0, 10.);
//    set_u(u_next, 1, 0, 0, 0.);
//    std::cout << get_u(u_next, 1, 0, 0) << std::endl;
//    swap(u_actual, u_next);
//    std::cout << get_u(u_next, 1, 0, 0) << std::endl;
//    std::cout << get_u(u_actual, 1, 0, 0) << std::endl;

//    clock_t tStart = clock();
   auto t = omp_get_wtime();
    for (int l = 1; l * dt <= 1; l++) {
//        print();
//        std::cout<<i<<std::endl;
        double u;
        for (int i = 0; i < n; ++i) {
#pragma omp parallel for
            for (int j = 0; j < m; ++j) {
                for (int k = 0; k < kk; ++k) {
                    if (i == 0) { set_u(u_actual, 0, j, k, 0.); }
                    else if (i == n - 1) { set_u(u_actual, i, j, k, 0.); }
                    else if (j == 0) { set_u(u_actual, i, 0, k, 0.); }
                    else if (j == m - 1) { set_u(u_actual, i, j, k, 0.); }

                    else if (k == 0) { set_u(u_actual, i, j, 0, 0.); }
                    else if (k == kk - 1) { set_u(u_actual, i, j, k, 0.); }
                    else {
                        u = get_u(u_actual, i, j, k) +
                            dt * (f(i * hx, j * hy, k * hz)
                                  + D0 * Lx(u_actual, i, j, k)
                                  + D1 * Ly(u_actual, i, j, k)
                                  + D2 * Lz(u_actual, i, j, k));
//                    std::cout << "dt: " << dt << " get_u: " << get_u(u_actual, i, j, k) << " " << " f(i, j, k): "
//                              << f(i * hx, j * hy, k * hz) <<
//                              " D0 * Lx(i, j, k): " << +D0 * Lx(u_actual, i, j, k) << " D1 * Ly(i, j, k): "
//                              << D1 * Ly(u_actual, i, j, k)
//                              << " D2 * Lz(i, j, k): " << D2 * Lz(u_actual, i, j, k) << " u: "
//                              << u << std::endl;

                        set_u(u_next, i, j, k, u);
                    }
                }
            }
        }
        swap(u_actual, u_next);
    }

 //   auto time = (double) (clock() - tStart) / CLOCKS_PER_SEC;
   t=omp_get_wtime() - t;
    std::cout << "time: " << t;

    double mist = 0.;
    for (int i = 0; i < n - 1; ++i) {
        for (int j = 0; j < m; ++j) {
            for (int k = 0; k < kk; ++k) {
//                std::cout <<"mist: "<< mist << std::endl;
                if (std::abs(get_u(u_actual, i, j, k) - get_analytic_u(i * hx, j * hy, k * hx)) > mist) {
                    mist = std::abs(get_u(u_actual, i, j, k) - get_analytic_u(i * hx, j * hy, k * hx));
                }
            }
        }
    }

    std::cout << "mist: " << mist << std::endl;
    std::cout << "dt: " << dt << std::endl;


    return 0;
}

