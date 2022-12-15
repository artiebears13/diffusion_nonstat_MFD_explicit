#include <iostream>
#include <cmath>
#include <chrono>
#include <array>
#include <omp.h>


#define Nx 120
#define Ny 60
#define Nz 20
#define n Nx
#define m Ny
#define kk Nz
#define len n * m * kk
#define hx 1./(n-1)*1.0
#define hy 1. / (m - 1) * 1.0
#define hz 1. / (kk - 1) * 1.0
#define D0 0.25
#define D1 0.15
#define D2 0.1
#define dt 0.9*(0.5 /(D0 / (hx*hx) + D1 / (hy*hy) + D2 / (hz*hz)))
//#define dt (0.45 /(D0 / (hx*hx) + D1 / (hy*hy) + D2 / (hz*hz)))
//#define dt 0.001
#define get_analytic_u(x, y, z) (sin(M_PI * x) * sin(M_PI * y) * sin(M_PI * z) *(1 - exp(-((D0 + D1 + D2) * M_PI*M_PI * 1))))
#define f(x, y, z) ((D0+D1+D2)*(M_PI*M_PI*sin(M_PI*x))* sin(M_PI * y) * sin(M_PI * z))


inline double get_u(double *const &u, int i, int j, int k) {
    return *(u + k * n * m + j * n + i);
}

inline void set_u(double * const &u, int i, int j, int k, double const &value) {
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

void swap(double * &u_actual, double * &u_next) {


    double *temp = u_actual;
    u_actual = u_next;
    u_next = temp;

}

int main() {

    double u;
    double u_act[len];
    double u_n[len];

    double *u_actual = &u_act[0];
    double *u_next = &u_n[0];


//    double* u_actual=u_act.data();
//    double *u_next=u_n.data();

    for (int k = 0; k < kk; ++k) {
        for (int j = 0; j < m; ++j) {
            for (int i = 0; i < n; ++i) {

                set_u(u_actual, i, j, k, 0.);

            }
        }
    }

//    std::cout<< get_u(u_actual,5,5,5)<<std::endl;

//    double *u_actual = u_act.();
//    double *u_next = u_n.data();




//    set_u(u_actual, 1, 0, 0, 10.);
//    set_u(u_next, 1, 0, 0, 0.);
//    std::cout << get_u(u_next, 1, 0, 0) << std::endl;
//    swap(u_actual, u_next);
//    std::cout << get_u(u_next, 1, 0, 0) << std::endl;
//    std::cout << get_u(u_actual, 1, 0, 0) << std::endl;
//-----------------------------------------------------------------------------
    clock_t tStart = clock();
//   auto t = omp_get_wtime();
    for (int l = 1; l * dt <= 1; l++) {
//        print();
//        std::cout<<i<<std::endl;

        for (int k = 0; k < kk; ++k) {
//#pragma omp parallel for
            for (int j = 0; j < m; ++j) {
                for (int i = 0; i < n; ++i) {
                    if (i == 0) { set_u(u_next, i, j, k, 0.); }
                    else if (i == n - 1) { set_u(u_next, i, j, k, 0.); }
                    else if (j == 0) { set_u(u_next, i, j, k, 0.); }
                    else if (j == m - 1) { set_u(u_next, i, j, k, 0.); }

                    else if (k == 0) { set_u(u_next, i, j, k, 0.); }
                    else if (k == kk - 1) { set_u(u_next, i, j, k, 0.); }
                    else {
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
                    }
                }
            }
        }
//        std::cout << "BEFORE: " << "u_actu: " << u_actual[93] << " u_next: " << u_next[93] << std::endl;

        swap(u_actual, u_next);
//        std::cout << " AFTER: " << "u_actu: " << u_actual[93] << " u_next: " << u_next[93] << std::endl;
    }

    double time = (double) (clock() - tStart) / CLOCKS_PER_SEC;
//   t=omp_get_wtime() - t;
    std::cout << "time: " << time << std::endl;

    double mist = 0.;
    for (int k = 0; k < kk; ++k) {
        for (int j = 0; j < m; ++j) {
            for (int i = 0; i < n; ++i) {
//                std::cout <<"mist: "<< mist << std::endl;



//                std::cout<<"numbering: "<< get_u(u_actual, i, j, k)<<" analytic: " <<get_analytic_u(i * hx, j * hy, k * hx)<<std::endl;

                if (std::abs(get_u(u_actual, i, j, k) - get_analytic_u(i * hx, j * hy, k * hz)) > mist) {
                    mist = std::abs(get_u(u_actual, i, j, k) - get_analytic_u(i * hx, j * hy, k * hz));
                }
            }
        }
    }

    std::cout << "mist: " << mist << std::endl;
    std::cout << "dt: " << dt << std::endl;


    return 0;
}

