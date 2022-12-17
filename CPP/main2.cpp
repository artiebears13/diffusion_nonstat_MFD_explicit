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
#define kmax n * m * kk - m*n
#define jmax m * n-n
#define K kmax+m*n
#define J jmax+n
#define hx 1./(n-1)*1.0
#define hy 1. / (m - 1) * 1.0
#define hz 1. / (kk - 1) * 1.0
#define D0 0.25
#define D1 0.15
#define D2 0.1
#define dt 0.9*(0.5 /(D0 / (hx*hx) + D1 / (hy*hy) + D2 / (hz*hz)))
#define get_analytic_u(x, y, z) (sin(M_PI * x) * sin(M_PI * y) * sin(M_PI * z) *(1 - exp(-((D0 + D1 + D2) * M_PI*M_PI * 1))))
#define f(x, y, z) ((D0+D1+D2)*(M_PI*M_PI*sin(M_PI*x))* sin(M_PI * y) * sin(M_PI * z))


inline double get_u(double *const &u, int i, int j, int k) {
    return *(u + k * n * m + j * n + i);
}

inline void set_u(double *const &u, int i, int j, int k, double const &value) {
    *(u + k * n * m + j * n + i) = value;
}


//inline double Lx(double * index) {
////    std::cout<<"1 "<<(i + 1) + j + k<<std::endl;
//    return
//}
//
//inline double Ly(double * index) {
////    std::cout<<"1 "<<(i) + (j + n) + k<<std::endl;
//    return
//}
//
//inline double Lz(double * index) {
////    std::cout<<"1 "<<(i) + (j) + k + n*m<<std::endl;
//    return ;
//}

void swap(double *&u_actual, double *&u_next) {


    double *temp = u_actual;
    u_actual = u_next;
    u_next = temp;

}

int main() {


    double u_act[len];
    double u_n[len];
    int index;

    double *u_actual = &u_act[0];
    double *u_next = &u_n[0];


    for (int k = 0; k < kk; ++k) {
        for (int j = 0; j < m; ++j) {
            for (int i = 0; i < n; ++i) {
                set_u(u_actual, i, j, k, 0.);
            }
        }
    }


//-----------------------------------------------------------------------------
    clock_t tStart = clock();
//   auto t = omp_get_wtime();
    for (int l = 1; l * dt <= 1; l++) {

        for(int k = 0; k < K; k += n * m){
//#pragma omp parallel for
            for(int j = 0; j <J ; j += n){
                for(int i = 0; i < n; ++i){
                    index=i+j+k;
//                     std::cout << i<<" "<<j<<" "<<k<<std::endl;
//                    double* index = u_next + i + j + k;
//                    double* index_a = u_actual + i + j + k;

//                    std::cout << " i=" << i << " j=" << j << " k=" << k << std::endl;
                     if (i == 0 or (i == n - 1) or (j == 0) or (j == jmax) or (k == 0) or (k == kmax)){ *(u_next+index) = 0.; }
                    else {
                        *(u_next + index) = *(u_actual+index) +
                                            dt * (f(i  * hx, j / n * hy, k /n/m* hz)
                                                  + D0 * ((*(u_actual+index-1) - 2 * (*(u_actual+index)) +
                                                          *(u_actual+index+1))/ (hx * hx)
                                                  + D1 * ((*(u_actual+index-n) - 2 * (*(u_actual+index))) +
                                                          *(u_actual+index+n)) / (hy * hy))
                                                  + D2 * (((*(u_actual+index- n*m)) - 2 * (*(u_actual+index)) +
                                                          (*(u_actual+index+ n*m))) / ((hz) * (hz))));
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

