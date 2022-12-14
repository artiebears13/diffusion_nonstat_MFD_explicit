#include <iostream>
#include <cmath>

#define N 11
#define n N
#define m N
#define kk N
#define hx 1./(n-1)*1.0
#define hy 1. / (m - 1) * 1.0
#define hz 1. / (kk - 1) * 1.0
#define D0 0.25
#define D1 0.15
#define D2 0.1
#define dt 0.9*(0.5 /(D0 / (hx*hx) + D1 / (hy*hy) + D2 / (hz*hz)))
#define get_analytic_u(x,y,z) (sin(M_PI * x) * sin(M_PI * y) * sin(M_PI * z) *(1 - exp(-((D0 + D1 + thisD2) * M_PI*M_PI * 1))))


inline double get_u(double *u, int i, int j, int k) {
    return *(u + i * n * m + j * kk + k);
}

inline void set_u(double *u, int i, int j, int k, double value) {
    *(u + i * n * m + j * kk + k) = value;
}

inline void fill_boundary_x(double *u_actual) {
    for (int j = 0; j < m; ++j) {
        for (int k = 0; k < kk; ++k) {
            set_u(u_actual, 0, j, k, 0.);
            set_u(u_actual, n - 1, j, k, 0.);
        }
    }
}

inline void fill_boundary_y(double *u_actual) {
    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < kk; ++k) {
            set_u(u_actual, i, 0, k, 0.);
            set_u(u_actual, i, m - 1, k, 0.);
        }
    }
}


void fill_boundary_z(double *u_actual) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            set_u(u_actual, i, j, 0, 0.);
            set_u(u_actual, i, j, kk - 1, 0.);
        }
    }
}

inline double Lx(double *u_actual, int i, int j, int k) {
    return (get_u(u_actual, i - 1, j, k) - 2 * get_u(u_actual, i, j, k) + get_u(u_actual, i + 1, j, k))
    / (hx * hx );
}

inline double Ly(double *u_actual, int i, int j, int k) {
    return (get_u(u_actual, i, j - 1, k) - 2 * get_u(u_actual, i, j, k) + get_u(u_actual, i, j + 1, k))
    / (hy*hy);
}

inline double Lz(double *u_actual, int i, int j, int k) {
    return (get_u(u_actual, i, j, k - 1) - 2 * get_u(u_actual, i, j, k) + get_u(u_actual, i, j, k + 1))
    / ((hz)*(hz));
}

int main() {
    int len = n * m * kk;
    double *u_actual = new double[len];
    double *u_next = new double[len];


    return 0;
}