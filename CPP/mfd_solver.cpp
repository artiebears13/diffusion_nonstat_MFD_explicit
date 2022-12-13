#include "mfd_solver.h"
#include <cmath>

Solver::Solver(int N, int M, int K, int T) {
    this->n = N;
    this->m = M;
    this->kk = K;
    this->t = T;
    this->hx = 1 / (n - 1);
    this->hy = 1 / (m - 1);
    this->hz = 1 / (kk - 1);
    this->dt = 1 / T;

    for (int i = 0; i < this->n; ++i) {
        for (int j = 0; j < this->m; ++j) {
            for (int k = 0; k < this->kk; ++k) {
                this->u_actual->push_back(0.);
                this->u_next->push_back(0.);
            }
        }
    }
}

double const &Solver::operator()(int i, int j, int k) {
    return this->u_actual->at(i * this->n * this->m + j * this->kk + k);
}


void Solver::set_value(int i, int j, int k, const float &value) {
    this->u_actual->at(i * this->n * this->m + j * this->kk + k) = value;
}

const double Solver::get_u(int i, int j, int k) {
    return this->u_actual->at(i * this->n * this->m + j * this->kk + k);
}


double Solver::Lx(int i, int j, int k) {
    return (get_u(i - 1, j, k) - 2 * get_u(i, j, k) + get_u(i + 1, j, k))
           / (this->hx * this->hx);
}

double Solver::Ly(int i, int j, int k) {
    return (get_u(i, j - 1, k) - 2 * get_u(i, j, k) + get_u(i, j + 1, k))
           / (this->hy * this->hy);
}

double Solver::Lz(int i, int j, int k) {
    return (get_u(i, j, k - 1) - 2 * get_u(i, j, k) + get_u(i, j, k + 1))
           / (this->hz * this->hz);
}

double Solver::f(int i, int j, int k) {
    double x=i*this->hx;
    double y=j*this->hy;
    double z=k*this->hz;
    return ((this->D[0] + this->D[1] + this->D[2])
            * (M_PI * M_PI) * sin(M_PI * x) * sin(M_PI * y)
            * sin(M_PI * z));
}

void Solver::swap(){
    u_actual->swap(*u_next);
}