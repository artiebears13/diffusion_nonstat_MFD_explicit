#include "mfd_solver.h"
#include <cmath>

Solver::Solver(int N, int M, int K, int T) {
    this->n = N;
    this->m = M;
    this->kk = K;
    this->t = T;
    this->hx = 1. / (n - 1) * 1.0;
    this->hy = 1. / (m - 1) * 1.0;
    this->hz = 1. / (kk - 1) * 1.0;
    this->dt = 1. / T * 1.0;

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

bool Solver::check_boundary(int I, int J, int K) {
    if ((I == 0) || (I == this->n) || (J == 0) || (J == this->m) ||
        (K == 0) || (K == this->kk))
        return true;
    else
        return false;
}


void Solver::check_t() {
    std::cout << "dt = " << this->dt << " < " << 0.5 /
                                                 (this->D[0] / (pow(this->hx, 2)) + this->D[1] / (pow(this->hy, 2)) +
                                                  this->D[2] / (pow(this->hz, 2))) << std::endl;
    std::string answer;
    std::cout << "Do you want to continue: (y/n):";
    std::cin >> answer;

    if (answer == "y") {

    } else if (answer == "n") {
        std::cout << "enter new t:";
        std::cin >> this->t;
        this->dt = 1. / this->t;
        this->check_t();
    } else {
        std::cerr << "wrong input" << std::endl;
        exit(1);
    }
}

std::vector<double> *Solver::solve() {
    this->check_t();
    for (int i = 1; i <= this->t; i++) {
//        print();
//        std::cout<<i<<std::endl;
        calculate_next_u();
    }
    return this->u_actual;
}

void Solver::fill_boundary_x() {
    for (int j = 0; j < this->m; ++j) {
        for (int k = 0; k < this->kk; ++k) {
            set_value(0, j, k, 0.);
            set_value(this->n - 1, j, k, 0.);
        }
    }
}

void Solver::fill_boundary_y() {
    for (int i = 0; i < this->n; ++i) {
        for (int k = 0; k < this->kk; ++k) {
            set_value(i, 0, k, 0.);
            set_value(i, this->m - 1, k, 0.);
        }
    }
}


void Solver::fill_boundary_z() {
    for (int i = 0; i < this->n; ++i) {
        for (int j = 0; j < this->m; ++j) {
            set_value(i, j, 0, 0.);
            set_value(i, j, this->kk - 1, 0.);
        }
    }
}


void Solver::fill_boundaries() {
    fill_boundary_x();
    fill_boundary_y();
    fill_boundary_z();
}


double Solver::get_analytic_u(int i, int j, int k) const {
    double x = i * this->hx;
    double y = j * this->hy;
    double z = k * this->hz;
    return sin(M_PI * x) * sin(M_PI * y) * sin(M_PI * z) *
           (1 - exp(-((this->D[0] + this->D[1] + this->D[2]) * pow(M_PI, 2) * 1)));
}


void Solver::set_value(int i, int j, int k, const double &value) {
    this->u_actual->at(i * this->n * this->m + j * this->kk + k) = value;
}

void Solver::set_value_next(int i, int j, int k, const double &value) {
    this->u_next->at(i * this->n * this->m + j * this->kk + k) = value;
}

double Solver::get_u(int i, int j, int k) const {
    return this->u_actual->at(i * this->n * this->m + j * this->kk + k);
}


inline double Solver::Lx(int i, int j, int k) {
    return (get_u(i - 1, j, k) - 2 * get_u(i, j, k) + get_u(i + 1, j, k))
           / (this->hx * this->hx);
}

inline double Solver::Ly(int i, int j, int k) {
    return (get_u(i, j - 1, k) - 2 * get_u(i, j, k) + get_u(i, j + 1, k))
           / (this->hy * this->hy);
}

inline double Solver::Lz(int i, int j, int k) {
    return (get_u(i, j, k - 1) - 2 * get_u(i, j, k) + get_u(i, j, k + 1))
           / (this->hz * this->hz);
}

inline double Solver::f(int i, int j, int k) {
    double x = i * this->hx;
    double y = j * this->hy;
    double z = k * this->hz;
    return ((this->D[0] + this->D[1] + this->D[2])
            * (M_PI * M_PI) * sin(M_PI * x) * sin(M_PI * y)
            * sin(M_PI * z));
}

void Solver::swap() {
    u_actual->swap(*u_next);
}

void Solver::get_next_u() {
    fill_boundaries();
    double u;
    for (int i = 1; i < this->n - 1; ++i) {
        for (int j = 1; j < this->m - 1; ++j) {
            for (int k = 1; k < this->kk - 1; ++k) {
                u = get_u(i, j, k) +
                    this->dt * (f(i, j, k)
                          + D[0] * Lx(i, j, k)
                          + D[1] * Ly(i, j, k)
                          + D[2] * Lz(i, j, k));
//                std::cout << "dt: " << this->dt << " get_u: " << get_u(i, j, k) << " " << " f(i, j, k): " << f(i, j, k) <<
//                          " D[0] * Lx(i, j, k): " << +D[0] * Lx(i, j, k) << " D[1] * Ly(i, j, k): "
//                          << D[1] * Ly(i, j, k)
//                          << " D[2] * Lz(i, j, k): " << D[2] * Lz(i, j, k) << " u: "
//                          << u << std::endl;
                set_value_next(i, j, k, u);
            }
        }
    }
}

inline void Solver::calculate_next_u() {
    get_next_u();
    swap();
}

double Solver::mistake() const {
    double mist = 0.;
    for (int i = 1; i < this->n - 1; ++i) {
        for (int j = 1; j < this->m - 1; ++j) {
            for (int k = 1; k < this->kk - 1; ++k) {
//                std::cout<<std::abs(this->get_u(i, j, k) - this->get_analytic_u(i, j, k))<< std::endl;
                if (std::abs(this->get_u(i, j, k) - this->get_analytic_u(i, j, k)) > mist) {
                    mist = std::abs(get_u(i, j, k) - get_analytic_u(i, j, k));
                }
            }
        }
    }
    return mist;
}


