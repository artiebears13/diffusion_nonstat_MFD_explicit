#include <iostream>
#include <vector>

class Solver {
public:
    double D[3] = {0.25, 0.15, 0.1};
    int n, m, kk, t;
    double hx, hy, hz, dt;
    std::vector<double> *u_actual=new std::vector<double>();
    std::vector<double> *u_next=new std::vector<double>();

    Solver(int n, int m, int k, int t);

    double const &operator()(int i, int j, int k);

    void set_value(int i, int j, int k, const float &value);

    const double get_u(int i, int j, int k);


    double Lx(int i, int j, int k);

    double Ly(int i, int j, int k);

    double Lz(int i, int j, int k);

    double f(int i, int j, int k);

    void swap();
};