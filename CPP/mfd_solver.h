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

    void set_value(int i, int j, int k, const double &value);

    double get_u(int i, int j, int k) const;


    inline double Lx(int i, int j, int k);

    inline double Ly(int i, int j, int k);

    inline double Lz(int i, int j, int k);

    inline double f(int i, int j, int k);

    inline void swap();

    bool check_boundary(int I, int J, int K);

    void check_t();

    void calculate_next_u();

    std::vector<double> * solve();

    void fill_boundary_x();

    void fill_boundary_y();

    void fill_boundary_z();

    void fill_boundaries();

    void get_next_u();

    void set_value_next(int i, int j, int k, const double &value);

    double mistake() const;

    double get_analytic_u(int i, int j, int k) const;

    void print();
};