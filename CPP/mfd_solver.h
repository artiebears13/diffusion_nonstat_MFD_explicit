#include <iostream>
#include <vector>

class Solver {
public:
    double D[3] = {0.25, 0.15, 0.1};
    int n, m, kk, t;
    double hx, hy, hz, dt;
    std::vector<double> *u_actual=new std::vector<double>();
    std::vector<double> *u_next=new std::vector<double>();

    //constructor
    Solver(int n, int m, int k, int t);

    //overload operator(): returns u_actual with i,j,k index
    double const &operator()(int i, int j, int k);

    //write to u_actual[i,j,k] value
    void set_value(int i, int j, int k, const double &value);

    //write to u_actual[i,j,k] value
    double get_u(int i, int j, int k) const;

    //count Lx
    inline double Lx(int i, int j, int k);

    //count Ly
    inline double Ly(int i, int j, int k);

    //count Lz
    inline double Lz(int i, int j, int k);

    //count f
    inline double f(int i, int j, int k);

    //swapping pointers of u_actual and u_next
    inline void swap();

    //return true if u_actual[i,j,k] is boundary
    bool check_boundary(int I, int J, int K);

    //check is dt is satisfy condition
    //if not - reread from keyboard
    void check_t();

    //calculate u on next layer and jump to next layer
    void calculate_next_u();

    //solver (found u on t=1)
    std::vector<double> * solve();

    //fill boundaries on x direction (x=0, x=N)
    void fill_boundary_x();

    //fill boundaries on y direction (y=0, y=N)
    void fill_boundary_y();

    //fill boundaries on z direction (z=0, z=N)
    void fill_boundary_z();

    //fill all boundary u with zeros
    void fill_boundaries();

    //calculates u on next layer
    void get_next_u();

    //write to u_next value
    void set_value_next(int i, int j, int k, const double &value);

    //count mistake versus analytic solve
    double mistake() const;

    //calculates analytic solve
    double get_analytic_u(int i, int j, int k) const;

};