#include <iostream>
#include "mfd_solver.h"

int main() {
    Solver a(3, 3, 3, 3);
    a.set_value(1,1,1,10.);
    a.swap();
    a.swap();
    std::cout << a(1, 1, 1) << " "<< a(0,0,0);


    return 0;
}