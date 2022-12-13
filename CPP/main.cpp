#include <iostream>
#include "mfd_solver.h"

int main() {
    Solver a(10, 10, 10, 1000);
//    a.set_value(1,1,1,10.);
//    a.swap();
//    a.swap();
//    std::cout << a(1, 1, 1) << " "<< a(0,0,0);
    a.solve();
    float mist = a.mistake();
    std::cout << mist << std::endl;

    return 0;
}