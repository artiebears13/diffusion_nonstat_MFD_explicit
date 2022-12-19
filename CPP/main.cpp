#include <iostream>
//#include "mfd_solver.h" //uncomment if compilate with cmakelists
//#include "mfd_solver.cpp" //uncomment if compilate manually
#include <chrono>

int main() {
    Solver a(11, 11, 11, 1650);
//    a.set_value(1,1,1,10.);
//    a.swap();
//    a.swap();
//    std::cout << a(1, 1, 1) << " "<< a(0,0,0);
    clock_t tStart = clock();
    /* Do your stuff here */

    a.solve();
    auto time = (double) (clock() - tStart) / CLOCKS_PER_SEC;
    std::cout << "time: " << time;

    float mist = a.mistake();
    std::cout << "mist: " << mist << std::endl;


    return 0;
}