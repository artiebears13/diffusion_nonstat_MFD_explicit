#include <iostream>
#include "mfd_solver.h"
#include <chrono>

int main() {
    Solver a(11, 11, 11, 1000);
//    a.set_value(1,1,1,10.);
//    a.swap();
//    a.swap();
//    std::cout << a(1, 1, 1) << " "<< a(0,0,0);
    clock_t tStart = clock();
    /* Do your stuff here */

    a.solve();
    printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

    float mist = a.mistake();
    std::cout << "mist: " << mist << std::endl;


    return 0;
}