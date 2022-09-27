#include <iostream>
#include <cstdlib>
#include "simulation.h"

// g++ *.cpp -o simulation -O3 -std=c++17 -I. -Lmfem -lmfem -lrt

int main() {
    Simulation simulation(0.2, 0.2, 0.0015, 0.18e6, 1020.0, 1e-2, 10.0, 0.01);
    simulation.run();
        
    return 0;
}