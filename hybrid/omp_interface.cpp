//
// Created by shubham on 08.06.21.
//

#include "omp_interface.h"
#include "utility.h"
#include <iostream>
#include <chrono>
#include "ompbh.h"
#include <omp.h>
#include "constant.h"

void Parallel::OmpSolve(int argc, char **argv, float &parallel_runtime) {
    int numBodies, numSteps;
    Utility::ParseInputs(argc, argv, &numBodies, &numSteps);
    Body *bodies = new Body[numBodies];

    Utility::Initialize(bodies, numBodies, 0);

    auto start = std::chrono::steady_clock::now();
    OMP::Solve(bodies, numBodies, numSteps);
    auto end = std::chrono::steady_clock::now();
    parallel_runtime = std::chrono::duration<float>(end - start).count();
    Utility::PrintBodies(bodies, numBodies);
    delete[] bodies;
}