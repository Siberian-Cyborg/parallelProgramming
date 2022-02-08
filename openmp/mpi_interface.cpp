//
// Created by shubham on 08.06.21.
//

#include "mpi_interface.h"
#include "utility.h"
#include <iostream>
#include <chrono>
#include "mpibh.h"
#include "constant.h"
#include "mpi.h"

void Parallel::MpiSolve(int argc, char **argv, float &parallel_runtime, int& print_rank) {
    int numBodies, numSteps;
    int rank, size;
    Body *bodies = nullptr;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (rank==0) {
        Utility::ParseInputs(argc, argv, &numBodies, &numSteps);
        bodies = new Body[numBodies];

        Utility::Initialize(bodies, numBodies, 0);
        print_rank = rank;
    }
    else{
        print_rank = 1;
    }

    auto start = std::chrono::steady_clock::now();
    MPI::Solve(bodies, numBodies, numSteps);
    if (rank==0) {
        auto end = std::chrono::steady_clock::now();
        parallel_runtime = std::chrono::duration<float>(end - start).count();
    }

    if (rank==0){
        Utility::PrintBodies(bodies, numBodies);
        delete[] bodies;
    }
    MPI_Finalize();
}