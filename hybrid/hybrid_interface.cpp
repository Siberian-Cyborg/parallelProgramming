//
// Created by shubham on 08.06.21.
//

#include "hybrid_interface.h"
#include "utility.h"
#include <iostream>
#include <chrono>
#include "hybridbh.h"
#include "constant.h"
#include "mpi.h"
#include <omp.h>

void Parallel::HybridSolve(int argc, char **argv, float &parallel_runtime, int& print_rank) {

    int numBodies, numSteps;
    int rank, size;
    Body *bodies = nullptr;
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    if(provided < MPI_THREAD_FUNNELED)
    {
        printf("The threading support level is lesser than that demanded.\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    else
    {
        printf("The threading support level corresponds to that demanded.\n");
    }
 
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
    Hybrid::Solve(bodies, numBodies, numSteps);
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