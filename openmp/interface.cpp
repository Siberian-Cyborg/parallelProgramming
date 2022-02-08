//
// Created by shubham on 08.06.21.
//

#include "interface.h"
#include "utility.h"
#include <iostream>
#include "build.h"
#include "serial.h"
#include <chrono>

#ifdef BUILDOMP
#include "omp_interface.h"
#endif

#ifdef BUILDMPI
#include "mpi_interface.h"
#endif

#ifdef BUILDHYBRID
#include "hybrid_interface.h"
#endif

void Parallel::Solve(int argc, char** argv,
                     float& sequential_runtime, float& parallel_runtime, int& print_rank){


#if defined(BUILDOMP)
    Parallel::OmpSolve(argc, argv, parallel_runtime);
#endif

#if defined(BUILDMPI)
    Parallel::MpiSolve(argc, argv, parallel_runtime, print_rank);
#endif

#if defined(BUILDHYBRID)
    Parallel::HybridSolve(argc, argv, parallel_runtime, print_rank);
#endif

#if defined(BUILDOMP) || (defined(BUILDSERIAL) && !defined(BUILDHYBRID) && !defined(BUILDMPI))
    print_rank = 0;
#endif

#ifdef BUILDSERIAL
    Serial::SerialSolve(argc, argv,sequential_runtime, print_rank);
#endif

}