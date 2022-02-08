//
// Created by shubham on 08.06.21.
//

#ifndef TEST_MPI_INTERFACE_H
#define TEST_MPI_INTERFACE_H

namespace Parallel {
    void MpiSolve(int argc, char **argv, float &parallel_runtime, int& print_rank);
}

#endif //TEST_MPI_INTERFACE_H
