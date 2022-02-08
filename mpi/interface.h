//
// Created by shubham on 08.06.21.
//

#ifndef BARNES_HUT_INTERFACE_H
#define BARNES_HUT_INTERFACE_H

namespace Parallel{
    void Solve(int argc, char** argv,
               float& sequential_runtime, float& parallel_runtime, int& print_rank);
}

#endif //BARNES_HUT_INTERFACE_H
