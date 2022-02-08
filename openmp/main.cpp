//
// Created by shubham on 28.05.21.
//

#include <iostream>
#include <iomanip>
#include "utility.h"
#include "interface.h"

int main(int argc, char* argv[]){
    float sequential_runtime, parallel_runtime=1;
    int print_rank;

    Parallel::Solve(argc, argv, sequential_runtime, parallel_runtime, print_rank);
    // This is to handle MPI case
    if (print_rank==0){
        std::cout<<"DONE"<<std::endl;
    }
    return 0;
}