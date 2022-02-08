# parallelProgramming
Parallelization of the Barnes Hut algorithm using OMP, MPI, SIMD and Hybrid approaches.

## Barnes Hut
The Barnes Hut algorithm is an efficient way of sloving the N-Body Simulation problem in which N bodies interact with each other resulting in a complexity of N^2. Barnes-Hut uses a Octree data structure to approximate the interaction between the bodies thus reducing the complexity to Nlog(N) (which is waaaay less than N^2).
Check out more detailed explanations on YouTube. Here is a nice animation of the Barnes Hut algorithm:

https://www.youtube.com/watch?v=rkPxkdBCvEA

## Parallelization
The point of this project was to accelerate the Barnes-Hut algorithm by running it in parallel on multiple threads/cores/sockets on a server. For that different frameworks were used. At first the parallelization was done using OpenMP (https://en.wikipedia.org/wiki/OpenMP). Afterwards MPI (https://en.wikipedia.org/wiki/Message_Passing_Interface) was used. Then a hybrid approach utilizing both frameworks was implemented. Finally some SIMD operations were added to further increase the speedup. 

These presentation slides summarize all techniques used:

[BarnesHutParallel.pdf](https://github.com/Siberian-Cyborg/parallelProgramming/files/8019563/BarnesHutParallel.pdf)
