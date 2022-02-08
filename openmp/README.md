# OpenMP

## serial optimizations

- replaced pow function calls in CaseSelector with bit-shift operations, which can be evaluated by the compiler at compile-time
- switched to custom octree data structure (OurOctreeNode and OurOctree)
- OurOctree contains a field of all nodes in the entire tree (stored in a std::vector)
- individual nodes do not store pointers, but indices to their children. this leads to less memory allocations during tree construction, and less memory fragmentation of the nodes in the octree
- the octree data structure (i.e. the contained std::vector) is reused between runs to avoid unnecessary memory allocations

## parallelization strategy

We chose to focus on parallelizing the treeInteract part of the algorithm, since according to our performance measurements around 80% of the time was spent there after our serial optimizations.

With 24 threads we achieved a maximum speedup of ~5 on the submission server. In the end we decided to use the guided schedule, since it outperformed the static variant in our measurements.
