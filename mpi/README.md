# MPI

## parallelization strategy

For MPI we parallelized all three parts of the algorithm: octree construction, force computation and position updates.

Before the computation begins, the root MPI process (rank #0) broadcasts the body positions to all other processes.

### octree construction

Octree construction happens in parallel on all MPI processes.

The octree has 1 entry on the level 0, 8 entries on level 1, 64 entries on level 2, etc.
We now select the level, which has at least as many entries as there are MPI processes (we call this the *root cell level*).

Now all nodes of the root cell level are split equally among the MPI processes, resulting in one set of assigned root nodes per MPI process.

During tree construction, all processes consider all bodies. 

However, we change the octree construction algorithm step #4. In addition to checking whether the child node is empty,
we also check if the child node's level is the root node level. If so, we only insert a new node, 
if the node is in the set of assigned nodes for the active MPI process. Otherwise, we abort the insertion, without inserting the node.

Starting from the root cell level the construction of the octree is spread equally among the MPI processes. 
Before that, each MPI process is building the exact same tree.

After tree construction, each MPI process broadcasts its octree to all other processes, so all processes then have the entire tree.
The broadcasting is done using non-blocking MPI_Isend, MPI_Iprobe and MPI_Irecv calls. 
By using MPI_Iprobe to dynamically query the data size, we avoid the overhead of an additional message for broadcasting the tree size.

### force computation

The bodies are assigned equally to all MPI processes. Then each process computes the forces for its assigned bodies 
using the full octree obtained by merging all partial octrees. This merging is done on-the-fly by selecting the appropriate partial octree,
once the tree recursion depth hits the root cell level.

### position updates

After the force computation, each MPI process then updates the positions of the bodies, for which it has computed the forces.

Following this each MPI process broadcasts the updated positions to all other MPI processes, 
so they are available in the next octree construction phase on all other nodes.
