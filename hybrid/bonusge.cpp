//
// Created by shubham on 08.06.21.
//

#include "bonusge.h"
#include <stdexcept>
#include <cmath>
#include <list>
#include <sched.h>
#include <unistd.h>
#include <immintrin.h>
#include "omp.h"
#include "mpi.h"
#include "utility.h"
#include "constant.h"

namespace Hybrid {
        ///// Implementation of helper function/class for solve here

OurOctreeNode::OurOctreeNode(uint32_t nodeOrd, uint32_t curLevelOrd, uint32_t parentId, Utility::OctreeCell cell) : nodeOrd_(nodeOrd), curLevelOrd_(curLevelOrd), parentId_(parentId), cell_(std::move(cell))
{
}

std::pair<int, Utility::OctreeCell> OurOctreeNode::selectCell(OurOctree& octree, NodeId ownId, const Body& b)
{
    auto& self = octree.getNode(ownId);
    bool vertUD = self.cell_.IsUp(const_cast<Body&>(b));
    bool longNS = self.cell_.IsNorth(const_cast<Body&>(b));
    bool latEW = self.cell_.IsEast(const_cast<Body&>(b));
    int pos = latEW * (int)(1<<0) + longNS * (int)(1<<1) + vertUD * (int)(1<<2);
    Utility::OctreeCell c;
    switch (pos) {
        case UNE: //Up North East
            c = self.cell_.UpNorthEast();
            break;
        case UNW: //Up North West
            c = self.cell_.UpNorthWest();
            break;
        case USE: //Up South East
            c = self.cell_.UpSouthEast();
            break;
        case USW: //Up South West
            c = self.cell_.UpSouthWest();
            break;
        case DNE: //Down North East
            c = self.cell_.DownNorthEast();
            break;
        case DNW: //Down North West
            c = self.cell_.DownNorthWest();
            break;
        case DSE: //Down South East
            c = self.cell_.DownSouthEast();
            break;
        default: //Down South West (case 0)
            c = self.cell_.DownSouthWest();
            break;
    }
    return std::make_pair(pos, c);
}

void OurOctreeNode::insertBody(OurOctree& octree, const Body& bodyToInsert, NodeId ownId, int32_t curDepth, const RankInfo& rankInfo)
{
    // printf("body(%d).insertBody([%f, %f, %f])\n", ownId, bodyToInsert.pos.x, bodyToInsert.pos.y, bodyToInsert.pos.z);
    if (octree.getNode(ownId).body_.m == 0) {
        octree.getNode(ownId).body_ = bodyToInsert;
    }
    else //if (!IsExternal())
    {
        bool wasLeafNode = octree.getNode(ownId).isLeafNode();
        const Body* updatedBod = nullptr;
        if (!wasLeafNode) {
            auto& self = octree.getNode(ownId);
            self.body_.pos.x = (bodyToInsert.pos.x * bodyToInsert.m + self.body_.pos.x * self.body_.m) / (bodyToInsert.m + self.body_.m);
            self.body_.pos.y = (bodyToInsert.pos.y * bodyToInsert.m + self.body_.pos.y * self.body_.m) / (bodyToInsert.m + self.body_.m);
            self.body_.pos.z = (bodyToInsert.pos.z * bodyToInsert.m + self.body_.pos.z * self.body_.m) / (bodyToInsert.m + self.body_.m);
            self.body_.m += bodyToInsert.m;
            updatedBod = &bodyToInsert;
        }
        auto selection = selectCell(octree, ownId, updatedBod ? *updatedBod : octree.getNode(ownId).body_);
        auto selectedChild = octree.getNode(ownId).children_[selection.first];
        if (!selectedChild) {
            uint32_t childOrder = 0;

            uint32_t parentId = ownId;
            int k = curDepth+1;
            childOrder += ((1<<(k*3))-1)/7;

            int i = 0;
            uint32_t curLevelOrd = selection.first;
            while (parentId != 0) {
                childOrder += (1<<(3*i)) * curLevelOrd;
                curLevelOrd = octree.getNode(parentId).getCurLevelOrd();
                parentId = octree.getNode(parentId).getParentId();
                i++;
            }

            selectedChild = octree.addNode(childOrder, selection.first, ownId, std::move(selection.second));
            octree.getNode(ownId).children_[selection.first] = selectedChild;
        }

        auto nodeOrd = octree.getNode(selectedChild).getNodeOrd();
        if (nodeOrd < rankInfo.firstRootCell || nodeOrd >= rankInfo.endRootCell || (rankInfo.ownRootCellLo <= nodeOrd && rankInfo.ownRootCellHi > nodeOrd))
        {
            octree.getNode(selectedChild).insertBody(octree, updatedBod ? *updatedBod : octree.getNode(ownId).body_, selectedChild, curDepth+1, rankInfo);
        }

        if (wasLeafNode) {
            // printf("self recurse triggered by call to body(%d).insertBody([%f, %f, %f])\n", ownId, bodyToInsert.pos.x, bodyToInsert.pos.y, bodyToInsert.pos.z);
            octree.getNode(ownId).insertBody(octree, bodyToInsert, ownId, curDepth, rankInfo);
            // printf("self recurse end\n");
        }
    }
}

void OurOctreeNode::treeInteract(OurOctree& octree, Body& bod, int curDepth, int ignoreBelowDepth, std::vector<Body*>& interactingBodies)
{
    //////////////////////////////////////////////////////////////////
    ///////// Computes the force on body due to all the others ///////
    /////////////////////////////////////////////////////////////////
    vec3 diff = body_.pos - bod.pos;
    double dist = Utility::magnitude(&diff);
    double theta = cell_.getRadius() / dist;
    if (isLeafNode()) {
        if (curDepth >= ignoreBelowDepth) {
            // Case when a octree cell contains only one element
            //Utility::BodyInteract(&bod, &body_, true);
            // diff = (body_.pos - bod.pos)*TO_METERS;
            // double dist2 = Utility::magnitude(&diff);
            if (dist != 0.0) {
                interactingBodies.push_back(&body_);
            }
        }
    }
    else if (theta < MAX_DISTANCE) {
        if (curDepth >= ignoreBelowDepth) {
            // Case where octree cell contains more than one element but is taken
            // in aggregated fashion
            //Utility::BodyInteract(&bod, &body_, false);
            interactingBodies.push_back(&body_);
        }
    }
    else {
        for (auto c : children_) {
            if (c) {
                octree.getNode(c).treeInteract(octree, bod, curDepth+1, ignoreBelowDepth, interactingBodies);
            }
        }
        // if (getChild(UNW)) octree.getNode(getChild(UNW)).treeInteract(octree, bod, curDepth+1, ignoreBelowDepth);
        // if (getChild(UNE)) octree.getNode(getChild(UNE)).treeInteract(octree, bod, curDepth+1, ignoreBelowDepth);
        // if (getChild(USW)) octree.getNode(getChild(USW)).treeInteract(octree, bod, curDepth+1, ignoreBelowDepth);
        // if (getChild(USE)) octree.getNode(getChild(USE)).treeInteract(octree, bod, curDepth+1, ignoreBelowDepth);
        // if (getChild(DNW)) octree.getNode(getChild(DNW)).treeInteract(octree, bod, curDepth+1, ignoreBelowDepth);
        // if (getChild(DNE)) octree.getNode(getChild(DNE)).treeInteract(octree, bod, curDepth+1, ignoreBelowDepth);
        // if (getChild(DSW)) octree.getNode(getChild(DSW)).treeInteract(octree, bod, curDepth+1, ignoreBelowDepth);
        // if (getChild(DSE)) octree.getNode(getChild(DSE)).treeInteract(octree, bod, curDepth+1, ignoreBelowDepth);
    }
}

size_t OurOctreeNode::size(OurOctree& octree) const
{
    size_t n {1};
    for (auto& c : children_) {
        if (c) {
            n += octree.getNode(c).size(octree);
        }
    }
    return n;
}

void PositionUpdate(Body* bodies, int bodyLo, int bodyHi)
{
    for (int i = bodyLo; i < bodyHi; i++) {
        Utility::Integrate(&bodies[i]);
    }
}

struct ProblemData {
    int rootNodeDepth;

    std::vector<uint32_t> treeSizes;
    std::vector<OurOctree> trees;
};

int GetCappedCommSize()
{
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#ifdef CPU_PINNING
    if (size > 2) {
        size = 2;
    }
#endif
    return size;
}

void BodyInteractSIMD(Body* target, Body** others, size_t numOthers)
{
    auto toMeters = _mm256_set1_pd(TO_METERS);
    auto targetX = _mm256_set1_pd(target->pos.x);
    auto targetY = _mm256_set1_pd(target->pos.y);
    auto targetZ = _mm256_set1_pd(target->pos.z);
    auto G_targetM = _mm256_set1_pd(G * target->m);
    auto SOFTENING_SQ = _mm256_set1_pd(SOFTENING * SOFTENING);
    auto invTargetM = _mm256_set1_pd(1.0 / target->m);

    auto targetAcc = _mm256_set_pd(0.0, target->acc.z, target->acc.y, target->acc.x);

    for (size_t i=0;i<numOthers/4*4;i+=4) {
        //vec3 diff;
        //diff = (target->pos - other->pos)*TO_METERS;
        //double dist = Utility::magnitude(&diff);

        auto& o0 = *others[i+0];
        auto& o1 = *others[i+1];
        auto& o2 = *others[i+2];
        auto& o3 = *others[i+3];

        auto othersM = _mm256_set_pd(o3.m, o2.m, o1.m, o0.m);

        auto othersX = _mm256_set_pd(o3.pos.x, o2.pos.x, o1.pos.x, o0.pos.x);
        auto othersY = _mm256_set_pd(o3.pos.y, o2.pos.y, o1.pos.y, o0.pos.y);
        auto othersZ = _mm256_set_pd(o3.pos.z, o2.pos.z, o1.pos.z, o0.pos.z);

        auto diffX = _mm256_mul_pd(_mm256_sub_pd(targetX, othersX), toMeters);
        auto diffY = _mm256_mul_pd(_mm256_sub_pd(targetY, othersY), toMeters);
        auto diffZ = _mm256_mul_pd(_mm256_sub_pd(targetZ, othersZ), toMeters);

        auto diffX2 = _mm256_mul_pd(diffX, diffX);
        auto diffY2 = _mm256_mul_pd(diffY, diffY);
        auto diffZ2 = _mm256_mul_pd(diffZ, diffZ);

        // diffX2 = d0.x d1.x d2.x d3.x
        // diffY2 = d0.y d1.y d2.y d3.y
        // diffZ2 = d0.z d1.z d2.z d3.z

        auto distSq = _mm256_add_pd(_mm256_add_pd(diffX2, diffY2), diffZ2);
        auto dist = _mm256_sqrt_pd(distSq);
        //auto distReSq = _mm256_mul_pd(dist, dist);

        auto alpha = _mm256_mul_pd(G_targetM, othersM);
        auto F = _mm256_div_pd(alpha, (_mm256_mul_pd(_mm256_add_pd(distSq, SOFTENING_SQ), dist)));

        auto diffX_F_invM = _mm256_mul_pd(_mm256_mul_pd(diffX, F), invTargetM);
        auto diffY_F_invM = _mm256_mul_pd(_mm256_mul_pd(diffY, F), invTargetM);
        auto diffZ_F_invM = _mm256_mul_pd(_mm256_mul_pd(diffZ, F), invTargetM);

        targetAcc = _mm256_sub_pd(targetAcc, _mm256_set_pd(0.0, diffZ_F_invM[0], diffY_F_invM[0], diffX_F_invM[0]) );
        targetAcc = _mm256_sub_pd(targetAcc, _mm256_set_pd(0.0, diffZ_F_invM[1], diffY_F_invM[1], diffX_F_invM[1]) );
        targetAcc = _mm256_sub_pd(targetAcc, _mm256_set_pd(0.0, diffZ_F_invM[2], diffY_F_invM[2], diffX_F_invM[2]) );
        targetAcc = _mm256_sub_pd(targetAcc, _mm256_set_pd(0.0, diffZ_F_invM[3], diffY_F_invM[3], diffX_F_invM[3]) );
    }

    target->acc.x = targetAcc[0];
    target->acc.y = targetAcc[1];
    target->acc.z = targetAcc[2];

    for (size_t i=numOthers/4*4;i<numOthers;++i) {
        Utility::BodyInteract(target, others[i], false);
    }
}

void UpdateStep(Body* bodies, const int numBodies, ProblemData& data, const RankInfo& rankInfo)
{
    /////////////////////////////////////////////////////////////
    /// The function sets up tree and update the ////////////////
    /// positions of bodies for one time step ///////////////////
    ////////////////////////////////////////////////////////////

    // Interaction with sun is individual
    Body* sun = &bodies[0];
    for(int bcount=1; bcount<numBodies; bcount++){
        Utility::BodyInteract(sun, &bodies[bcount]);
    }

    Utility::OctreeCell root = Utility::OctreeCell(0, //// x center of root
                                                    0,       //// y center of root
                                                    0.1374, //// z center of root
                                                    60*SYSTEM_SIZE
    );

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    size = GetCappedCommSize();

    OurOctree& tree = data.trees[rank];
    tree.clear();
    tree.reserve(numBodies * 2);
    auto rootNodeId = tree.addNode(0, 0, 0, std::move(root));

    // Now compute the interaction due to foces between the bodies
    // This makes sure that the bodies that are just too far are
    // not included in the computation makes job easier

    for (int bcount=1; bcount<numBodies; bcount++){
        // Check if the body lies in the system
        if(tree.getNode(rootNodeId).getOctreeCell().contains(bodies[bcount].pos)){
            tree.getNode(rootNodeId).insertBody(tree, bodies[bcount], rootNodeId, 0, rankInfo);
        }
    }

    std::vector<MPI_Request> sendRequests;
    sendRequests.resize(size, MPI_REQUEST_NULL);
    uint32_t localTreeSize = tree.nodes.size();
    //like a non-blocking broadcast that also provides a request object for later probing
    for (int i=0;i<size;++i) {
        if (i != rank) {
            MPI_Isend(tree.nodes.data(), localTreeSize * sizeof(OurOctreeNode), MPI_BYTE, i, rank, MPI_COMM_WORLD, &sendRequests[i]);
        }
    }

    std::vector<bool> transfersPending;
    int numTransfersPending { size - 1 };
    transfersPending.resize(size, true);
    transfersPending[rank] = false;

    std::vector<MPI_Request> recvRequests;
    recvRequests.resize(size, MPI_REQUEST_NULL);

    //like a non-blocking allgather where order does not matter, thats why its faster
    //also the treeSize can be determined in the same step. No extra broadcast needed.
    while (numTransfersPending > 0) {
      MPI_Status probeStatus;
      int ready;
      MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &ready, &probeStatus);

      if (ready && transfersPending[probeStatus.MPI_SOURCE]) {
        int i = probeStatus.MPI_SOURCE;
        int treeDataSize;
        MPI_Get_count(&probeStatus, MPI_BYTE, &treeDataSize);
        auto treeSize = treeDataSize / sizeof(OurOctreeNode);
        data.trees[i].nodes.resize(treeSize);

        MPI_Irecv(data.trees[i].nodes.data(), treeDataSize, MPI_BYTE, i, i, MPI_COMM_WORLD, &recvRequests[i]);

        transfersPending[i] = false;
        numTransfersPending--;
      }
    }

    MPI_Waitall(recvRequests.size(), recvRequests.data(), MPI_STATUSES_IGNORE);
    MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUSES_IGNORE);

    thread_local std::vector<Body*> interactingBodies;

    #pragma omp parallel for schedule(guided) num_threads(24)
    for (int bcount = rankInfo.bodyLo; bcount < rankInfo.bodyHi; bcount++) {
        interactingBodies.clear();
        // Check if the body lies in the system
        for (int i=0;i<size;++i) {
            auto& curTree = data.trees[i];
            int ignoreBelowDepth = i == rank ? 0 : data.rootNodeDepth;
            if (curTree.getNode(rootNodeId).getOctreeCell().contains(bodies[bcount].pos)) {
                curTree.getNode(rootNodeId).treeInteract(curTree, bodies[bcount], 0, ignoreBelowDepth, interactingBodies);
            }
        }
        if (!interactingBodies.empty()) {
            //printf("interactingBodies count: %d\n", interactingBodies.size());
            // if (interactingBodies.size() > 1) {
            //     throw std::logic_error { "logic error" };
            // }
            BodyInteractSIMD(&bodies[bcount], interactingBodies.data(), interactingBodies.size());
        }
    }


    /*
    for (int i=0;i<size;++i) {
            // Check if the body lies in the system
            for (int bcount = rankInfo.bodyLo; bcount < rankInfo.bodyHi; bcount++) {
                auto& curTree = data.trees[i];
                int ignoreBelowDepth = i == rank ? 0 : data.rootNodeDepth;
                if (curTree.getNode(rootNodeId).getOctreeCell().contains(bodies[bcount].pos)) {
                    curTree.getNode(rootNodeId).treeInteract(curTree, bodies[bcount], 0, ignoreBelowDepth);
                }
            }
        }
    */

    PositionUpdate(bodies, rankInfo.bodyLo, rankInfo.bodyHi);

    std::vector<int> countsRecv;
    std::vector<int> displacementsRecv;
    countsRecv.resize(size);
    displacementsRecv.resize(size);

    int bodiesPerNode = 1 + ((numBodies - 2) / size);
    for (int i=0;i<size;++i) {
        int nodeBodyLo = 1 + bodiesPerNode * i;
        int nodeBodyHi = 1 + std::min(bodiesPerNode * (i+1), numBodies-1);
        int nodeBodyCount = nodeBodyHi - nodeBodyLo;
        int chunkSize = nodeBodyCount * sizeof(Body);
        countsRecv[i] = chunkSize;
        displacementsRecv[i] = nodeBodyLo * sizeof(Body);
    }

    MPI_Request request;

    MPI_Iallgatherv(MPI_IN_PLACE, countsRecv[rank], MPI_BYTE,
        bodies, countsRecv.data(), displacementsRecv.data(), MPI_BYTE, MPI_COMM_WORLD, &request);
    MPI_Wait(&request, MPI_STATUS_IGNORE);
    /*
    MPI_Allgatherv(MPI_IN_PLACE, countsRecv[rank], MPI_BYTE,
        bodies, countsRecv.data(), displacementsRecv.data(), MPI_BYTE, MPI_COMM_WORLD);
    */
}

std::pair<std::vector<std::pair<int,int>>,int> CalculateRootCells(int size)
{
    int rootCellOrderOffset = 1;
    int rootCellLevel = 1;
    int rootCellsTotal = 8;
    while (size > rootCellsTotal) {
        rootCellLevel++;
        rootCellOrderOffset += rootCellsTotal;
        rootCellsTotal *= 8;
    }

    int cellsPerNode = 1 + ((rootCellsTotal - 1) / size);
    std::vector<std::pair<int,int>> rootCells;
    for (int nodeId=0;nodeId<size;++nodeId) {
        int lowerCell = rootCellOrderOffset + cellsPerNode * nodeId;
        int upperCell = rootCellOrderOffset + std::min(cellsPerNode * (nodeId+1), rootCellsTotal);
        //printf("rank %d: [%d,%d]\n", nodeId, lowerCell, upperCell);
        rootCells.emplace_back(lowerCell, upperCell);
    }

    return std::make_pair(rootCells, rootCellLevel);
}

}//end namespace
void Hybrid::Solve(Body *bodies, int numBodies, int numSteps) {
#ifdef VISUALIZATION
    char* image = new char[IMAGEWIDTH * IMAGEHEIGHT * 3];
    double* imageArray = new double[IMAGEWIDTH * IMAGEHEIGHT * 3];
    Visualizer::SaveFrame(image, imageArray, bodies, numBodies, 1);
#endif
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    size = GetCappedCommSize();

#ifdef CPU_PINNING
    cpu_set_t mask;
    CPU_ZERO(&mask);

    if (rank == 0) {
        for (int i=0;i<=7;++i) {
            CPU_SET(i, &mask);
        }
        for (int i=16;i<=23;++i) {
            CPU_SET(i, &mask);
        }
    }
    else if (rank == 1) {
        for (int i=8;i<=15;++i) {
            CPU_SET(i, &mask);
        }
        for (int i=24;i<=31;++i) {
            CPU_SET(i, &mask);
        }
    }
    else {
        return;
    }

    sched_setaffinity(0, sizeof(mask), &mask);
    omp_set_num_threads(16);
#endif

    MPI_Bcast(&numBodies, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&numSteps, 1, MPI_INT, 0, MPI_COMM_WORLD);
    std::vector<Body> localBodies;
    if (bodies == nullptr) {
        localBodies.resize(numBodies);
        bodies = localBodies.data();
    }
    MPI_Bcast(bodies, numBodies * sizeof(Body), MPI_BYTE, 0, MPI_COMM_WORLD);

    RankInfo rankInfo {};

    auto rootCellResult = CalculateRootCells(size);
    auto cellRangesByRank = rootCellResult.first;

    rankInfo.firstRootCell = cellRangesByRank[0].first;
    rankInfo.endRootCell = cellRangesByRank.back().second;
    auto cellRange = cellRangesByRank[rank];
    rankInfo.ownRootCellLo = cellRange.first;
    rankInfo.ownRootCellHi = cellRange.second;

    int bodiesPerNode = 1 + ((numBodies - 2) / size);
    rankInfo.bodyLo = 1 + bodiesPerNode * rank;
    rankInfo.bodyHi = 1 + std::min(bodiesPerNode * (rank+1), numBodies-1);

    ProblemData problemData;
    problemData.rootNodeDepth = rootCellResult.second;
    problemData.treeSizes.resize(size);
    problemData.trees.resize(size);
    for (int step = 1; step < numSteps; step++) {
        Hybrid::UpdateStep(bodies, numBodies, problemData, rankInfo);
#ifdef VISUALIZATION
        if (step % RENDER_INTERVAL == 0) {
            Visualizer::SaveFrame(image, imageArray, bodies, numBodies, step + 1);
        }
#endif
    }

#ifdef VISUALIZATION
    delete[] image;
    delete[] imageArray;
#endif
}

