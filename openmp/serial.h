//
// Created by shubham on 08.06.21.
//

#ifndef BARNES_HUT_SERIAL_H
#define BARNES_HUT_SERIAL_H
#include "utility.h"
namespace Serial{
    class Octree{
    private:
        Body b;
        Utility::OctreeCell thisNode;
        Octree* upNorthWest;
        Octree* upNorthEast;
        Octree* upSouthWest;
        Octree* upSouthEast;
        Octree* downNorthWest;
        Octree* downNorthEast;
        Octree* downSouthWest;
        Octree* downSouthEast;
        void CaseSelector(Body* b, Utility::OctreeCell& casecell, Octree*** casenode);

    public:
        Octree(Utility::OctreeCell&& o): thisNode(std::move(o))
        {
            upNorthWest = nullptr;
            upNorthEast = nullptr;
            upSouthWest = nullptr;
            upSouthEast = nullptr;
            downNorthWest = nullptr;
            downNorthEast = nullptr;
            downSouthWest = nullptr;
            downSouthEast = nullptr;
        }
        Octree(const Utility::OctreeCell& o): thisNode(o)
        {
            upNorthWest = nullptr;
            upNorthEast = nullptr;
            upSouthWest = nullptr;
            upSouthEast = nullptr;
            downNorthWest = nullptr;
            downNorthEast = nullptr;
            downSouthWest = nullptr;
            downSouthEast = nullptr;
        }
        const Utility::OctreeCell& getOctreeCell () const;
        ~Octree();
        bool IsExternal();
        void InsertBody(Body* insertBod);
        void TreeInteract(Body* bod);
    };

    void Solve(Body* bodies, const int numBodies, const int numSteps);
    void SerialSolve(int argc, char **argv, float &sequential_runtime, int rank);

}
#endif //BARNES_HUT_SERIAL_H
