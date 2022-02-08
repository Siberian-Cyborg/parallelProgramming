//
// Created by shubham on 08.06.21.
//

#include <cassert>
#include <iostream>
#include "constant.h"
#include <utility>
#include <cmath>
#include <chrono>
#include <vector>
#include "utility.h"
#include "serial.h"
#include "visualizer.h"

const Utility::OctreeCell& Serial::Octree::getOctreeCell () const {
    return thisNode;
}

Serial::Octree::~Octree(){
    if (upNorthWest!= nullptr) delete upNorthWest;
    if (upNorthEast!= nullptr) delete upNorthEast;
    if (upSouthWest!= nullptr) delete upSouthWest;
    if (upSouthEast!= nullptr) delete upSouthEast;
    if (downNorthWest!= nullptr) delete downNorthWest;
    if (downNorthEast!= nullptr) delete downNorthEast;
    if (downSouthWest!= nullptr) delete downSouthWest;
    if (downSouthEast!= nullptr) delete downSouthEast;
}
bool Serial::Octree::IsExternal()
{
    return (upNorthWest== nullptr) && (upNorthEast== nullptr) && (upSouthWest== nullptr) &&
           (upSouthEast== nullptr) && (downNorthWest== nullptr) && (downNorthEast== nullptr) &&
           (downSouthWest== nullptr) && (downSouthEast== nullptr);
}
void Serial::Octree::InsertBody(Body* insertBod)
{
    if (b.m == 0)
    {
        b = *insertBod;
    } else //if (!IsExternal())
    {
        bool isExtern = IsExternal();
        Body* updatedBod;
        if (!isExtern)
        {
            b.pos.x = (insertBod->pos.x * insertBod->m + b.pos.x * b.m)/(insertBod->m + b.m);
            b.pos.y = (insertBod->pos.y * insertBod->m + b.pos.y * b.m)/(insertBod->m + b.m);
            b.pos.z = (insertBod->pos.z * insertBod->m + b.pos.z * b.m)/(insertBod->m + b.m);
            b.m += insertBod->m;
            updatedBod = insertBod;
        } else{
            updatedBod = &b;
        }
        Utility::OctreeCell&& casecell = thisNode.DownNorthEast();
        Serial::Octree** casenode;
        CaseSelector(updatedBod, casecell, &casenode);
        if (casecell.contains(updatedBod->pos)){
            if ((*casenode)== nullptr){
                (*casenode) = new Serial::Octree(std::move(casecell));
            }
            (*casenode)->InsertBody(updatedBod);
        }
        if (isExtern) {
            InsertBody(insertBod);
        }
    }
}


void Serial::Octree::CaseSelector(Body* b, Utility::OctreeCell& c, Serial::Octree*** casenode){
    /////////////////////////////////////////////////////////////////////////////////////
    //////// Uses 3 digit binary representation to find the correct octant //////////////
    /////////////////////////////////////////////////////////////////////////////////////

    bool vertUD = thisNode.IsUp(*b);
    bool longNS = thisNode.IsNorth(*b);
    bool latEW = thisNode.IsEast(*b);
    constexpr int UNE = 1 + 2 + 4;
    constexpr int UNW = 2 + 4;
    constexpr int USE = 1 + 4;
    constexpr int USW = 4;
    constexpr int DNE = 1 + 2;
    constexpr int DNW =  2;
    constexpr int DSE = 1;
//    constexpr int DSW = 0;
    int pos = latEW*(int)pow(2,0) + longNS*(int)pow(2,1) + vertUD*(int)pow(2,2);
    switch (pos) {
        case UNE: //Up North East
            c = thisNode.UpNorthEast();
            *casenode = &upNorthEast;
            break;

        case UNW: //Up North West
            c = thisNode.UpNorthWest();
            *casenode = &upNorthWest;
            break;
        case USE: //Up South East
            c = thisNode.UpSouthEast();
            *casenode = &upSouthEast;
            break;
        case USW: //Up South West
            c = thisNode.UpSouthWest();
            *casenode = &upSouthWest;
            break;
        case DNE: //Down North East
            c = thisNode.DownNorthEast();
            *casenode = &downNorthEast;
            break;
        case DNW: //Down North West
            c = thisNode.DownNorthWest();
            *casenode = &downNorthWest;
            break;
        case DSE: //Down South East
            c = thisNode.DownSouthEast();
            *casenode = &downSouthEast;
            break;
        default: //Down South West (case 0)
            c = thisNode.DownSouthWest();
            *casenode = &downSouthWest;
            break;
    }
}

void Serial::Octree::TreeInteract(Body* bod)
{
    //////////////////////////////////////////////////////////////////
    ///////// Computes the force on body due to all the others ///////
    /////////////////////////////////////////////////////////////////
    vec3 diff = b.pos - bod->pos;
    double dist = Utility::magnitude(&diff);
    double theta = thisNode.getRadius() / dist;
    if (IsExternal())
    {
        // Case when a octree cell contains only one element
        Utility::BodyInteract(bod, &b, true);
    }
    else if (theta < MAX_DISTANCE)
    {
        // Case where octree cell contains more than one element but is taken
        // in aggregated fashion
        Utility::BodyInteract(bod, &b, false);
    } else {
        // Case where further refinement is necessary
        if (upNorthWest!=NULL) upNorthWest->TreeInteract(bod);
        if (upNorthEast!=NULL) upNorthEast->TreeInteract(bod);
        if (upSouthWest!=NULL) upSouthWest->TreeInteract(bod);
        if (upSouthEast!=NULL) upSouthEast->TreeInteract(bod);
        if (downNorthWest!=NULL) downNorthWest->TreeInteract(bod);
        if (downNorthEast!=NULL) downNorthEast->TreeInteract(bod);
        if (downSouthWest!=NULL) downSouthWest->TreeInteract(bod);
        if (downSouthEast!=NULL) downSouthEast->TreeInteract(bod);
    }
}

namespace Serial{
    void PositionUpdate(Body* bodies, const int numBodies){
        for (int i=0; i<numBodies; i++){
            Utility::Integrate(&bodies[i]);
        }
    }
    void UpdateStep(Body* bodies, const int numBodies){
        /////////////////////////////////////////////////////////////
        /// The function sets up tree and update the ////////////////
        /// positions of bodies for one time step ///////////////////
        ////////////////////////////////////////////////////////////

        // Interaction with sun is individual
        Body* sun = &bodies[0];
        for(int bcount=1; bcount<numBodies; bcount++){
            Utility::BodyInteract(sun, &bodies[bcount]);
        }

        Utility::OctreeCell&& root = Utility::OctreeCell(0, //// x center of root
                0,       //// y center of root
                0.1374, //// z center of root
                60*SYSTEM_SIZE
                );

        Serial::Octree* tree = new Serial::Octree(std::move(root));

        // Now compute the interaction due to foces between the bodies
        // This makes sure that the bodies that are just too far are
        // not included in the computation makes job easier
        for (int bcount=1; bcount<numBodies; bcount++){
            // Check if the body lies in the system
            if(tree->getOctreeCell().contains(bodies[bcount].pos)){
                tree->InsertBody(&bodies[bcount]);
            }
        }

        for (int bcount=1; bcount<numBodies; bcount++){
            // Check if the body lies in the system
            if(tree->getOctreeCell().contains(bodies[bcount].pos)){
                tree->TreeInteract(&bodies[bcount]);
            }
        }

        // remove the tree
        delete tree;
        Serial::PositionUpdate(bodies, numBodies);
    }
}

void Serial::SerialSolve(int argc, char **argv, float &sequential_runtime, int rank){
    if (rank==0){
        int numBodies, numSteps;
        Utility::ParseInputs(argc, argv, &numBodies, &numSteps);
        Body *bodies = new Body[numBodies];
        Utility::Initialize(bodies, numBodies, 0);
        auto start = std::chrono::steady_clock::now();
        Serial::Solve(bodies, numBodies, numSteps);
        auto end = std::chrono::steady_clock::now();
        sequential_runtime = std::chrono::duration<float>(end - start).count();

        Utility::PrintBodies(bodies, numBodies);
        delete[] bodies;
    }
}

void Serial::Solve(Body* bodies, const int numBodies, const int numSteps){
#ifdef VISUALIZATION
    char *image = new char[IMAGEWIDTH*IMAGEHEIGHT*3];
    double *imageArray = new double[IMAGEWIDTH*IMAGEHEIGHT*3];
    Visualizer::SaveFrame(image, imageArray, bodies, numBodies, 1);
#endif
    for(int step=1; step<numSteps; step++){
        Serial::UpdateStep(bodies, numBodies);
#ifdef VISUALIZATION
        if(step%RENDER_INTERVAL==0){
            Visualizer::SaveFrame(image, imageArray, bodies, numBodies, step+1);
        }
#endif

    }

#ifdef VISUALIZATION
    delete[] image;
    delete[] imageArray;
#endif

}
