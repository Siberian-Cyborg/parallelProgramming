//
// Created by shubham on 08.06.21.
//

#include <stdexcept>
#include <cmath>
#include <omp.h>
#include "utility.h"
#include "ompbh.h"
#include "constant.h"

namespace OMP {

const Utility::OctreeCell& Octree::getOctreeCell() const
{
    return thisNode;
}

Octree::~Octree()
{
    if (upNorthWest != nullptr) delete upNorthWest;
    if (upNorthEast != nullptr) delete upNorthEast;
    if (upSouthWest != nullptr) delete upSouthWest;
    if (upSouthEast != nullptr) delete upSouthEast;
    if (downNorthWest != nullptr) delete downNorthWest;
    if (downNorthEast != nullptr) delete downNorthEast;
    if (downSouthWest != nullptr) delete downSouthWest;
    if (downSouthEast != nullptr) delete downSouthEast;
}
bool Octree::IsExternal()
{
    return (upNorthWest == nullptr) && (upNorthEast == nullptr) && (upSouthWest == nullptr) && (upSouthEast == nullptr) && (downNorthWest == nullptr) && (downNorthEast == nullptr) && (downSouthWest == nullptr) && (downSouthEast == nullptr);
}

void Octree::InsertBody(Body* insertBod)
{
    if (b.m == 0) {
        b = *insertBod;
    }
    else //if (!IsExternal())
    {
        bool isExtern = IsExternal();
        Body* updatedBod;
        if (!isExtern) {
            b.pos.x = (insertBod->pos.x * insertBod->m + b.pos.x * b.m) / (insertBod->m + b.m);
            b.pos.y = (insertBod->pos.y * insertBod->m + b.pos.y * b.m) / (insertBod->m + b.m);
            b.pos.z = (insertBod->pos.z * insertBod->m + b.pos.z * b.m) / (insertBod->m + b.m);
            b.m += insertBod->m;
            updatedBod = insertBod;
        }
        else {
            updatedBod = &b;
        }
        Utility::OctreeCell&& casecell = thisNode.DownNorthEast();
        Octree** casenode;
        CaseSelector(updatedBod, casecell, &casenode);
        //if (casecell.contains(updatedBod->pos)) {
            if ((*casenode) == nullptr) {
                (*casenode) = new Octree(std::move(casecell));
            }
            (*casenode)->InsertBody(updatedBod);
        //}
        //else {
        //    throw std::runtime_error { "this should never happen" };
        //}
        if (isExtern) {
            InsertBody(insertBod);
        }
    }
}

void Octree::CaseSelector(Body* b, Utility::OctreeCell& c, Octree*** casenode)
{
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
    constexpr int DNW = 2;
    constexpr int DSE = 1;
    //    constexpr int DSW = 0;
    //
    int pos = latEW * (int)(1<<0) + longNS * (int)(1<<1) + vertUD * (int)(1<<2);
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

void Octree::TreeInteract(Body* bod)
{
    //////////////////////////////////////////////////////////////////
    ///////// Computes the force on body due to all the others ///////
    /////////////////////////////////////////////////////////////////
    vec3 diff = b.pos - bod->pos;
    double dist = Utility::magnitude(&diff);
    double theta = thisNode.getRadius() / dist;
    if (IsExternal()) {
        // Case when a octree cell contains only one element
        Utility::BodyInteract(bod, &b, true);
    }
    else if (theta < MAX_DISTANCE) {
        // Case where octree cell contains more than one element but is taken
        // in aggregated fashion
        Utility::BodyInteract(bod, &b, false);
    }
    else {
        // Case where further refinement is necessary
        if (upNorthWest != NULL) upNorthWest->TreeInteract(bod);
        if (upNorthEast != NULL) upNorthEast->TreeInteract(bod);
        if (upSouthWest != NULL) upSouthWest->TreeInteract(bod);
        if (upSouthEast != NULL) upSouthEast->TreeInteract(bod);
        if (downNorthWest != NULL) downNorthWest->TreeInteract(bod);
        if (downNorthEast != NULL) downNorthEast->TreeInteract(bod);
        if (downSouthWest != NULL) downSouthWest->TreeInteract(bod);
        if (downSouthEast != NULL) downSouthEast->TreeInteract(bod);
    }
}

OurOctreeNode::OurOctreeNode(Utility::OctreeCell cell) : cell_(std::move(cell))
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

void OurOctreeNode::insertBody(OurOctree& octree, const Body& bodyToInsert, NodeId ownId)
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
            selectedChild = octree.addNode(std::move(selection.second));
            octree.getNode(ownId).children_[selection.first] = selectedChild;
        }
        octree.getNode(selectedChild).insertBody(octree, updatedBod ? *updatedBod : octree.getNode(ownId).body_, selectedChild);
        if (wasLeafNode) {
            // printf("self recurse triggered by call to body(%d).insertBody([%f, %f, %f])\n", ownId, bodyToInsert.pos.x, bodyToInsert.pos.y, bodyToInsert.pos.z);
            octree.getNode(ownId).insertBody(octree, bodyToInsert, ownId);
            // printf("self recurse end\n");
        }
    }
}

void OurOctreeNode::treeInteract(OurOctree& octree, Body& bod)
{
    //////////////////////////////////////////////////////////////////
    ///////// Computes the force on body due to all the others ///////
    /////////////////////////////////////////////////////////////////
    vec3 diff = body_.pos - bod.pos;
    double dist = Utility::magnitude(&diff);
    double theta = cell_.getRadius() / dist;
    if (isLeafNode()) {
        // Case when a octree cell contains only one element
        Utility::BodyInteract(&bod, &body_, true);
    }
    else if (theta < MAX_DISTANCE) {
        // Case where octree cell contains more than one element but is taken
        // in aggregated fashion
        Utility::BodyInteract(&bod, &body_, false);
    }
    else {
        for (auto c : children_) {
            if (c) {
                octree.getNode(c).treeInteract(octree, bod);
            }
        }
        // if (getChild(UNW)) octree.getNode(getChild(UNW)).treeInteract(octree, bod);
        // if (getChild(UNE)) octree.getNode(getChild(UNE)).treeInteract(octree, bod);
        // if (getChild(USW)) octree.getNode(getChild(USW)).treeInteract(octree, bod);
        // if (getChild(USE)) octree.getNode(getChild(USE)).treeInteract(octree, bod);
        // if (getChild(DNW)) octree.getNode(getChild(DNW)).treeInteract(octree, bod);
        // if (getChild(DNE)) octree.getNode(getChild(DNE)).treeInteract(octree, bod);
        // if (getChild(DSW)) octree.getNode(getChild(DSW)).treeInteract(octree, bod);
        // if (getChild(DSE)) octree.getNode(getChild(DSE)).treeInteract(octree, bod);
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

void PositionUpdate(Body* bodies, const int numBodies)
{
    for (int i = 0; i < numBodies; i++) {
        Utility::Integrate(&bodies[i]);
    }
}

void UpdateStep(Body* bodies, const int numBodies)
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

    Utility::OctreeCell&& root = Utility::OctreeCell(0, //// x center of root
                                                     0,       //// y center of root
                                                     0.1374, //// z center of root
                                                     60*SYSTEM_SIZE
    );

    Octree* tree = new Octree(std::move(root));

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
    PositionUpdate(bodies, numBodies);
}

void UpdateStepOurs(Body* bodies, const int numBodies, OurOctree& tree)
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

    tree.clear();
    tree.reserve(numBodies * 2);
    auto rootNodeId = tree.addNode(std::move(root));

    // Now compute the interaction due to foces between the bodies
    // This makes sure that the bodies that are just too far are
    // not included in the computation makes job easier
    for (int bcount=1; bcount<numBodies; bcount++){
        // Check if the body lies in the system
        if(tree.getNode(rootNodeId).getOctreeCell().contains(bodies[bcount].pos)){
            tree.getNode(rootNodeId).insertBody(tree, bodies[bcount], rootNodeId);
        }
    }

    #pragma omp parallel for schedule(guided) num_threads(24)
    for (int bcount = 1; bcount < numBodies; bcount++) {
        // Check if the body lies in the system
        if (tree.getNode(rootNodeId).getOctreeCell().contains(bodies[bcount].pos)) {
            tree.getNode(rootNodeId).treeInteract(tree, bodies[bcount]);
        }
    }

    // remove the tree
    //tree.reset();
    PositionUpdate(bodies, numBodies);
}

} // namespace OMP

void OMP::Solve(Body* bodies, const int numBodies, const int numSteps)
{
#ifdef VISUALIZATION
    char* image = new char[IMAGEWIDTH * IMAGEHEIGHT * 3];
    double* imageArray = new double[IMAGEWIDTH * IMAGEHEIGHT * 3];
    Visualizer::SaveFrame(image, imageArray, bodies, numBodies, 1);
#endif
    OurOctree tree;
    for (int step = 1; step < numSteps; step++) {
        OMP::UpdateStepOurs(bodies, numBodies, tree);
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
