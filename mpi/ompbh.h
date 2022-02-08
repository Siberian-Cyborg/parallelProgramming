//
// Created by shubham on 08.06.21.
//

#ifndef TEST_OMPBH_H
#define TEST_OMPBH_H

#include <algorithm>
#include <array>
#include <memory>
#include <vector>
#include "constant.h"

namespace OMP {
class Octree {
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
    Octree(Utility::OctreeCell&& o) : thisNode(std::move(o))
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
    Octree(const Utility::OctreeCell& o) : thisNode(o)
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
    const Utility::OctreeCell& getOctreeCell() const;
    ~Octree();
    bool IsExternal();
    void InsertBody(Body* insertBod);
    void TreeInteract(Body* bod);

    const Body& getBody() const
    {
        return b;
    }

    Octree* getChild(int pos) const
    {
        constexpr int UNE = 1 + 2 + 4;
        constexpr int UNW = 2 + 4;
        constexpr int USE = 1 + 4;
        constexpr int USW = 4;
        constexpr int DNE = 1 + 2;
        constexpr int DNW = 2;
        constexpr int DSE = 1;
        constexpr int DSW = 0;
        switch (pos) {
            case UNE: //Up North East
                return upNorthEast;
            case UNW: //Up North West
                return upNorthWest;
            case USE: //Up South East
                return upSouthEast;
            case USW: //Up South West
                return upSouthWest;
            case DNE: //Down North East
                return downNorthEast;
            case DNW: //Down North West
                return downNorthWest;
            case DSE: //Down South East
                return downSouthEast;
            case DSW: //Down South West (case 0)
                return downSouthWest;
            default:
                return nullptr;
        }
    }
};

using NodeId = uint32_t;

class OurOctree;

class OurOctreeNode {
public:
    enum {
        UNE = 1 + 2 + 4,
        UNW = 2 + 4,
        USE = 1 + 4,
        USW = 4,
        DNE = 1 + 2,
        DNW = 2,
        DSE = 1,
        DSW = 0
    };

    explicit OurOctreeNode(Utility::OctreeCell cell);

    const Utility::OctreeCell& getOctreeCell() const
    {
        return cell_;
    }

    const Body& getBody() const
    {
        return body_;
    }

    NodeId getChild(int i) const
    {
        return children_[i];
    }

    static std::pair<int, Utility::OctreeCell> selectCell(OurOctree& octree, NodeId ownId, const Body& b);

    static void insertBody(OurOctree& octree, const Body& bodyToInsert, NodeId ownId);

    void treeInteract(OurOctree& octree, Body& bod);

    bool isLeafNode() const
    {
        //return std::all_of(children_.begin(), children_.end(), [] (auto& c) { return !c; } );
        return !children_[0] && !children_[1] && !children_[2] && !children_[3] && !children_[4]
            && !children_[5] && !children_[6] && !children_[7];
    }

    size_t size(OurOctree& octree) const;

private:
    Body body_ {};
    Utility::OctreeCell cell_ {};
    std::array<NodeId, 8> children_ {};
};

class OurOctree {
public:
    OurOctreeNode& getNode(NodeId id)
    {
        return nodes_[id - 1];
    }

    const OurOctreeNode& getNode(NodeId id) const
    {
        return nodes_[id - 1];
    }

    NodeId addNode(Utility::OctreeCell cell)
    {
        nodes_.emplace_back(std::move(cell));
        return nodes_.size();
    }

    void reserve(size_t numNodes)
    {
        nodes_.reserve(numNodes);
    }

    void clear()
    {
        nodes_.clear();
    }

private:
    std::vector<OurOctreeNode> nodes_;
};

void Solve(Body* bodies, const int numBodies, const int numSteps);
} // namespace OMP

#endif //TEST_OMPBH_H
