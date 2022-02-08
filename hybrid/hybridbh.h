//
// Created by shubham on 08.06.21.
//

#ifndef TEST_HYBRIDBH_H
#define TEST_HYBRIDBH_H

#include "constant.h"
#include <algorithm>
#include <array>
#include <memory>
#include <vector>
#include "utility.h"

namespace Hybrid{
    
using NodeId = uint32_t;

class OurOctree;

struct RankInfo {
    uint32_t firstRootCell;
    uint32_t endRootCell;
    uint32_t ownRootCellLo;
    uint32_t ownRootCellHi;

    uint32_t bodyLo;
    uint32_t bodyHi;
};

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

    OurOctreeNode() = default;

    explicit OurOctreeNode(uint32_t nodeOrd, uint32_t curLevelOrd, uint32_t parentId, Utility::OctreeCell cell);

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

    static void insertBody(OurOctree& octree, const Body& bodyToInsert,
        NodeId ownId, int32_t curDepth,
        const RankInfo& rankInfo);

    void treeInteract(OurOctree& octree, Body& bod, int curDepth, int ignoreBelowDepth);

    bool isLeafNode() const
    {
        //return std::all_of(children_.begin(), children_.end(), [] (auto& c) { return !c; } );
        return !children_[0] && !children_[1] && !children_[2] && !children_[3] && !children_[4]
            && !children_[5] && !children_[6] && !children_[7];
    }

    size_t size(OurOctree& octree) const;

    uint32_t getParentId() const
    {
        return parentId_;
    }

    uint32_t getCurLevelOrd() const
    {
        return curLevelOrd_;
    }

    uint32_t getNodeOrd() const
    {
        return nodeOrd_;
    }

private:
    Body body_ {};
    Utility::OctreeCell cell_ {};
    std::array<NodeId, 8> children_ {};
    uint32_t nodeOrd_ {};
    uint32_t curLevelOrd_ {};
    uint32_t parentId_ {};
};

class OurOctree {
public:
    OurOctreeNode& getNode(NodeId id)
    {
        return nodes[id - 1];
    }

    const OurOctreeNode& getNode(NodeId id) const
    {
        return nodes[id - 1];
    }

    NodeId addNode(uint32_t nodeOrd, uint32_t curLevelOrd, uint32_t parentId, Utility::OctreeCell cell)
    {
        nodes.emplace_back(nodeOrd, curLevelOrd, parentId, std::move(cell));
        return nodes.size();
    }

    void reserve(size_t numNodes)
    {
        nodes.reserve(numNodes);
    }

    void clear()
    {
        nodes.clear();
    }

    size_t getNodeCount() const
    {
        return nodes.size();
    }

public:
    std::vector<OurOctreeNode> nodes;
};

std::pair<std::vector<std::pair<int,int>>,int> CalculateRootCells(int size);

void Solve(Body* bodies, int numBodies, int numSteps);
}

#endif //TEST_HYBRIDBH_H
