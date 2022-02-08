#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "mpibh.h"

using namespace MPI;
using namespace ::testing;

// NOLINTNEXTLINE
TEST(test_mpi, RootCellCalculationTest)
{
    EXPECT_THAT(CalculateRootCells(4).first,
        ElementsAre(Pair(1, 3), Pair(3, 5), Pair(5, 7), Pair(7, 9)));

    EXPECT_THAT(CalculateRootCells(8).first,
                ElementsAre(Pair(1, 2), Pair(2, 3), Pair(3, 4),
                            Pair(4, 5), Pair(5, 6), Pair(6, 7), Pair(7, 8), Pair(8, 9)));

    auto level2 = CalculateRootCells(9).first;
    EXPECT_THAT(level2,
                ElementsAre(Pair(9, 17), Pair(17, 25), Pair(25, 33), Pair(33, 41), Pair(41, 49), Pair(49, 57), Pair(57, 65), Pair(65, 73), Pair(73, 73)));

    auto level3 = CalculateRootCells(65).first;
    level3.resize(5);
    EXPECT_THAT(level3,
                ElementsAre(Pair(73, 81), Pair(81, 89), Pair(89, 97), Pair(97, 105), Pair(105, 113)));
}

class MpiOctree : public ::testing::Test {
public:
    void compareOctree(OurOctree& leftTree, MPI::NodeId leftNodeId, OurOctree& rightTree, MPI::NodeId rightNodeId, const std::string& path = "R")
    {
        auto left = &leftTree.getNode(leftNodeId);
        auto right = &rightTree.getNode(rightNodeId);
        EXPECT_DOUBLE_EQ(left->getBody().m, right->getBody().m) << path;

        EXPECT_DOUBLE_EQ(left->getBody().pos.x, right->getBody().pos.x) << path;
        EXPECT_DOUBLE_EQ(left->getBody().pos.y, right->getBody().pos.y) << path;
        EXPECT_DOUBLE_EQ(left->getBody().pos.z, right->getBody().pos.z) << path;

        EXPECT_DOUBLE_EQ(left->getBody().vel.x, right->getBody().vel.x) << path;
        EXPECT_DOUBLE_EQ(left->getBody().vel.y, right->getBody().vel.y) << path;
        EXPECT_DOUBLE_EQ(left->getBody().vel.z, right->getBody().vel.z) << path;

        EXPECT_DOUBLE_EQ(left->getBody().acc.x, right->getBody().acc.x) << path;
        EXPECT_DOUBLE_EQ(left->getBody().acc.y, right->getBody().acc.y) << path;
        EXPECT_DOUBLE_EQ(left->getBody().acc.z, right->getBody().acc.z) << path;

        EXPECT_DOUBLE_EQ(left->getOctreeCell().getMidpoint().x, right->getOctreeCell().getMidpoint().x) << path;
        EXPECT_DOUBLE_EQ(left->getOctreeCell().getMidpoint().y, right->getOctreeCell().getMidpoint().y) << path;
        EXPECT_DOUBLE_EQ(left->getOctreeCell().getMidpoint().z, right->getOctreeCell().getMidpoint().z) << path;
        EXPECT_DOUBLE_EQ(left->getOctreeCell().getRadius(), right->getOctreeCell().getRadius()) << path;

        for (int i = 0; i < 8; ++i) {
            auto leftChild = left->getChild(i);
            auto rightChild = right->getChild(i);
            EXPECT_EQ(!!leftChild, !!rightChild) << path;
            if (leftChild && rightChild) {
                compareOctree(leftTree, leftChild, rightTree, rightChild, path + "." + std::to_string(i));
            }
        }
    }
};

// NOLINTNEXTLINE
TEST_F(MpiOctree, SplitConstructionYieldsSameResult)
{
    int numBodies = 1000;
    std::vector<Body> bodies;
    bodies.resize(numBodies);

    Utility::Initialize(bodies.data(), numBodies, 0);

    Utility::OctreeCell root = Utility::OctreeCell(0, //// x center of root
                                                   0, //// y center of root
                                                   0.1374, //// z center of root
                                                   60 * SYSTEM_SIZE);

    OurOctree fullOctree;
    RankInfo fullOctreeRank { };
    fullOctreeRank.firstRootCell = std::numeric_limits<uint32_t>::max();
    fullOctreeRank.endRootCell = std::numeric_limits<uint32_t>::max();
    auto ourRootNodeId = fullOctree.addNode(0, 0, 0, root);

    for (int bcount = 1; bcount < numBodies; bcount++) {
        // Check if the body lies in the system
        if (fullOctree.getNode(ourRootNodeId).getOctreeCell().contains(bodies[bcount].pos)) {
            fullOctree.getNode(ourRootNodeId).insertBody(fullOctree, bodies[bcount], ourRootNodeId, 0, fullOctreeRank);
        }
    }

    int size = 32;

    std::vector<OurOctree> partialTrees;
    partialTrees.resize(size);
    for (int rank=0;rank<size;++rank) {
        OurOctree& partialTree = partialTrees[rank];

        RankInfo rankInfo;
        auto cellRangesByRank = CalculateRootCells(size).first;
        rankInfo.firstRootCell = cellRangesByRank[0].first;
        rankInfo.endRootCell = cellRangesByRank.back().second;
        auto cellRange = cellRangesByRank[rank];
        rankInfo.ownRootCellLo = cellRange.first;
        rankInfo.ownRootCellHi = cellRange.second;

        auto partialRootNodeId = partialTree.addNode(0, 0, 0, root);
        for (int bcount = 1; bcount < numBodies; bcount++) {
            // Check if the body lies in the system
            if (partialTree.getNode(partialRootNodeId).getOctreeCell().contains(bodies[bcount].pos)) {
                partialTree.getNode(partialRootNodeId).insertBody(partialTree, bodies[bcount], partialRootNodeId, 0, rankInfo);
            }
        }

        //compareOctree(partialTree, partialRootNodeId, fullOctree, ourRootNodeId);

        for (size_t pnid=0;pnid<partialTree.nodes.size();++pnid) {
            auto& partialNode = partialTree.nodes[pnid];
            uint32_t nodeOrd = partialNode.getNodeOrd();
            if (nodeOrd >= cellRange.first && nodeOrd < cellRange.second) {
                auto matchingNode = std::find_if(fullOctree.nodes.begin(), fullOctree.nodes.end(),
                    [nodeOrd](const auto& n) {
                        return n.getNodeOrd() == nodeOrd;
                });
                ASSERT_NE(matchingNode, fullOctree.nodes.end());
                compareOctree(fullOctree, std::distance(fullOctree.nodes.begin(), matchingNode)+1,
                    partialTree, pnid+1);
            }
        }
    }
}
