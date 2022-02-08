#include <gtest/gtest.h>
#include <vector>

#include "utility.h"
#include "ompbh.h"
#include "serial.h"

using OMP::Octree;
using OMP::OurOctree;
using OMP::OurOctreeNode;

namespace {

class OmpOctree : public ::testing::Test {
public:
    template <class OctreeA, class OctreeB>
    void compareOctree(OctreeA* left, OctreeB* right, const std::string& path = "R")
    {
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
            auto* leftChild = left->getChild(i);
            auto* rightChild = right->getChild(i);
            EXPECT_EQ(!!leftChild, !!rightChild) << path;
            if (leftChild && rightChild) {
                compareOctree(leftChild, rightChild, path + "." + std::to_string(i));
            }
        }
    }

    template <class OctreeA>
    void compareOctree(OctreeA* left, OurOctree& rightTree, OMP::NodeId nodeId, const std::string& path = "R")
    {
        auto right = &rightTree.getNode(nodeId);
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
            auto* leftChild = left->getChild(i);
            auto rightChild = right->getChild(i);
            EXPECT_EQ(!!leftChild, !!rightChild) << path;
            if (leftChild && rightChild) {
                compareOctree(leftChild, rightTree, rightChild, path + "." + std::to_string(i));
            }
        }
    }
};

// NOLINTNEXTLINE
TEST_F(OmpOctree, CanCompareOctreeWithOctree)
{
    int numBodies = 1024;
    std::vector<Body> bodies;
    bodies.resize(numBodies);

    Utility::Initialize(bodies.data(), numBodies, 0);

    Utility::OctreeCell&& root = Utility::OctreeCell(0, //// x center of root
        0, //// y center of root
        0.1374, //// z center of root
        60 * SYSTEM_SIZE);
    Octree* tree = new Octree(std::move(root));

    Utility::OctreeCell root2 = Utility::OctreeCell(0, //// x center of root
        0, //// y center of root
        0.1374, //// z center of root
        60 * SYSTEM_SIZE);
    Octree* tree2 = new Octree(root2);

    for (int bcount = 1; bcount < numBodies; bcount++) {
        // Check if the body lies in the system
        if (tree->getOctreeCell().contains(bodies[bcount].pos)) {
            tree->InsertBody(&bodies[bcount]);
        }
    }
    for (int bcount = 1; bcount < numBodies; bcount++) {
        // Check if the body lies in the system
        if (tree2->getOctreeCell().contains(bodies[bcount].pos)) {
            tree2->InsertBody(&bodies[bcount]);
        }
    }

    compareOctree(tree, tree2);
}

// NOLINTNEXTLINE
TEST_F(OmpOctree, CanCompareOctreeWithOurOctree)
{
    int numBodies = 1024;
    std::vector<Body> bodies;
    bodies.resize(numBodies);

    Utility::Initialize(bodies.data(), numBodies, 0);

    Utility::OctreeCell root = Utility::OctreeCell(0, //// x center of root
        0, //// y center of root
        0.1374, //// z center of root
        60 * SYSTEM_SIZE);
    Octree* tree = new Octree(root);

    OurOctree ourOctree;
    auto ourRootNodeId = ourOctree.addNode(root);

    for (int bcount = 1; bcount < numBodies; bcount++) {
        // Check if the body lies in the system
        if (tree->getOctreeCell().contains(bodies[bcount].pos)) {
            tree->InsertBody(&bodies[bcount]);
        }
    }
    for (int bcount = 1; bcount < numBodies; bcount++) {
        // Check if the body lies in the system
        if (ourOctree.getNode(ourRootNodeId).getOctreeCell().contains(bodies[bcount].pos)) {
            ourOctree.getNode(ourRootNodeId).insertBody(ourOctree, bodies[bcount], ourRootNodeId);
        }
    }

    compareOctree(tree, ourOctree, ourRootNodeId);
}

class OmpSolution : public ::testing::Test {
public:
    int numBodies { 8192 };
    int numSteps { 500 };
    int seed { 0 };

    double maxDelta { 1.0 };

    void compareBody(const Body& left, const Body& right, const std::string& name)
    {
        ASSERT_NEAR(left.m, right.m, maxDelta) << name;

        ASSERT_NEAR(left.pos.x, right.pos.x, maxDelta) << name;
        ASSERT_NEAR(left.pos.y, right.pos.y, maxDelta) << name;
        ASSERT_NEAR(left.pos.z, right.pos.z, maxDelta) << name;

        ASSERT_NEAR(left.vel.x, right.vel.x, maxDelta) << name;
        ASSERT_NEAR(left.vel.y, right.vel.y, maxDelta) << name;
        ASSERT_NEAR(left.vel.z, right.vel.z, maxDelta) << name;

        ASSERT_NEAR(left.acc.x, right.acc.x, maxDelta) << name;
        ASSERT_NEAR(left.acc.y, right.acc.y, maxDelta) << name;
        ASSERT_NEAR(left.acc.z, right.acc.z, maxDelta) << name;
    }
};

// NOLINTNEXTLINE
TEST_F(OmpSolution, OmpSolutionMatchesSerialImplementation)
{
    GTEST_SKIP();

    std::vector<Body> bodiesSerial, bodiesOmp;
    bodiesSerial.resize(numBodies);
    bodiesOmp.resize(numBodies);

    Utility::Initialize(bodiesSerial.data(), numBodies, seed);
    Utility::Initialize(bodiesOmp.data(), numBodies, seed);

    Serial::Solve(bodiesSerial.data(), numBodies, numSteps);
    OMP::Solve(bodiesOmp.data(), numBodies, numSteps);

    //for (int k=0;k<numSteps;++k) {

    for (int i = 0; i < numBodies; ++i) {
        compareBody(bodiesSerial[i], bodiesOmp[i], "body " + std::to_string(i));
    }
    //}
}

} // namespace
