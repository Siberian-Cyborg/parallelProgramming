//
// Created by shubham on 08.06.21.
//

#ifndef BARNES_HUT_UTILITY_H
#define BARNES_HUT_UTILITY_H

#include "constant.h"
#include <utility>

namespace Utility{
    double magnitude(vec3* v);
    double magnitude( double x, double y, double z);
    void Initialize(Body* bods, const int& numBodies, int seed=0);
    void PrintBodies(Body* bodies, int numBodies);
    ////////////////////////////////////////////////////////////////
    /// You can modify the below function but you should ///////////
    /// implement it in the your submission file and use ///////////
    /// corresponding namespace for distinction ////////////////////
    ////////////////////////////////////////////////////////////////

    void BodyInteract(Body* a, Body* b);
    void BodyInteract(Body* target, Body* other, bool singlePart);
    void Integrate(Body* b);
    void ParseInputs(int argc, char* argv[], int* numBodies, int* numSteps);

    ///////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////

    class OctreeCell{
    ////////////////////////////////////////////////////////////////////
    /// You are allowed to modify the implementation of this class /////
    /// If you wish to do so, make an implementation in your ///////////
    /// parallel file under the same namespace /////////////////////////
    ////////////////////////////////////////////////////////////////////
    private:
        vec3 midpoint;
        double radius;

    public:
        OctreeCell(){
            midpoint = vec3(0,0,0);
            radius = 0.;
        }
        OctreeCell(double x, double y, double z, double rad): midpoint(x,y,z), radius(rad){}
        OctreeCell(OctreeCell&& o): midpoint(std::move(o.midpoint)), radius(std::move(o.radius)){}
        OctreeCell(const OctreeCell& o): midpoint(o.midpoint), radius(o.radius){}
        ~OctreeCell(){}
        OctreeCell& operator=(const OctreeCell&);
        double getRadius()const;
        vec3 getMidpoint()const;
        bool contains(const vec3& p) const;
        OctreeCell UpNorthWest() const;
        OctreeCell UpNorthEast() const;

        OctreeCell UpSouthWest() const;
        OctreeCell UpSouthEast() const;
        OctreeCell DownNorthWest() const;
        OctreeCell DownNorthEast() const;
        OctreeCell DownSouthWest() const;
        OctreeCell DownSouthEast() const;
        bool IsUp(Body& b)const;
        bool IsEast(Body& b)const;
        bool IsNorth(Body& b)const;

    };
}

#endif //BARNES_HUT_UTILITY_H
