//
// Created by shubham on 08.06.21.
//

#include <cmath>
#include "utility.h"
#include <random>
#include <iostream>

double Utility::magnitude(vec3* v)
{
    return Utility::magnitude( v->x, v->y, v->z);
}

double Utility::magnitude( double x, double y, double z)
{
    return sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
}

void Utility::Initialize(Body* bods, const int& numBodies, int seed)
{
    using std::uniform_real_distribution;
    uniform_real_distribution<double> randAngle (0.0, 200.0*PI);
    uniform_real_distribution<double> randRadius (INNER_BOUND, SYSTEM_SIZE);
    uniform_real_distribution<double> randHeight (0.0, SYSTEM_THICKNESS);
    std::default_random_engine gen (seed);
    double angle;
    double radius;
    double velocity;
//    double binary_velocity;
    Body *current;

    //STARS
//    binary_velocity = 0.67*sqrt((G*SOLAR_MASS)/(4*BINARY_SEPARATION*TO_METERS));

    //STAR 1
    current = &bods[0];
    current->pos.x = 0.0;///-BINARY_SEPARATION;
    current->pos.y = 0.0;
    current->pos.z = 0.0;
    current->vel.x = 0.0;
    current->vel.y = 0.0;//binary_velocity;
    current->vel.z = 0.0;
    current->m = SOLAR_MASS;
    //STAR 2
    /*
    current = bods + 1;
    current->position.x = BINARY_SEPARATION;
    current->position.y = 0.0;
    current->position.z = 0.0;
    current->velocity.x = 0.0;
    current->velocity.y = -binary_velocity;
    current->velocity.z = 0.0;
    current->mass = SOLAR_MASS;
    */

    ///STARTS AT NUMBER OF STARS///
    double totalExtraMass = 0.0;
    for (int index=1; index<numBodies; index++)
    {
        angle = randAngle(gen);
        radius = sqrt(SYSTEM_SIZE)*sqrt(randRadius(gen));
        velocity = pow(((G*(SOLAR_MASS+((radius-INNER_BOUND)/SYSTEM_SIZE)*EXTRA_MASS*SOLAR_MASS))
                        / (radius*TO_METERS)), 0.5);
        current = &bods[index];
        current->pos.x =  radius*cos(angle);
        current->pos.y =  radius*sin(angle);
        current->pos.z =  randHeight(gen)-SYSTEM_THICKNESS/2.;
        current->vel.x =  velocity*sin(angle);
        current->vel.y = -velocity*cos(angle);
        current->vel.z =  0.0;
        current->m = (EXTRA_MASS*SOLAR_MASS)/numBodies;
        totalExtraMass += (EXTRA_MASS*SOLAR_MASS)/numBodies;
    }

}

Utility::OctreeCell& Utility::OctreeCell::operator=(const OctreeCell& o){
    if(this!= &o){
        midpoint = vec3(o.midpoint);
        radius = o.radius;
    }
    return *this;
}

double Utility::OctreeCell::getRadius()const{
    return radius;
}

vec3 Utility::OctreeCell::getMidpoint()const{
    return midpoint;
}

bool Utility::OctreeCell::contains(const vec3& p) const {
    bool xcond = p.x<=(midpoint.x+(radius/2.0)) && p.x>=(midpoint.x-(radius/2.0));
    bool ycond = p.y<=(midpoint.y+(radius/2.0)) && p.y>=(midpoint.y-(radius/2.0));
    bool zcond = p.z<=(midpoint.z+(radius/2.0)) && p.z>=(midpoint.z-(radius/2.0));
    return xcond && ycond && zcond;
}

bool Utility::OctreeCell::IsUp(Body& b)const {
    return b.pos.z>=midpoint.z;
}
bool Utility::OctreeCell::IsEast(Body& b)const{
    return b.pos.x>=midpoint.x;
}
bool Utility::OctreeCell::IsNorth(Body& b)const{
    return b.pos.y>=midpoint.y;
}

Utility::OctreeCell Utility::OctreeCell::UpNorthWest() const {
    double halfRadius = radius/4.0;
    return OctreeCell(midpoint.x - halfRadius, midpoint.y + halfRadius,
                      midpoint.z + halfRadius, 2*halfRadius);
}
Utility::OctreeCell Utility::OctreeCell::UpNorthEast() const {
    double halfRadius = radius/4.0;
    return OctreeCell(midpoint.x + halfRadius, midpoint.y + halfRadius,
                      midpoint.z + halfRadius, 2*halfRadius);
}

Utility::OctreeCell Utility::OctreeCell::UpSouthEast() const {
    double halfRadius = radius/4.0;
    return OctreeCell(midpoint.x + halfRadius, midpoint.y - halfRadius,
                      midpoint.z + halfRadius, 2*halfRadius);
}

Utility::OctreeCell Utility::OctreeCell::UpSouthWest() const {
    double halfRadius = radius/4.0;
    return OctreeCell(midpoint.x - halfRadius, midpoint.y - halfRadius,
                      midpoint.z + halfRadius, 2*halfRadius);
}

Utility::OctreeCell Utility::OctreeCell::DownNorthWest() const {
    double halfRadius = radius/4.0;
    return OctreeCell(midpoint.x - halfRadius, midpoint.y + halfRadius,
                      midpoint.z - halfRadius, 2*halfRadius);
}
Utility::OctreeCell Utility::OctreeCell::DownNorthEast() const {
    double halfRadius = radius/4.0;
    return OctreeCell(midpoint.x + halfRadius, midpoint.y + halfRadius,
                      midpoint.z - halfRadius, 2*halfRadius);
}

Utility::OctreeCell Utility::OctreeCell::DownSouthEast() const {
    double halfRadius = radius/4.0;
    return OctreeCell(midpoint.x + halfRadius, midpoint.y - halfRadius,
                      midpoint.z - halfRadius, 2*halfRadius);
}

Utility::OctreeCell Utility::OctreeCell::DownSouthWest() const {
    double halfRadius = radius/4.0;
    return OctreeCell(midpoint.x - halfRadius, midpoint.y - halfRadius,
                      midpoint.z - halfRadius, 2*halfRadius);
}

void Utility::BodyInteract(Body* a, Body* b)
{
    //////////////////////////////////////////////////////
    //////// This is for pairwise interaction ////////////
    //////////////////////////////////////////////////////
    vec3 diff = a->pos - b->pos;
    diff = diff*TO_METERS;

    double dist = Utility::magnitude(&diff);
    double F = (G*a->m*b->m) / ((dist*dist + SOFTENING*SOFTENING) * dist);

    a->acc = a->acc - diff*(F/a->m);

    b->acc = b->acc + diff*(F/b->m);
}

void Utility::BodyInteract(Body* target, const Body* other, bool singlePart)
{
    //////////////////////////////////////////////////////////////////
    //////// This is for Tree interaction ////////////////////////////
    //////// Note the difference here only target is updated /////////
    //////////////////////////////////////////////////////////////////
    vec3 diff;
    diff = (target->pos - other->pos)*TO_METERS;
    double dist = Utility::magnitude(&diff);

    // this test can be true only when singlePart is true
    // (is always false when singlePart is false = when target is not external)
    if (dist == 0)
        // avoid division by 0
        return;

    double F = (G * target->m * other->m) / (( dist*dist + SOFTENING*SOFTENING) * dist);

    target->acc = target->acc - (diff*F)/target->m;

    //Friction
#if ENABLE_FRICTION
    if (singlePart)
    {
        double friction = 0.5/pow(2.0,FRICTION_FACTOR*(
                ((dist+SOFTENING))/(TO_METERS)));
        //	cout << friction << "\n";
        if (friction>0.0001 && ENABLE_FRICTION)
        {
            target->acc = target->acc + (other->vel-target->vel)*(friction/2);
        }
    }
#else
    (void)singlePart;
#endif
}

void Utility::Integrate(Body *b) {
    ////////////////////////////////////////////////
    ///// Update the velocity and //////////////////
    /// position of the body ///////////////////////
    ////////////////////////////////////////////////
    b->vel = b->vel + (b->acc)*TIME_STEP;
    b->acc = 0.0;
    b->pos = b->pos + (b->vel)*(TIME_STEP / TO_METERS);
}

void Utility::ParseInputs(int argc, char* argv[], int* numBodies, int* numSteps){
    if (argc>=4)
    {
        if (std::string(argv[1]) == "-n")
        {
            *numBodies = std::stoi(argv[2]);
        }
        else
        {
            fprintf(stderr, "Wrong argument - %s\n", argv[1]);
            exit(EXIT_FAILURE);
        }
        if (std::string(argv[3]) == "-t"){
            *numSteps = std::stoi(argv[4]);
        }
        else
        {
            fprintf(stderr, "Wrong argument - %s\n", argv[3]);
            exit(EXIT_FAILURE);
        }
    }
    else
    {
        fprintf(stderr, "Usage:./BarnesHut -n <numBodies> -t <numSteps>\n");
        exit(EXIT_FAILURE);
    }
}

void Utility::PrintBodies(Body* bodies, int numBodies){
    std::cout<<"# Bodies = "<<numBodies<<std::endl;
    std::cout<<"Positions (x, y, z)"<<std::endl;
    for(int i=0; i<numBodies; ++i){
        std::cout<<bodies[i].pos.x<<", "<<bodies[i].pos.y<<", "<<bodies[i].pos.z<<std::endl;
    }
}