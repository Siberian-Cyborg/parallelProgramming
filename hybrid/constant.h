#ifndef CONSTANT_H_
#define CONSTANT_H_

#include <string>
#include <sstream>

#define PI      3.14159265358979323846
#define TO_METERS 1.496e11 // Meters in an AU
#define SYSTEM_SIZE 3.5    // Farthest particles in AU
#define SYSTEM_THICKNESS 0.08  //  Thickness in AU
#define INNER_BOUND 0.3    // Closest particles to center in AU
#define SOFTENING (0.015*TO_METERS) // Softens particles interactions at close distances
#define SOLAR_MASS (2.0e30)  // in kg
#define BINARY_SEPARATION 0.07 // AU (only applies when binary code uncommented)
#define EXTRA_MASS 1.5 // 0.02 Disk mask as a portion of center star/black hole mass
#define ENABLE_FRICTION 0 // For experimentation only. Will probably cause weird results
#define FRICTION_FACTOR 25.0 // Only applies if friction is enabled
#define MAX_DISTANCE 0.75 //2.0  Barnes-Hut Distance approximation factor
#define G (6.67408e-11) // The gravitational constant
#define TIME_STEP (1*32*1024) //(1*128*1024) Simulated time between integration steps, in seconds

/// For visualization
#define RENDER_INTERVAL 100 // Choose higher value when simulating for lot of steps
#define SIZE_SCALE 2.5 //1.1
#define IMAGEWIDTH 1024
#define IMAGEHEIGHT 1024
#define PIXEL_SIZE 8
#define MAX_VEL_BOUND 40000.0  // Both in km/s
#define MIN_VEL_BOUND 14000.0
#define PARTICLE_BRIGHTNESS 0.52 //can be varied for different resolution
#define PARTICLE_SHARPNESS 1.0

struct vec3{
    vec3(){}
    vec3(double x, double y, double z): x(x), y(y), z(z){}
    vec3 operator+(vec3 a){
        return vec3(this->x + a.x, this->y + a.y, this->z + a.z);
    }
    vec3 operator-(vec3 a){
        return vec3(this->x - a.x, this->y - a.y, this->z - a.z);
    }
    vec3 operator*(double a){
        return vec3(this->x * a, this->y * a, this->z * a);
    }
    vec3 operator*(int a){
        return vec3(this->x * (double)a, this->y * (double)a, this->z * (double)a);
    }
    vec3 operator/(double a){
        return vec3(this->x / a, this->y / a, this->z / a);
    }
    vec3& operator=(double a){
        this->x = a;
        this->y = a;
        this->z = a;
        return *this;
    }

    bool operator==(const vec3& a){
        if(this==&a){
            return true;
        } else{
            if((x==a.x)&(y==a.y)&&(z==a.z)){
                return true;
            }
            else{
                return false;
            }
        }
    }
    double x, y, z;
};

//Body Structure containing 7 variables
struct Body
{
    Body(): m(0){}
    vec3 pos, vel, acc;
    float m;
    std::string print(){
        std::ostringstream ss;
        ss<<"Mass = "<<m<<"\n";
        ss<<"Pos = ("<<pos.x<<", "<<pos.y<<", "<<pos.z<<")"<<"\n";
        ss<<"Vel = ("<<vel.x<<", "<<vel.y<<", "<<vel.z<<")"<<"\n";
        ss<<"Acc = ("<<acc.x<<", "<<acc.y<<", "<<acc.z<<")"<<"\n";
        std::string str(ss.str());
        return str;
    }
};

//// For visualization
struct color
{
    double r, g, b;
};


#endif /* CONSTANTS_H_ */
