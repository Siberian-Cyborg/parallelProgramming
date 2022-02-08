//
// Created by shubham on 08.06.21.
//
#include "visualizer.h"
#include "constant.h"
#include "utility.h"
#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>

void Visualizer::SaveFrame(char* image, double* imageArray, Body* b, int numBodies, int step){
    std::cout<<"Writing Frame "<<step<<std::endl;
    Visualizer::FrameSetup(image, imageArray);
    Visualizer::PopulateImage(b, imageArray, numBodies);
    Visualizer::WriteFrame(image, imageArray, step);
}

double Visualizer::MeterToPixel(double p, int size){
    return (size/2.0)*(1.0+p/(SYSTEM_SIZE*SIZE_SCALE));
}

void Visualizer::FrameSetup(char* image, double* imageArray){
    std::memset(image, 0, IMAGEWIDTH*IMAGEHEIGHT*3);
    std::memset(imageArray, 0, IMAGEWIDTH*IMAGEHEIGHT*3*sizeof(double));
}

void Visualizer::PopulateImage(Body* b, double* imageArray, int numBodies){
    for(int index=0; index<numBodies; index++)
    {
        int x = Visualizer::MeterToPixel((&b[index])->pos.x, IMAGEWIDTH);
        int y = Visualizer::MeterToPixel((&b[index])->pos.y, IMAGEHEIGHT);
        if (x>PIXEL_SIZE && x<IMAGEWIDTH-PIXEL_SIZE &&
            y>PIXEL_SIZE && y<IMAGEHEIGHT-PIXEL_SIZE)
        {
            double vMag = Utility::magnitude(&((&b[index])->vel));
            ImageColorSetup((&b[index])->pos.x, (&b[index])->pos.y, vMag, imageArray);
        }
    }
}
void Visualizer::ImageColorSetup(double x, double y, double vMag, double* imageArray){
    constexpr double velocityMax = MAX_VEL_BOUND;
    double velocityMin = sqrt(0.8*(G*(SOLAR_MASS+EXTRA_MASS*SOLAR_MASS))/
                                        (SYSTEM_SIZE*TO_METERS));
    if (vMag < velocityMin)
        return;
    const double vPortion = sqrt((vMag-velocityMin) / velocityMax);
    color c;
    c.r = Visualizer::RangeAdjuster(4*(vPortion-0.333));
    c.g = Visualizer::RangeAdjuster(fmin(4*vPortion,4.0*(1.0-vPortion)));
    c.b = Visualizer::RangeAdjuster(4*(0.5-vPortion));
    for (int i=-PIXEL_SIZE/2; i<PIXEL_SIZE/2; i++)
    {
        for (int j=-PIXEL_SIZE/2; j<PIXEL_SIZE/2; j++)
        {
            double xP = floor(Visualizer::MeterToPixel(x, IMAGEWIDTH));
            double yP = floor(Visualizer::MeterToPixel(y, IMAGEHEIGHT));
            double cFactor = PARTICLE_BRIGHTNESS /
                             (pow(exp(pow(PARTICLE_SHARPNESS*
                                          (xP+i-Visualizer::MeterToPixel(x, IMAGEWIDTH)),2.0))
                                  + exp(pow(PARTICLE_SHARPNESS*
                                            (yP+j-Visualizer::MeterToPixel(y, IMAGEHEIGHT)),2.0)),/*1.25*/0.75)+1.0);
            Visualizer::PixelColor(int(xP+i),int(yP+j),c, cFactor, imageArray);
        }
    }
}

void Visualizer::PixelColor(int x, int y, const struct color& c, double f, double* imageArray){
    int pix = 3*(x+IMAGEWIDTH*y);
    imageArray[pix+0] += c.r*f;
    imageArray[pix+1] += c.g*f;
    imageArray[pix+2] += c.b*f;
}
double Visualizer::RangeAdjuster(double x){
    return fmax(fmin(x,1.0),0.0);
}
void Visualizer::WriteFrame(char* data, double* imageArray, int step){
    for (int i=0; i<IMAGEWIDTH*IMAGEHEIGHT*3; i++)
    {
        data[i] = int(255.0*Visualizer::RangeAdjuster(imageArray[i]));
    }

    int frame = step/RENDER_INTERVAL + 1;
    char name[128];
    sprintf(name, "images/Step%05i.ppm", frame);
    std::ofstream file (name, std::ofstream::binary);

    if (file.is_open())
    {
        file << "P6\n" << IMAGEWIDTH << " " << IMAGEHEIGHT << "\n" << "255\n";
        file.write(data, IMAGEWIDTH*IMAGEHEIGHT*3);
        file.close();
    }else{
        std::cout<<"Make an image directory at the location"<<std::endl;
    }
}


