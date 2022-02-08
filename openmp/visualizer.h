//
// Created by shubham on 08.06.21.
//

#ifndef BARNES_HUT_VISUALIZER_H
#define BARNES_HUT_VISUALIZER_H
#include "constant.h"

namespace Visualizer{
    void SaveFrame(char* image, double* imageArray, Body* b, int numBodies, int step);
    double MeterToPixel(double p, int size);
    void FrameSetup(char* image, double* imageArray);
    void PopulateImage(Body* b, double* imageArray, int numBodies);
    void ImageColorSetup(double x, double y, double vMag, double* imageArray);
    void PixelColor(int x, int y, const struct color& c, double f, double* imageArray);
    double RangeAdjuster(double x);
    void WriteFrame(char* data, double* imageArray, int step);
}

#endif //BARNES_HUT_VISUALIZER_H
