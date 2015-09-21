#ifndef MAT_IMAGE_H
#define MAT_IMAGE_H
#include <opencv2/opencv.hpp>
#include "image.h"
#include "segment-graph.h"
#include "misc.h"
using namespace cv;
using namespace std;

/*
* Given an image in Mat_<Vec3b>, turn it into image<rgb>
*/
image<rgb>* MatToImage(const Mat_<Vec3b>& input){
    int height = input.rows;
    int width = input.cols;
    image<rgb>* output = new image<rgb>(width, height);

    for (int x = 0; x < height; x++){
        for (int y = 0; y < width; y++){
            rgb color;
            color.b = input(x, y)[0];
            color.g = input(x, y)[1];
            color.r = input(x, y)[2];
            output->access[x][y] = color;
        }
    }
    return output;
}

/*
 * Given an image stored in image<rgb>, turn it into Mat_<Vec3b>.
*/
Mat_<Vec3b> ImageToMat(const image<rgb>& input){
    int height = input.height();
    int width = input.width();
    Mat_<Vec3b> output(height, width, Vec3b(0, 0, 0));

    for (int x = 0; x < height; x++){
        for (int y = 0; y < width; y++){
            output(x, y)[0] = input.access[x][y].b;
            output(x, y)[1] = input.access[x][y].g;
            output(x, y)[2] = input.access[x][y].r;
        }
    }
    return output;
}

#endif