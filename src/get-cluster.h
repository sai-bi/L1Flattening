#ifndef GET_CLUSTER_H
#define GET_CLUSTER_H
#include "cv.h"
#include <opencv2/opencv.hpp>
#include "segment-image.h"
#include <queue>
using namespace std;
using namespace cv;
/*
* Segment an image
*
* Given an image, return its cluster result.
*
* im: image to segment.
* sigma: to smooth the image.
* c: constant for treshold function.
* min_size: minimum component size (enforced by post-processing stage).
* num_ccs: number of connected components in the segmentation.
* pixel_label: the label for each pixel, pixels in the same cluster have the same label
*/
void GetReflectanceCluster(const Mat_<Vec3d>& im,
                           float sigma, 
                           float c, 
                           int min_size,
                           int *num_ccs, 
                           Mat_<Vec3b>& output,
                           Mat_<int>& pixel_label,
                           int expected_cluster_num,
                           const Mat_<int>& region) {
    int width = im.cols;
    int height = im.rows;

    image<float> *r = new image<float>(width, height);
    image<float> *g = new image<float>(width, height);
    image<float> *b = new image<float>(width, height);
    pixel_label = Mat_<int>(height, width);

    // smooth each color channel  
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            r->access[y][x] = im(y, x)[0];
            r->access[y][x] = im(y, x)[1];
            r->access[y][x] = im(y, x)[2];
        }
    }
    image<float> *smooth_r = smooth(r, sigma);
    image<float> *smooth_g = smooth(g, sigma);
    image<float> *smooth_b = smooth(b, sigma);
    delete r;
    delete g;
    delete b;

    // build graph
    edge *edges = new edge[width*height * 4];
    int num = 0;
    int vertex_num = 0;
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            vertex_num++;

            if (x < width - 1) {
                edges[num].a = y * width + x;
                edges[num].b = y * width + (x + 1);
                edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x + 1, y);
                num++;
            }

            if (y < height - 1) {
                edges[num].a = y * width + x;
                edges[num].b = (y + 1) * width + x;
                edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x, y + 1);
                num++;
            }

            if ((x < width - 1) && (y < height - 1)) {
                edges[num].a = y * width + x;
                edges[num].b = (y + 1) * width + (x + 1);
                edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x + 1, y + 1);
                num++;
            }

            if ((x < width - 1) && (y > 0)) {
                edges[num].a = y * width + x;
                edges[num].b = (y - 1) * width + (x + 1);
                edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x + 1, y - 1);
                num++;
            }
        }
    }
    delete smooth_r;
    delete smooth_g;
    delete smooth_b;
    
    // restrict the maximum number of pixels in a cluster
    int max_pixel_number = vertex_num / (double)expected_cluster_num + 20;

    // segment
    universe *u = segment_graph(vertex_num, num, edges, c, max_pixel_number, region);

    // post process small components
    for (int i = 0; i < num; i++) {
        int a = u->find(edges[i].a);
        int b = u->find(edges[i].b);
        if ((a != b) && ((u->size(a) < min_size) || (u->size(b) < min_size)))
            u->join(a, b);
    }
    delete[] edges;
    *num_ccs = u->num_sets();

    // get the pixels in each cluster 
    map<int, int> index;
    int count = 1;
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int comp = u->find(y * width + x);
            if (index.count(comp) > 0){
                int temp = index[comp];
                pixel_label(y, x) = temp;
            }
            else{
                index[comp] = count;
                pixel_label(y, x) = count;
                count++;
			}
        }
    }

    /*
    output = Mat_<Vec3b>(height, width, Vec3b(0, 0, 0));
    rgb *colors = new rgb[width*height];
    for (int i = 0; i < width*height; i++){
        colors[i] = random_rgb();
    }

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int comp = u->find(y * width + x);
            // imRef(output, x, y) = colors[comp];
            output(y, x)[0] = colors[comp].b;
            output(y, x)[1] = colors[comp].g;
            output(y, x)[2] = colors[comp].r;
        }
    }
    */

    delete u;
}


#endif
