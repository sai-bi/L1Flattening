#ifndef SHADING_H
#define SHADING_H
#include "reflectance-cluster.h"
using namespace cv;
using namespace std;


/*
 * Given the sovled reflectance, get the real reflectance and shading.
 */
void GetShadingReflectance(const Mat_<Vec3b>& image, 
                           const Mat_<double>& solved_reflectance,
                           vector<ReflectanceCluster>& clusters,
                           Mat_<Vec3b>& reflectance, Mat_<double>& shading){
    int width = image.cols;
    int height = image.rows;
    double max_reflectance = *max_element(solved_reflectance.begin(), solved_reflectance.end());
    cout << "Max reflectance: " << max_reflectance << endl;
    reflectance = Mat_<Vec3b>(height, width);
    shading = Mat_<double>(height, width, 0.0);
    
    for (int i = 0; i < solved_reflectance.rows; ++i) {
        double r = 0;
        double g = 0;
        double b = 0;
        vector<Point2i> pixels_in_cluster = clusters[i].GetPixelLocations();
        for (int j = 0; j < pixels_in_cluster.size(); ++j){
            r = r + image(pixels_in_cluster[j].x, pixels_in_cluster[j].y)[2];
            g = g + image(pixels_in_cluster[j].x, pixels_in_cluster[j].y)[1];
            b = b + image(pixels_in_cluster[j].x, pixels_in_cluster[j].y)[0];
        }
        r = r / pixels_in_cluster.size();
        g = g / pixels_in_cluster.size();
        b = b / pixels_in_cluster.size();

        double temp = exp(solved_reflectance(i));
        double temp_1 = exp(solved_reflectance(i));
        double scale_ratio = 1;
        double reflectance_r = scale_ratio * 3.0 * r * temp / (r + g + b);
        double reflectance_g = scale_ratio * 3.0 * g * temp / (r + g + b);
        double reflectance_b = scale_ratio * 3.0 * b * temp / (r + g + b);
        double shading_r = (r + g + b) / (3.0 * temp);
        for (int j = 0; j < pixels_in_cluster.size(); ++j) {
            int x = pixels_in_cluster[j].x;
            int y = pixels_in_cluster[j].y;
            reflectance(x, y) = Vec3b((uchar)reflectance_b, (uchar)reflectance_g, (uchar)reflectance_r);
            // shading(x, y) = shading_r;
            shading(x, y) = (image(x, y)[0] + image(x, y)[1] + image(x, y)[2]) / (3.0 * temp);
        }
    }
}


#endif
