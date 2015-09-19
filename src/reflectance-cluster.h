#ifndef REFLECTANCE_H
#define REFLECTANCE_H
#include <vector>
#include <opencv2/opencv.hpp>
#include "segment-image.h"
#include <algorithm>
using namespace std;
using namespace cv;

class ReflectanceCluster{
    public:
        ReflectanceCluster(){
            cluster_size_ = 0;
            cluster_center_ = cv::Point2i(-1,-1); 
        }
        cv::Point2i GetClusterCenter(){
            if(cluster_center_.x >= 0 && cluster_center_.y >= 0){
                return cluster_center_; 
            } 
            double x = 0;
            double y = 0;
            for(int i = 0;i < pixel_locations_.size();i++){
                x += pixel_locations_[i].x;
                y += pixel_locations_[i].y;
            } 
            assert(pixel_locations_.size()!=0);
            cluster_center_.x = (int)(x / pixel_locations_.size());
            cluster_center_.y = (int)(y / pixel_locations_.size());
            return cluster_center_;
        }
        int GetClusterSize(){
            return cluster_size_;
        }
        void AddPixel(cv::Point2i pixel){
            pixel_locations_.push_back(pixel); 
            cluster_size_++;
        } 
        std::vector<cv::Point2i> GetPixelLocations(){
            return pixel_locations_;
        }
    private:
        int cluster_size_;
        cv::Point2i cluster_center_;
        std::vector<cv::Point2i> pixel_locations_; 
};

/*
 * Get the weight for Reflectance Sparse Matrix between each two clusters
 */
void GetPairwiseWeight(vector<ReflectanceCluster>& clusters,
                               const Mat_<Vec3b>& image,
                               const Mat_<int>& cluster_same_region,
                               Mat_<double>& weight){
    int cluster_num = clusters.size();
    weight = Mat_<double>(cluster_num, cluster_num, 0.0);
    Mat_<double> feature(cluster_num, 3, 0.0);

    double theta_l = 0.1;
    double theta_c = 0.0025;
    for(int i = 0;i < cluster_num; ++i){
        vector<Point2i> pixels_in_clusters = clusters[i].GetPixelLocations();
        for(int j = 0;j < pixels_in_clusters.size(); ++j){
            double r = image(pixels_in_clusters[j].x, pixels_in_clusters[j].y)[2];
            double g = image(pixels_in_clusters[j].x, pixels_in_clusters[j].y)[1];
            double b = image(pixels_in_clusters[j].x, pixels_in_clusters[j].y)[0];
            double shrink = 0.0001;
            feature(i, 0) += ((r + g + b + shrink) / 3.0);
            feature(i, 1) += (r / (r + g + b + shrink));
            feature(i, 2) += (g / (r + g + b + shrink)); 
        }    
        feature(i, 0) /= (pixels_in_clusters.size() * theta_l);
        feature(i, 1) /= (pixels_in_clusters.size() * theta_c);
        feature(i, 2) /= (pixels_in_clusters.size() * theta_c);
    }

    Mat_<double> row_sum(1, 3, 0.0);
    row_sum(0) = sum(feature.col(0))[0];
    row_sum(1) = sum(feature.col(1))[0];
    row_sum(2) = sum(feature.col(2))[0];
    row_sum = (1.0 / feature.rows) * row_sum;

    double variance = 0;
    for (int i = 0; i < row_sum.rows; ++i){
        double temp = pow(norm(feature.row(i) - row_sum), 2.0);
        variance += temp;
    }
    variance = variance / (cluster_num - 1);
    for (int i = 0; i < cluster_num; ++i){
        for (int j = i + 1; j < cluster_num; ++j) {
            Mat_<double> row_1 = feature.row(i);
            Mat_<double> row_2 = feature.row(j);
            double temp = pow(norm(row_1 - row_2), 2.0);
            double w = exp(-0.5 * temp / variance);
            // if two clusters belong to different regions, this weight should be smaller
            /*
            if (cluster_same_region(i, j) == 0){
                w = 0.01 * w;
            }
            */
            weight(i, j) = w;
            weight(j, i) = w;
        }
    }

}

/*
 * Construct adjacent matrix for clusters.
 * The value in (i, j) represents the number of pair of pixels that on boundary of 
 * cluster i and cluster j, 0 represents they are not adjacent.
 * 
 * Input:
 *      pixel_label: the label of superpixels for each pixel
 *      region: the region label for each pixel
 *      cluster_num: number of superpixels
 * Output: 
 *      adjacent: adjacent matrix for superpixels
 *      cluster_same_region: whether two superpixels belong to same region
 */
void GetClusterAdjacentMatrix(const Mat_<int>& pixel_label, 
                              int cluster_num, 
                              const Mat_<int>& region,
                              Mat_<int>& adjacent,
                              Mat_<int>& cluster_same_region){
    int width = pixel_label.cols;
    int height = pixel_label.rows;
    adjacent = Mat_<int>(cluster_num, cluster_num, 0);
    cluster_same_region = Mat_<int>(cluster_num, cluster_num, 1);

    for (int i = 0; i < height - 1; ++i) {
        for(int j = 0; j < width - 1; ++j) {
            int label_1 = pixel_label(i, j); 
            int label_2; 
            
            if(j + 1 < width){
                label_2 = pixel_label(i, j + 1);
                if (label_1 != label_2){
                    adjacent(label_1, label_2) += 1;
                    adjacent(label_2, label_1) += 1;
                    if (region(i, j) != region(i, j + 1)){
                        cluster_same_region(label_1, label_2) = 0;
                        cluster_same_region(label_2, label_1) = 0;
                    }
                }
            }
            if (i + 1 < height){
                label_2 = pixel_label(i + 1, j);
                if (label_1 != label_2){
                    adjacent(label_1, label_2) += 1;
                    adjacent(label_2, label_1) += 1;
                    if (region(i, j) != region(i, j + 1)){
                        cluster_same_region(label_1, label_2) = 0;
                        cluster_same_region(label_2, label_1) = 0;
                    }
                }
            }
            if (i - 1 >= 0){
                label_2 = pixel_label(i - 1, j);
                if (label_1 != label_2){
                    adjacent(label_1, label_2) += 1;
                    adjacent(label_2, label_1) += 1;
                    if (region(i, j) != region(i, j + 1)){
                        cluster_same_region(label_1, label_2) = 0;
                        cluster_same_region(label_2, label_1) = 0;
                    }
                }
            }
            if (j - 1 >= 0){
                label_2 = pixel_label(i, j-1);
                if (label_1 != label_2){
                    adjacent(label_1, label_2) += 1;
                    adjacent(label_2, label_1) += 1;
                    if (region(i, j) != region(i, j + 1)){
                        cluster_same_region(label_1, label_2) = 0;
                        cluster_same_region(label_2, label_1) = 0;
                    }
                }
            }
        }
    }
}

#endif
