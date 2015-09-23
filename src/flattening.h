#ifndef FLATTENING_H
#define FLATTENING_H
#include <Eigen/sparse>
#include "cv.h"
#include <opencv2/opencv.hpp>
#include <vector>

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> Triple;

SpMat& BuildWindowVariationMatrix(const cv::Mat_<cv::Vec3b>& image,
	int window_size);

#endif