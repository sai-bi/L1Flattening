#include "flattening.h"
using namespace cv;
using namespace Eigen;

SpMat BuildWindowVariationMatrix(const cv::Mat_<cv::Vec3b>& image,
	int window_size){
	double mu = 10;
	double ga = 120;
	double sigma = 0.5;

	Mat_<Vec3b> lab_image;
	cvtColor(image, lab_image, cv::COLOR_BGR2Lab);

	int height = image.rows;
	int width = image.cols;



}