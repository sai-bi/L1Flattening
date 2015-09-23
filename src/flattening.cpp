#include "flattening.h"
#include <cmath>
using namespace cv;
using namespace Eigen;
using namespace std;

SpMat& BuildWindowVariationMatrix(const cv::Mat_<cv::Vec3b>& image,
	int window_size){
	double mu = 10;
	double ga = 120;
	double sigma = 0.5;

	Mat_<Vec3d> lab_image;
	cvtColor(image, lab_image, cv::COLOR_BGR2Lab);
	lab_image = 1.0 / 255 * lab_image;

	int height = image.rows;
	int width = image.cols;

	vector<Triple> triples;
	int row_count = 0;
	int pixel_num = height * width;
	for (int i = 0; i < width; i++){
		for (int j = 0; j < height; j++){
			for (int p = 0; p < i + window_size && p < width; p++){
				for (int q = 0; q < j + window_size && q < height; q++){
					if (i == p && j == q){
						continue;
					}
					int center_id = (i - 1) * height + j;
					int neigh_id = (p - 1)*height + q;
					
					double l1 = lab_image(j, i)[0];
					double a1 = lab_image(j, i)[1];
					double b1 = lab_image(j, i)[2];
					double l2 = lab_image(q, p)[0];
					double a2 = lab_image(q, p)[1];
					double b2 = lab_image(q, p)[2];

					double value = pow(mu*(l1 - l2), 2.0) + pow(ga*(a1 - a2), 2.0) + pow(ga*(b1 - b2), 2.0);
					value = exp(-0.5 * value);

					triples.push_back(Triple(row_count,center_id,value));
					triples.push_back(Triple(row_count,neigh_id,-value));
					triples.push_back(Triple(row_count + 1, center_id + pixel_num, value));
					triples.push_back(Triple(row_count + 1, neigh_id + pixel_num, -value));
					triples.push_back(Triple(row_count + 2, center_id + 2 * pixel_num, value));
					triples.push_back(Triple(row_count + 2, neigh_id + 2 * pixel_num, -value));
					row_count = row_count + 3;
				}
			}
		}
	}
	SpMat sp_mat(row_count,3 * pixel_num);
	sp_mat.setFromTriplets(triples.begin(), triples.end());
	return sp_mat;
}