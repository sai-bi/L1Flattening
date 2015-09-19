#ifndef OPTIMIZE_H
#define OPTIMIZE_H
#include <iostream>
#include <vector>
#include <cmath>
#include "reflectance-cluster.h"
#include <cstdio>
#include <fstream>
#include "shading.h"
using namespace std;
using namespace cv;

/*
 * Shrinkage process
 */
Mat_<double> Shrink(const Mat_<double>& input, double lambda){
    Mat_<double> output(input.rows, input.cols, 0.0);
    for (int i = 0; i < input.rows; i++){
        for (int j = 0; j < input.cols; j++){
            double temp = input(i, j);
            if (temp > lambda){
                output(i, j) = temp - lambda;
            }
            else if (temp < -lambda){
                output(i, j) = temp + lambda;
            }
        }
    }
    return output;
}


/*
 * Construct the reflectance sparseness matrix, shading smooth matrix
 * 
 * Input:
 *      adjacent: adjacent matrix for clusters.
 *      pairwise_weight: weight betwee a pair of clusters
 *      intensity: average intensity for each cluster
 *      cluster_same_region: whether two clusters belong to same region
 * Ouput:
 *      A, D, C: see energy function
 */
void GetReflectanceSparseMatrix(const Mat_<int>& adjacent, 
                                const Mat_<double>& pairwise_weight,
                                const Mat_<double>& intensity,
                                const Mat_<int>& cluster_same_region,
                                Mat_<double>& A, 
                                Mat_<double>& D, 
                                Mat_<double>& C){
    int cluster_num = pairwise_weight.rows;
    vector<Point2i> cluster_pairs;
    for (int i = 0; i < adjacent.rows; ++i){
        for (int j = i + 1; j < adjacent.cols; ++j){
            if (adjacent(i, j) == 0){
                continue;
            }
            cluster_pairs.push_back(Point2i(i,j));
        }
    }

    int pair_num = cluster_pairs.size();
    cout << "Number of adjacent pairs: " << pair_num << endl;
    cout << "Sum of adjacent pixels: " << sum(adjacent)[0] << endl;

    A = Mat_<double>(pair_num, cluster_num, 0.0);
    D = Mat_<double>(pair_num, cluster_num, 0.0);
    C = Mat_<double>(pair_num, 1, 0.0);
    double ratio = 1.0;

    for (int i = 0; i < cluster_pairs.size(); ++i){
        int label_1 = cluster_pairs[i].x;
        int label_2 = cluster_pairs[i].y;
        double w = pairwise_weight(label_1, label_2);
        A(i, label_1) = w;
        A(i, label_2) = -w;
        
        // here if two clusters are adjacent but belong to different regions, increase weight for 
        // shading smooth.
        if (cluster_same_region(label_1, label_2) == 0){
            ratio = 1.0;
        }
        else{
            ratio = 1.0;
        }
        D(i, label_1) = 1 * ratio;
        D(i, label_2) = -1 * ratio;
        C(i) = ratio * (intensity(label_2) - intensity(label_1));
    }
}


void GetSmoothMatrix(const Mat_<Vec3b>& image,
                     const Mat_<int>& label,
                     int cluster_num,
                     Mat_<double>& DTD,
                     Mat_<double>& DTC){
    int width = image.cols;
    int height = image.rows;
    vector<Point2i> pixel_pair_1;
    vector<Point2i> pixel_pair_2;

    for (int i = 0; i < height - 1; ++i){
        for (int j = 0; j < width - 1; ++j){
            if (label(i, j) != label(i, j + 1)){
                pixel_pair_1.push_back(Point2i(i, j));
                pixel_pair_2.push_back(Point2i(i, j + 1));
            }
            if (label(i, j) != label(i + 1, j)){
                pixel_pair_1.push_back(Point2i(i, j));
                pixel_pair_2.push_back(Point2i(i + 1, j));
            }
        }
    }
    int pair_num = pixel_pair_1.size();
    cout << "Number of adjacent pixel pairs: " << pair_num << endl;
    Mat_<double> pair_label(pair_num,2);
    Mat_<double> C(pair_num, 1, 0.0);
    for (int i = 0; i < pixel_pair_1.size(); ++i){
        Point2i p_1 = pixel_pair_1[i];
        Point2i p_2 = pixel_pair_2[i];
        int label_1 = label(p_1.x, p_1.y);
        int label_2 = label(p_2.x, p_2.y);
        pair_label(i, 0) = label_1;
        pair_label(i, 1) = label_2;
        double r_1 = image(p_1.x, p_1.y)[0];
        double g_1 = image(p_1.x, p_1.y)[1];
        double b_1 = image(p_1.x, p_1.y)[2];
        double r_2 = image(p_2.x, p_2.y)[0];
        double g_2 = image(p_2.x, p_2.y)[1];
        double b_2 = image(p_2.x, p_2.y)[2];
        C(i) = log((r_2 + g_2 + b_2 + 1) / (r_1 + g_1 + b_1 + 1));
    }
    cout << "Write pair matrix..." << endl;
    ofstream fout;
    fout.open("C:\\Users\\BiSai\\Documents\\GitHub\\Intrinsic\\IntrinsicImage\\matrix-pair.txt");
    fout << format(pair_label, "csv") << endl;
    fout.close();
    fout.open("C:\\Users\\BiSai\\Documents\\GitHub\\Intrinsic\\IntrinsicImage\\matrix-c2.txt");
    fout << format(C, "csv") << endl;
    fout.close();

    cout << "Read DTD and DTC..." << endl;
    DTD = Mat_<double>(cluster_num, cluster_num, 0.0); 
    DTC = Mat_<double>(cluster_num, 1.0, 0.0);
    ifstream fin;
    fin.open("C:\\Users\\BiSai\\Documents\\GitHub\\Intrinsic\\IntrinsicImage\\matrix-dtd.txt");
    for (int i = 0; i < cluster_num; ++i){
        for (int j = 0; j < cluster_num; ++j){
            fin >> DTD(i, j);
        }
    }
    fin.close();
    fin.open("C:\\Users\\BiSai\\Documents\\GitHub\\Intrinsic\\IntrinsicImage\\matrix-dtc.txt");
    for (int i = 0; i < cluster_num; ++i){
        fin >> DTC(i);
    }
    fin.close();
}


vector<map<int, double> > GetB(int cluster_num, double weight){
    vector<map<int, double> > columns(cluster_num);
    int count = 0;
    for (int i = 0; i < cluster_num; ++i){
        for (int j = i + 1; j < cluster_num; ++j){
            columns[i][count] = weight;
            columns[j][count] = -weight;
            count++;
        }
    };
    return columns;
}

/*
 * Get B^T * B;
 * Here we cannot represent the matrix directly, due to the size of the matrix.
 * This is a way to represent sparse matrix.
 */
Mat_<double> GetBTB(vector<map<int, double> >& columns){
	// ifstream fin;
	// fin.open("matrix-btb.txt");
    int cluster_num = columns.size();
	Mat_<double> BTB(cluster_num, cluster_num, 0.0);
	for (int i = 0; i < cluster_num; ++i){
        for (int j = 0; j < cluster_num; ++j){
            // fin >> BTB(i, j);
			if(i != j){
				BTB(i,j) = -1;
			}
			else{
				BTB(i,j) = cluster_num-1;
			}
        }
    }
    return BTB;

	
	/*
    for (int i = 0; i < cluster_num; ++i){
        for (int j = i; j < cluster_num; ++j){
            double temp = 0;
            for (auto itr = columns[i].begin(); itr != columns[i].end(); ++itr){
                int index = itr->first;
                if (columns[j].count(index) > 0){
                    temp += (itr->second * columns[j][index]);
                }
            }
            BTB(i, j) = temp;
            BTB(j, i) = temp;
        }
    }
    return BTB;
	*/
}

/*
 * Calculate B^T * v
 */
Mat_<double> MultiplyBT(vector<map<int, double> >& columns, 
                        const Mat_<double>& v){
    Mat_<double> result(columns.size(), 1, 0.0);
    
    for (int i = 0; i < columns.size(); ++i){
        double temp = 0;
        for (auto itr = columns[i].begin(); itr != columns[i].end(); ++itr){
            int index = itr->first;
            temp += (itr->second * v(index));
        }
        result(i) = temp;
    }
    return result;
}



/* 
 * Calculate the difference of elements in v.
 * v: v is a vector.
 */
Mat_<double> GlobalSparseMulti(double weight, 
                               const Mat_<double>& v){
    int row_num = (v.rows * (v.rows - 1)) / 2;
    Mat_<double> result(row_num, 1, 0.0);
    int count = 0;
    for (int i = 0; i < v.rows; ++i){
        for (int j = i + 1; j < v.rows; ++j){
            result(count) = weight * (v(i) - v(j));
            count++;
        }
    }
    return result;
}

/*
* Show the sparsity result, with each reflectance value indicated by a random color;
*/
void ShowSparsity(const Mat_<double>& reflectance,
                  const Mat_<int>& label){
    map<int, Vec3b> color_index;
    int width = label.cols;
    int height = label.rows;
    Mat_<Vec3b> image(height, width);
    RNG random_generator(getTickCount());
    for (int i = 0; i < reflectance.rows; ++i){
        int temp = exp(reflectance(i));
        if (color_index.count(temp) > 0){
            continue;
        }
        Vec3b color;
        color[0] = random_generator.uniform(0, 256);
        color[1] = random_generator.uniform(0, 256);
        color[2] = random_generator.uniform(0, 256);

        color_index[temp] = color;
    }

    for (int i = 0; i < height; ++i){
        for (int j = 0; j < width; ++j){
            int curr_label = label(i, j);
            double temp = exp(reflectance(curr_label));
            image(i, j) = color_index[temp];
        }
    }

    cout << "Number of reflectance values: " << color_index.size() << endl;
    imshow("Reflectance Sparsity", image);
}

void GetClusterIntensity(vector<ReflectanceCluster>& clusters,
                         const Mat_<Vec3b>& image,
                         Mat_<double>& intensity,
                         Mat_<double>& average_rgb){
    int cluster_num = clusters.size();
    intensity = Mat_<double>(cluster_num, 1, 0.0);
    average_rgb = Mat_<double>(cluster_num, 3, 0.0);
    for (int i = 0; i < clusters.size(); ++i){
        vector<Point2i> pixels_in_cluster = clusters[i].GetPixelLocations();
        double temp = 0;
        double r = 0;
        double g = 0;
        double b = 0;
        for (int j = 0; j < pixels_in_cluster.size(); ++j){
            int x = pixels_in_cluster[j].x;
            int y = pixels_in_cluster[j].y;
            temp = temp + (image(x, y)[0] + image(x, y)[1] + image(x, y)[2]) / 3.0;
            r += image(x, y)[2];
            g += image(x, y)[1];
            b += image(x, y)[0];
        }
        assert(pixels_in_cluster.size() != 0);
        temp = temp / pixels_in_cluster.size();
        r = r / pixels_in_cluster.size();
        g = g / pixels_in_cluster.size();
        b = b / pixels_in_cluster.size();
        intensity(i) = log(temp + 1);
        average_rgb(i, 0) = r;
        average_rgb(i, 1) = g;
        average_rgb(i, 2) = b;
    }
}

void OptimizeReflectance(vector<ReflectanceCluster>& clusters,
                         const Mat_<Vec3b>& image,
                         const Mat_<int>& pixel_label,
                         double alpha,
                         double beta,
                         double theta,
                         double lambda,
                         const Mat_<int>& region,
                         Mat_<double>& curr_r){
    int width = image.cols;
    int height = image.rows;
    int cluster_num = clusters.size();
    Mat_<double> intensity(cluster_num, 1, 0.0);
    Mat_<int> adjacent;
    Mat_<int> cluster_same_region;
    Mat_<double> pairwise_weight;
    Mat_<double> A;
    Mat_<double> BTB;
    Mat_<double> C;
    Mat_<double> D;
    Mat_<double> average_rgb;

    GetClusterIntensity(clusters, image, intensity, average_rgb);
    GetClusterAdjacentMatrix(pixel_label, cluster_num, region, adjacent, cluster_same_region);
    GetPairwiseWeight(clusters, image, cluster_same_region, pairwise_weight);
    // C and D is useless here
    GetReflectanceSparseMatrix(adjacent, pairwise_weight, intensity, cluster_same_region, A, D, C);

    A = alpha / 2.0 * A;

    cout << "Get matrix B & BTB..." << endl;
    vector<map<int, double> > B = GetB(cluster_num, beta / 2.0);
    BTB = GetBTB(B);
	BTB = (beta * beta / 4.0) * BTB;

    cout << "Calculate left hand matrix..." << endl;
    Mat_<double> identity_matrix = Mat_<double>::eye(cluster_num, cluster_num);
    Mat_<double> left_hand = theta * identity_matrix + lambda * A.t() * A + lambda * BTB;

    double stop_condition = 0.001;
    int a_row_num = A.rows;
    int b_row_num = (cluster_num * (cluster_num - 1)) / 2;
    int d_row_num = D.rows;
    Mat_<double> b_1(a_row_num, 1, 0.0);
    Mat_<double> b_2(b_row_num, 1, 0.0);
    Mat_<double> d_1(a_row_num, 1, 0.0);
    Mat_<double> d_2(b_row_num, 1, 0.0);

    // find the value that is closest to the average intensity
    double average_intensity = mean(intensity)[0];
    int closest_index = 0;
    double closest_value = 1e7;
    for (int i = 0; i < intensity.rows; ++i){
        double temp = abs(intensity(i) - average_intensity);
        if (temp < closest_value){
            closest_value = temp;
            closest_index = i;
        }
    }
    cout << "Average intensity: " << average_intensity << endl;
    double fix_r = average_intensity;
    // intensity(closest_index) = fix_r;
    Mat_<double> old_r(cluster_num, 1, 0.0);
    curr_r = intensity.clone();
    Mat_<double> final_r = intensity.clone();
    int iteration_num = 20;

    for (int i = 0; i < iteration_num; ++i){
        if (norm(curr_r - old_r) < stop_condition){
            curr_r = final_r.clone();
            break;
        }
        old_r = curr_r.clone();

		Mat_<double> right_hand = lambda * A.t() * (d_1 - b_1) + lambda * MultiplyBT(B, d_2 - b_2) + theta * intensity;
        solve(left_hand, right_hand, curr_r);

        // set hard constraint
        final_r = curr_r.clone();

        Mat_<double> temp_1 = A * curr_r;
        Mat_<double> temp_2 = GlobalSparseMulti(beta / 2.0, curr_r);
        d_1 = Shrink(temp_1 + b_1, 1.0 / lambda);
        d_2 = Shrink(temp_2 + b_2, 1.0 / lambda);
        b_1 = b_1 + temp_1 - d_1;
        b_2 = b_2 + temp_2 - d_2;

        double part_1 = sum(abs(temp_1))[0];
        double part_2 = sum(abs(temp_2))[0];
        Mat_<double> temp_3;
        pow(curr_r - intensity, 2.0, temp_3);
        double part_3 = sum(temp_3)[0];
        double total = part_1 + part_2 + theta / 2.0 * part_3;
        cout << "Iter: " << i << " " << total << " " << part_1 * (2.0 / alpha) << " " << part_2 * (2.0 / beta)
            << " " << part_3 << endl;
    }

    // show reflectance sparsify result
    // cout << "Show reflectance sparsity..." << endl;
    // ShowSparsity(curr_r, pixel_label);

    cout << "Get shading and reflectance image..." << endl;
    Mat_<Vec3b> reflectance;
    Mat_<double> shading;
    GetShadingReflectance(image, curr_r, clusters, reflectance, shading);
    cout << "Get maximum shading value..." << endl;
    double max_shading_value = *max_element(shading.begin(), shading.end());
    cout << "Max shading: " << max_shading_value << endl;
	shading = (200.0 / max_shading_value) * shading;
	// shading = shading  * 200; 

    imshow("Result", reflectance);
    imshow("Shading", Mat_<uchar>(shading));
    waitKey(0);
    return;
}

/*
 * Solve the l1-regularization problem using Bregman Iteration.
 * Return the reflectance value.
 * For details, please see "The Split Bregman Method for L1 Regularized Problems"
 * The Engergy function is: alpha / 2 * |AR| + beta /2 * |BR| + theta / 2 * (DR+C)^2
 * 
 * clusters: superpixels of the image
 * image: original image
 * pixel_label: a matrix with each pixel indicating its cluster index
 * region: region result got by Berkely EdgeSaliency
 * curr_r: output reflectance
 */
void SolveReflectance(vector<ReflectanceCluster>& clusters, 
                      const Mat_<Vec3b>& image, 
                      const Mat_<int>& pixel_label,
                      double alpha,
                      double beta,
                      double theta,
                      double lambda,
                      const Mat_<int>& region,
                      Mat_<double>& curr_r){
    int width = image.cols;
    int height = image.rows;
    int cluster_num = clusters.size();
    Mat_<double> intensity(cluster_num, 1, 0.0);
    Mat_<double> average_rgb(cluster_num, 3, 0.0);
    Mat_<int> adjacent;
    Mat_<int> cluster_same_region;
    Mat_<double> pairwise_weight;
    Mat_<double> A;
    Mat_<double> BTB;
    Mat_<double> C;
    Mat_<double> D;

    GetClusterIntensity(clusters, image, intensity, average_rgb);
    GetClusterAdjacentMatrix(pixel_label, cluster_num, region, adjacent, cluster_same_region);

    // Mat_<Vec3b> r_image = imread("C:\\Users\\BiSai\\Documents\\GitHub\\Intrinsic\\IntrinsicImage\\r.png");
    GetPairwiseWeight(clusters, image, cluster_same_region, pairwise_weight);
    GetReflectanceSparseMatrix(adjacent, pairwise_weight, intensity, cluster_same_region, A, D, C);
    cout << "Get D and C..." << endl;

    
    ofstream fout;
    fout.open("label.txt");
    fout << format(pixel_label, "csv") << endl;
    fout.close();
    exit(-1);
    fout.open("./matrix-a.txt");
    fout << format(A, "csv") << endl;
    fout.close();
    fout.open("./matrix-d.txt");
    fout << format(D, "csv") << endl;
    fout.close();
    fout.open("matrix-c.txt");
    fout << format(C.t(), "csv") << endl;
    fout.close();
    fout.open("intensity.txt");
    fout << format(intensity.t(), "csv") << endl;
    fout.close();
    fout.open("average_rgb.txt");
    fout << format(average_rgb.t(), "csv") << endl;
    fout.close();
    fout.open("pairweight.txt");
    fout << format(pairwise_weight, "csv") << endl;
    fout.close();
    fout.open("adjacent.txt");
    fout << format(adjacent, "csv") << endl;
    fout.close();
	
	// Here I just need the input for matlab
	exit(-1);

    A = alpha / 2.0 * A;
    D = theta / 2.0 * D;
    C = theta / 2.0 * C;

    Mat_<double> DTD = D.t() * D;
    Mat_<double> DTC;
    // GetSmoothMatrix(image, pixel_label, cluster_num, DTD, DTC);

    cout << "Get matrix B & BTB..." << endl;
    vector<map<int, double> > B = GetB(cluster_num, beta / 2.0);
    BTB = GetBTB(B);

    cout << "Calculate left hand matrix..." << endl;
    Mat_<double> left_hand = D.t() * D + lambda * A.t() * A + lambda * BTB;
    // Mat_<double> left_hand = theta * DTD + lambda * A.t() * A + lambda * BTB;

    double stop_condition = 0.001;
    int a_row_num = A.rows;
    int b_row_num = (cluster_num * (cluster_num - 1)) / 2;
    int d_row_num = D.rows;
    Mat_<double> b_1(a_row_num, 1, 0.0);
    Mat_<double> b_2(b_row_num, 1, 0.0);
    Mat_<double> d_1(a_row_num, 1, 0.0);
    Mat_<double> d_2(b_row_num, 1, 0.0);
    
    // find the value that is closest to the average intensity
    double average_intensity = mean(intensity)[0];
    int closest_index = 0;
    double closest_value = 1e10;
    for (int i = 0; i < intensity.rows; ++i){
        double temp = abs(intensity(i) - average_intensity);
        if (temp < closest_value){
            closest_value = temp;
            closest_index = i;
        }
    }
    cout << "Average intensity: " << average_intensity << endl;
    double fix_r = average_intensity;
    intensity(closest_index) = fix_r;

    Mat_<double> old_r(cluster_num, 1, 0.0);
    // solve(DTD, -1.0 * DTC, curr_r);
    // curr_r(closest_index) = fix_r;
    // Mat_<double> final_r = curr_r.clone();
    curr_r = intensity.clone();
    curr_r(closest_index) = fix_r;
    Mat_<double> final_r = intensity.clone();
    
    int iteration_num = 10;

    for (int i = 0; i < iteration_num; ++i){
        if (norm(curr_r - old_r) < stop_condition){
            curr_r = final_r.clone();
            break;
        }
        old_r = curr_r.clone();

        cout << "NAN number: " << sum(curr_r != curr_r)[0] << endl;

        Mat_<double> right_hand = lambda * A.t() * (d_1 - b_1) + lambda * MultiplyBT(B, d_2 - b_2) - D.t() * C;
        // Mat_<double> right_hand = lambda * A.t() * (d_1 - b_1) + lambda * MultiplyBT(B, d_2 - b_2) - theta * DTC;
        solve(left_hand, right_hand, curr_r);
        cout << "Min index: " << curr_r(closest_index) << endl;
        
        // set hard constraint
        final_r = curr_r.clone();
        curr_r(closest_index) = fix_r;

        Mat_<double> temp_1 = A * curr_r;
        Mat_<double> temp_2 = GlobalSparseMulti(beta / 2.0, curr_r);
        d_1 = Shrink(temp_1 + b_1, 1.0 / lambda);
        d_2 = Shrink(temp_2 + b_2, 1.0 / lambda);
        b_1 = b_1 + temp_1 - d_1;
        b_2 = b_2 + temp_2 - d_2;

        /*
        double part_1 = sum(abs(temp_1))[0];
        double part_2 = sum(abs(temp_2))[0];
        Mat_<double> temp_3;
        pow(D * curr_r + C, 2.0, temp_3);
        double part_3 = sum(temp_3)[0];
        double total = part_1 + part_2 + part_3;
        cout <<"Iter: "<<i<<" "<<total << " " << part_1 * (2.0 / alpha) << " " << part_2 * (2.0 / beta)
             << " " << part_3 * (4 / (theta * theta)) << endl;
        */
        cout << "Iter " << i << endl;
    }

    // show reflectance sparsify result
    cout << "Show reflectance sparsity..." << endl;
    ShowSparsity(curr_r, pixel_label);

    cout << "Get shading and reflectance image..." << endl;
    Mat_<Vec3b> reflectance;
    Mat_<double> shading;
    GetShadingReflectance(image, curr_r, clusters, reflectance, shading);
    double max_shading_value = *max_element(shading.begin(), shading.end());
    double min_shading_value = *min_element(shading.begin(), shading.end());
    cout << "Max shading: " << max_shading_value << endl;
    cout << "Min shading: " << min_shading_value << endl;
    shading = (255.0 / min_shading_value) * shading;

    imshow("Result", reflectance);
    imshow("Shading", Mat_<uchar>(shading));
    waitKey(0);
    return;
}
#endif


