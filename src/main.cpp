#include "get-cluster.h"
#include <fstream>
using namespace cv;
using namespace std;

int main(int argc, char* argv[]){
    if (argc < 3){
        cout << "Usage: " << endl
            << "IntrinsicImage.exe [-i input_image_path] [-r region_file_path] [-m mask_file_path] [-o output_path]" << endl;
        exit(-1);
    }
    string image_path; 
    string region_file_path;
    string mask_file_path;
    string output_path;
    for (int i = 1; i < argc; i = i + 2){
        if (strcmp(argv[i], "-i") == 0){
            image_path = string(argv[i+1]);
            cout << "image_path: " << image_path << endl;
        }
        else if (strcmp(argv[i], "-r") == 0){
            region_file_path = string(argv[i + 1]);
            cout << "region_file_path: " << region_file_path << endl;
        }
        else if (strcmp(argv[i], "-m") == 0){
            mask_file_path = string(argv[i + 1]);
        }
        else if (strcmp(argv[i], "-o") == 0){
            output_path = string(argv[i + 1]);
            cout << "output_path: " << output_path << endl;
        }
    }

    Mat_<Vec3b> image = imread(image_path);
    int image_width = image.cols;
    int image_height = image.rows;
    Mat_<Vec3b> lab_image(image_height, image_width);    
    Mat_<int> region(image_height, image_width, 1);
    int expected_cluster_num = 500;
    Mat_<int> mask(image_height, image_width, 1);

    // transform into lab color space
    cvtColor(image, lab_image, CV_BGR2Lab);

    if(region_file_path.empty() == false){
        ifstream fin;
        fin.open(region_file_path);
        for (int i = 0; i < image_height; i++){
            for (int j = 0; j < image_width; j++){
                fin >> region(i, j);
            }
        }
        expected_cluster_num = 3000;
    }

    if(mask_file_path.empty() == false){
        Mat_<uchar> temp = imread(mask_file_path);
        mask = temp > 200;
    }

    Mat_<Vec3d> project_image = Mat_<Vec3d>(lab_image);
    double l_ratio = 0.3;
    for (int i = 0; i < image_height; i++){
        for (int j = 0; j < image_width; j++){
            project_image(i, j)[0] *= 0.3;
        }
    }

    double sigma = 0;
    double c = 500; // 1.0
    int min_size = 5;
    int cluster_num = 0;
    Mat_<Vec3b> output;
    Mat_<int> label;
    // vector<ReflectanceCluster> clusters;

    cout << "Segment the image..." << endl;
    GetReflectanceCluster(project_image, sigma, c, min_size, &cluster_num,
            output, label, expected_cluster_num, region);
    cout << "Super-pixel num: " << cluster_num << endl;

	
    ofstream fout;
    fout.open(output_path);
    fout << cv::format(label, cv::Formatter::FMT_CSV) << endl;
	
    return 0;
}
