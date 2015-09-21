#ifndef GET_REGION_H
#define GET_REGION_H
#include <opencv2/opencv.hpp>
#include <string>
#include <list>
#include <stack>
#include <vector>
#include <fstream>
#define VALID(x, y, w, h, r) (x >= 0 && y >= 0 && x < h && y < w && r(x,y) < 0) 
using namespace std;
using namespace cv;

/*
 * Given a boundary map, get the regions.
 */
Mat_<int> GetRegion(string region_map_path){
    Mat_<uchar> image = imread(region_map_path, CV_LOAD_IMAGE_GRAYSCALE);
    int width = image.cols;
    int height = image.rows;
    Mat_<int> region(height, width, -1);
    list<Point2i> unvisited_points;
    for (int i = 0; i < height; ++i){
        for (int j = 0; j < width; ++j){
            if (image(i, j) == 0){
                unvisited_points.push_back(Point2i(i, j));
            }
        }
    }

    Mat_<Vec3b> region_image(height, width, Vec3b(0, 0, 0));
    int region_index = 0;
    RNG random_g(getTickCount());
    while (!unvisited_points.empty()){
        Point2i curr = unvisited_points.front();
        unvisited_points.pop_front();
        if (region(curr.x, curr.y) >= 0){
            continue;
        }
        stack<Point2i> my_stack;
        my_stack.push(curr);
        Vec3b color;
        color(0) = random_g.uniform(0, 256);
        color(1) = random_g.uniform(0, 256);
        color(2) = random_g.uniform(0, 256);

        while (true){
            if (my_stack.empty()){
                break;
            }
            Point2i p = my_stack.top();
            my_stack.pop();
            region(p.x, p.y) = region_index;
            region_image(p.x, p.y) = color;
            if (image(p.x, p.y) > 10){
                // cout << (int)image(p.x,p.y) << endl;
                continue;
            }
            if (VALID(p.x + 1, p.y, width, height, region)){
                my_stack.push(Point2i(p.x + 1, p.y));
            }
            if (VALID(p.x - 1, p.y, width, height, region)){
                my_stack.push(Point2i(p.x - 1, p.y));
            }
            if (VALID(p.x, p.y + 1, width, height, region)){
                my_stack.push(Point2i(p.x, p.y + 1));
            }
            if (VALID(p.x, p.y - 1, width, height, region)){
                my_stack.push(Point2i(p.x, p.y - 1));
            }
        }
        region_index++;
    }
    imshow("Region", region_image);

	/*
    ofstream fout;
    fout.open("region.txt");
    fout << format(region, "csv");
    fout.close();
	*/
    return region;
}



#endif
