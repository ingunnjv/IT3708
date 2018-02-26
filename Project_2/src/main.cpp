#include <iostream>
#include <math.h>
#include <vector>
#include <Eigen/Dense>
#include "image_loader.h"
#include <opencv2/highgui.hpp>
#include <opencv2/core/eigen.hpp>


using namespace std;
using Eigen::MatrixXd;

int main() {
    ImageLoader test = ImageLoader();
    test.LoadImagesFromFolder("353013");
    test.ExtractRGBChannels();


    MatrixXd Red;
    cv::cv2eigen(test.r_image, Red); // convert from cv::Mat to Eigen::MatrixXd
    cout << "Rows of eigen image = " << Red.rows() << endl;
    cout << "Cols of eigen image = " << Red.cols() << endl;

    cv::imshow("Red", test.r_image);
    cv::waitKey(0);


    MatrixXd m(2,2);
    m(0,0) = 3;
    m(1,0) = 2.5;
    m(0,1) = -1;
    m(1,1) = m(1,0) + m(0,1);
    std::cout << m << std::endl;
    return 0;
}
