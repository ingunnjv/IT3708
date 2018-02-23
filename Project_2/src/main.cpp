#include <iostream>
#include <math.h>
#include <vector>
#include <Eigen/Dense>
#include "image_loader.h"

#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
using namespace std;
using Eigen::MatrixXd;

int main() {
    cout << "Hello, World!" << endl;

    cv::Mat image;
    image = cv::imread("../Test Images/1/s6.jpg", cv::IMREAD_COLOR);
    cv::imshow("Displayed image", image);
    cv::waitKey(0);

    ImageLoader test = ImageLoader();
    cout << test.B << endl;

    MatrixXd m(2,2);
    m(0,0) = 3;
    m(1,0) = 2.5;
    m(0,1) = -1;
    m(1,1) = m(1,0) + m(0,1);
    std::cout << m << std::endl;
    return 0;
}
