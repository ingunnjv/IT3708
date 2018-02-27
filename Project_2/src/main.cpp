#include <iostream>
#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include "image_loader.h"
#include "genotype.h"
#include "nsga2.h"
#include <opencv2/highgui.hpp>
#include <opencv2/core/eigen.hpp>


using namespace std;
using Eigen::MatrixXd;

int main(int argc, char *argv[]) {
    ImageLoader test = ImageLoader();
    test.LoadImagesFromFolder("353013");
    test.ExtractRGBChannels();

    Genotype g = Genotype(test.r_channel, test.g_channel, test.b_channel);

    //cout << "Rows of eigen image = " << red.rows() << endl;
    //cout << "Cols of eigen image = " << red.cols() << endl;

    //cv::imshow("Red", test.r_image);
    //cv::waitKey(0);

    Nsga2 nsga = Nsga2();
    cout << "START\n";
    nsga.primMST(test.r_channel, test.g_channel, test.b_channel);
    cout << "END\n";

    MatrixXd m(2,2);
    m(0,0) = 3;
    m(1,0) = 2.5;
    m(0,1) = -1;
    m(1,1) = m(1,0) + m(0,1);
    std::cout << m << std::endl;
    return 0;
}
