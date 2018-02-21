#include <iostream>
#include <math.h>
#include <vector>
#include <Eigen/Dense>
#include "image_loader.h"
using namespace std;
using Eigen::MatrixXd;

int main() {
    cout << "Hello, World!" << endl;

    ImageLoader test = ImageLoader();
    cout << test.B;

    MatrixXd m(2,2);
    m(0,0) = 3;
    m(1,0) = 2.5;
    m(0,1) = -1;
    m(1,1) = m(1,0) + m(0,1);
    std::cout << m << std::endl;
    return 0;
}
