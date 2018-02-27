#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include "image_loader.h"
#include "genotype.h"
#include "nsga2.h"
#include "utils.h"
#include <opencv2/highgui.hpp>
#include <opencv2/core/eigen.hpp>
#include <set>

#pragma once



using namespace std;
using Eigen::MatrixXd;

int main(int argc, char *argv[]) {
    ImageLoader image = ImageLoader();
    image.LoadImagesFromFolder("353013");
    image.ExtractRGBChannels();

    Genotype g = Genotype(image.r_channel, image.g_channel, image.b_channel);

    //cout << "Rows of eigen image = " << red.rows() << endl;
    //cout << "Cols of eigen image = " << red.cols() << endl;

    //cv::imshow("Red", image.r_image);
    //cv::waitKey(0);

    double mutation_rate, crossover_Rate, tournament_size, time_limit, generation_limit, population_size;
    setUserArgs(argc, argv, mutation_rate, crossover_Rate, tournament_size,
                time_limit, generation_limit, population_size);
    Nsga2 ga = Nsga2(mutation_rate, crossover_Rate, tournament_size, time_limit, generation_limit);


    //cout << "START\n";
    //ga.primMST(image.r_channel, image.g_channel, image.b_channel);
    //cout << "END\n";


    MatrixXd m(2,2);
    m(0,0) = 3;
    m(1,0) = 2.5;
    m(0,1) = -1;
    m(1,1) = m(1,0) + m(0,1);
    std::cout << m << std::endl;
    return 0;
}
