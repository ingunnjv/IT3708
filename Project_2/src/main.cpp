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

// A utility function to print the constructed MST stored in parent[]
void printMST(vector<int> parent, int num_pixels, const Eigen::MatrixXi &red, const Eigen::MatrixXi &green, const Eigen::MatrixXi &blue)
{
    pixel_t x, y;
    auto cols = uint16_t(red.cols());
    printf("Edge   Weight\n");
    for (int i = 1; i < num_pixels; i++) {
        x.row = i / cols;
        x.col = i % cols;
        y.row = parent[i] / cols;
        y.col = parent[i] % cols;
        printf("%d - %d    %f \n", parent[i], i, rgbDistance(y, x, red, green, blue));
    }
}

int main(int argc, char *argv[]) {
    ImageLoader image = ImageLoader();
    image.LoadImagesFromFolder("353013");
    image.ExtractRGBChannels();

    //Genotype g = Genotype(image.r_channel, image.g_channel, image.b_channel);

    //cout << "Rows of eigen image = " << red.rows() << endl;
    //cout << "Cols of eigen image = " << red.cols() << endl;

    //cv::imshow("Red", image.r_image);
    //cv::waitKey(0);

    double mutation_rate, crossover_Rate, time_limit;
    uint16_t tournament_size, generation_limit, population_size;
    setUserArgs(argc, argv, mutation_rate, crossover_Rate, tournament_size,
                time_limit, generation_limit, population_size);
    Nsga2 ga = Nsga2(mutation_rate, crossover_Rate, tournament_size,
                     time_limit, generation_limit, population_size);

    cout << "START\n";
    ga.initializePopulation(image.r_channel, image.g_channel, image.b_channel);
    cout << "END\n";
    //vector<int> parent_graph = ga.primMST(image.r_channel, image.g_channel, image.b_channel);
    //
    //printMST(parent_graph, image.r_channel.rows() * image.r_channel.cols(), image.r_channel, image.g_channel, image.b_channel);


    MatrixXd m(2,2);
    m(0,0) = 3;
    m(1,0) = 2.5;
    m(0,1) = -1;
    m(1,1) = m(1,0) + m(0,1);
    std::cout << m << std::endl;
    return 0;
}
