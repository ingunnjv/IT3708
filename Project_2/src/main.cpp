#pragma once
#include <iostream>
#include <Eigen/Dense>
#include "image_loader.h"
#include "nsga2.h"

using namespace std;

int main(int argc, char *argv[]) {
    // Set hyper parameters
    double mutation_rate, crossover_Rate, time_limit;
    uint16_t tournament_size, generation_limit, population_size;
    int problem_num;
    setUserArgs(argc, argv, mutation_rate, crossover_Rate, tournament_size,
                time_limit, generation_limit, population_size, problem_num);

    // Load the test image and the solutions
    ImageLoader image = ImageLoader();
    image.LoadImagesFromFolder(to_string(problem_num));
    image.ExtractRGBChannels();

    // Create GA
    Nsga2 ga = Nsga2(mutation_rate, crossover_Rate, tournament_size,
                     time_limit, generation_limit, population_size);
    // Initialize a population
    cout << "Initializating a population with " << population_size << " individuals..." << endl;
    ga.initializePopulation(image.r_channel, image.g_channel, image.b_channel);

    // Run evolutionary process
    ga.runMainLoop(image.r_channel, image.g_channel, image.b_channel);

    // Test decoding and viz
    //ga.population[4].genotypeToPhenotypeDecoding(image.r_channel.rows(), image.r_channel.cols());
    //ga.population[4].visualize(image.b_channel, image.g_channel, image.r_channel, image.r_channel.rows(), image.r_channel.cols());

    //

    // Run the algorithm
    //ga.runMainLoop();
    //cv::namedWindow( "Test image", cv::WINDOW_AUTOSIZE );
    //cv::imshow( "Test image", image.test_image );
    //cv::waitKey(0);

    return 0;
}
