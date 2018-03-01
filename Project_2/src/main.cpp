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
    ga.initializePopulation(image.r_channel, image.g_channel, image.b_channel);
    cout << "Start decoding" << endl;
    //printMST(parent, image.r_channel.rows() * image.r_channel.cols(), image.r_channel, image.g_channel, image.b_channel);
    ga.population[4].genotypeToPhenotypeDecoding(image.r_channel.rows(), image.r_channel.cols());
    cout << "End decoding" << endl;

    // Run the algorithm
    //ga.runMainLoop();


    return 0;
}
