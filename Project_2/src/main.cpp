#pragma once
#include <iostream>
#include <Eigen/Dense>
#include "image_loader.h"
#include "nsga2.h"

using namespace std;

int main(int argc, char *argv[]) {
    // TEST running python script
    string command = "python \"..\\src\\testScript.py\"";
    string args = "";
    command += args;
    system(command.c_str());

    // Set hyper parameters
    double mutation_rate, crossover_rate, time_limit;
    uint16_t tournament_size, generation_limit, population_size;
    int problem_num;
    setUserArgs(argc, argv, mutation_rate, crossover_rate, tournament_size,
                time_limit, generation_limit, population_size, problem_num);

    // Load the test image and the solutions
    ImageLoader image = ImageLoader();
    image.LoadImagesFromFolder(to_string(problem_num));
    image.ExtractRGBChannels();

    // Create GA
    Nsga2 ga = Nsga2(mutation_rate, crossover_rate, tournament_size,
                     time_limit, generation_limit, population_size);
    // Initialize a population
    vector<Genotype> initial_pop (population_size);
    printf("Initializating a population with %d individuals..\n", population_size);
    ga.initializePopulation(image.r_channel, image.g_channel, image.b_channel, initial_pop);


    // TEST SEGMENT DECODING


    // Run evolutionary process
    printf("Starting evolutionary process (NSGA-II algorithm)..\n");
    ga.runMainLoop(image.r_channel, image.g_channel, image.b_channel, initial_pop);

    // Test decoding and viz
    //ga.population[4].genotypeToPhenotypeDecoding(image.r_channel.rows(), image.r_channel.cols());
    //ga.population[4].visualize(image.b_channel, image.g_channel, image.r_channel, image.r_channel.rows(), image.r_channel.cols());


    printf("Exiting program\n");

    return 0;
}
