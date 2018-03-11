#pragma once
#include <iostream>
#include <Eigen/Dense>
#include "image_loader.h"
#include "nsga2.h"
#include "segmentation.h"
#include <time.h>

using namespace std;

int main(int argc, char *argv[]) {
    clock_t t = clock();

    // Set hyper parameters
    double mutation_rate, crossover_rate, time_limit;
    uint16_t tournament_size, generation_limit, population_size;
    int problem_num;
    setUserArgs(argc, argv, mutation_rate, crossover_rate, tournament_size,
                time_limit, generation_limit, population_size, problem_num);

    // Load the test image and the solutions
    ImageLoader image = ImageLoader();
    image.loadImagesFromFolder(to_string(problem_num));
    //image.segmentation(to_string(problem_num));
    image.extractRGBChannels();

    // Create GA
    Nsga2 ga = Nsga2(mutation_rate, crossover_rate, tournament_size,
                     time_limit, generation_limit, population_size);

    // Initialize a population
    vector<Genotype> initial_pop (population_size);
    printf("Initializating a population with %d individuals..\n", population_size);
    ga.initializePopulationFromMst(image.r_channel, image.g_channel, image.b_channel, initial_pop);

    // Run evolutionary process
    printf("Starting evolutionary process (NSGA-II algorithm)..\n");
    ga.runMainLoop(image.r_channel, image.g_channel, image.b_channel, initial_pop, image.test_image);

    // Run MPRI score check
    printf("Run MPRI scoring on the solutions..\n");
    string command = "python \"..\\src\\run.py\"";
    string args = " " + to_string(problem_num);
    command += args;
    system(command.c_str());

    t = clock() - t;
    printf ("Time consumed: %d clicks (%f seconds).\n", int(t), ((float)int(t))/CLOCKS_PER_SEC);
    printf("Exiting program\n");

    return 0;
}
