#pragma once
#include <iostream>
#include <time.h>
#include "jssp.h"
#include "aco.h"
#include "abc.h"

using namespace std;

int main(int argc, char *argv[]) {
    clock_t t = clock();

    /* Create Job shop scheduling problem instance */
    string optimizer = "aco";
    JSSP jssp = JSSP();
    if(jssp.readInputData("5")) return 1;
    double optimal_solution_val = 1451;

    if(optimizer == "aco"){
        /* Create parameters and ant colony optimization object */
        int swarm_size = 20;
        int cycles = 1000;
        double alpha = 10;
        double beta = 20;
        double rho = 0.95;
        double initial_pheromone = 10;
        double max_pheromone = 100;
        double min_pheromone = 1;
        double Q = 2000;
        ACO aco = ACO(jssp, swarm_size, cycles, alpha, beta, rho,
                      initial_pheromone, Q, max_pheromone, min_pheromone, optimal_solution_val);

        /* Run optimization */
        aco.runOptimization();
    }
    else if(optimizer == "abc"){
        /* Create parameters and artificial bee colony optimization object */
        int num_food_sources = 20;
        int abandonment_limit = 100;
        int cycles = 500;
        int nl_length = 200;
        double p_local_search = 0.1;
        ABC abc = ABC(jssp, num_food_sources, abandonment_limit, cycles,
                      nl_length, p_local_search, optimal_solution_val);
        /* Run optimization */
        abc.runOptimization();
    }

    t = clock() - t;
    printf("Time consumed: %d clicks (%f seconds).\n", int(t), ((float)int(t))/CLOCKS_PER_SEC);
    printf("Exiting program\n");
    return 0;
}
