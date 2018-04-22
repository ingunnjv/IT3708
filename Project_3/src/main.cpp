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
    JSSP jssp = JSSP();
    if(jssp.readInputData("6")) return 1;

    /* Create parameters and ant colony optimization object */
    int swarm_size = jssp.getNumJobs()/2;
    int cycles = 1000;
    double alpha = 1;
    double beta = 5;
    double rho = 0.8;
    double initial_pheromone = 0.5;
    double max_pheromone = 1;
    double min_pheromone = 0.01;
    double Q = 100;
    ACO aco = ACO(jssp, swarm_size, cycles, alpha, beta, rho, initial_pheromone, Q,
    max_pheromone, min_pheromone);

    /* Run optimization */
    aco.runOptimization();
    aco.printPheromoneTrailsTable();

    /* Create parameters and artificial bee colony optimization object */
  //  int num_food_sources = 50;
//    int abandonment_limit = 100;
//    int cycles = 5000;
//    ABC abc = ABC(jssp, num_food_sources, abandonment_limit, cycles);



    /* Create gantt chart from python script */
    string solutionFile = "Best_bee";
    printf("Print Gantt chart of best solution..\n");
    string command = "python \"..\\src\\run_gantt.py\"";
    string args = " " + solutionFile;
    command += args;
    system(command.c_str());


    t = clock() - t;
    printf("Time consumed: %d clicks (%f seconds).\n", int(t), ((float)int(t))/CLOCKS_PER_SEC);
    printf("Exiting program\n");
    return 0;
}
