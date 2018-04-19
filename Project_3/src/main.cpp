#pragma once
#include <iostream>
#include <Eigen/Dense>
#include <time.h>
#include "jssp.h"
#include "aco.h"

using namespace std;

int main(int argc, char *argv[]) {
    clock_t t = clock();


    /* Create Job shop scheduling problem instance */
    JSSP jssp = JSSP();
    if(jssp.readInputData("1")) return 1;

    /* Create ant colony optimization object */
    int swarm_size = 10;
    int cycles = 10000;
    double alpha = 0.2;
    double beta = 0.8;
    double rho = 0.7;
    double initial_pheromone = 0.1;
    double max_pheromone = 1;
    double min_pheromone = 0.01;
    double Q = 5;
    ACO aco = ACO(jssp, swarm_size, cycles, alpha, beta, rho, initial_pheromone, Q,
    max_pheromone, min_pheromone);


    // Delete old solution files
//    for (int c = 0; c < cycles; c++){
//        string file = "../solutions/Cycle_" + to_string(c) + ".csv";
//        const char * filename = file.c_str();
//        if( remove( filename ) != 0 )
//            perror( "Error deleting file" );
//    }
    aco.runOptimization();
    aco.printPheromoneTrailsTable();

//    for (int c = 0; c < cycles; c++){
//        string solutionFile = "Cycle_" + to_string(c);
//        printf("Print Gantt chart of cycle %d solution..\n", c);
//        string command = "python \"..\\src\\run_gantt.py\"";
//        string args = " " + solutionFile;
//        command += args;
//        system(command.c_str());
//    }
    string solutionFile = "Best";
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
