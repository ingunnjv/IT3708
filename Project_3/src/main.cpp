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
    if(jssp.readInputData("test_3x3")) return 1;

    /* Create ant colony optimization object */
    int swarm_size = 10;
    int cycles = 1000;
    double alpha = 0.2;
    double beta = 0.8;
    double rho = 0.7;
    double initial_pheromone = 0.1;
    ACO aco = ACO(jssp, swarm_size, cycles, alpha, beta, rho, initial_pheromone);
    aco.printPheromoneTrailsTable();

    // Delete old solution files
    for (int k = 0; k < swarm_size; k++){
        string file = "../solutions/Schedule_" + to_string(k) + ".csv";
        const char * filename = file.c_str();
        if( remove( filename ) != 0 )
            perror( "Error deleting file" );
    }
    aco.runOptimization();

//    string solutionFile = "Schedule_9";
//    printf("Print Gantt chart of solution..\n");
//    string command = "python \"..\\src\\run_gantt.py\"";
//    string args = " " + solutionFile;
//    command += args;
//    system(command.c_str());

    t = clock() - t;
    printf("Time consumed: %d clicks (%f seconds).\n", int(t), ((float)int(t))/CLOCKS_PER_SEC);
    printf("Exiting program\n");
    return 0;
}
