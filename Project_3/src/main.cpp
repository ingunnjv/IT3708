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
    ACO aco = ACO(jssp, swarm_size, cycles, 0.2, 0.8, 0.7, 0.1);
    aco.runOptimization();
    aco.printPheromoneTrailsTable();

//    string solutionFile = "test_gantt";
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
