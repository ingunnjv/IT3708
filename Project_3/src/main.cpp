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
    ACO aco = ACO(jssp, 0, 0);
    aco.initializePheromoneTrails();
    aco.printPheromoneTrailsTable();




    t = clock() - t;
    printf("Time consumed: %d clicks (%f seconds).\n", int(t), ((float)int(t))/CLOCKS_PER_SEC);
    printf("Exiting program\n");
    return 0;
}
