#pragma once
#include <iostream>
#include <Eigen/Dense>
#include <time.h>
#include "jssp.h"
#include "aco.h"

using namespace std;

int main(int argc, char *argv[]) {
    clock_t t = clock();

    JSSP jssp = JSSP();
    if(jssp.readInputData("test_2x2")) return 1;
    ACO aco = ACO(jssp, 0, 0);
    aco.initializePheromoneTrails();

    t = clock() - t;
    printf("Time consumed: %d clicks (%f seconds).\n", int(t), ((float)int(t))/CLOCKS_PER_SEC);
    printf("Exiting program\n");

    return 0;
}
