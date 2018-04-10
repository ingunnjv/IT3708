#pragma once
#include <iostream>
#include <Eigen/Dense>
#include <time.h>
#include "jssp.h"
using namespace std;

int main(int argc, char *argv[]) {
    clock_t t = clock();

    JSSP jssp = JSSP();
    if(jssp.readInputData("1")) return 1;

    t = clock() - t;
    printf("Time consumed: %d clicks (%f seconds).\n", int(t), ((float)int(t))/CLOCKS_PER_SEC);
    printf("Exiting program\n");

    return 0;
}
