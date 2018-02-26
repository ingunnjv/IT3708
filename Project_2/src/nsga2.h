//
// Created by Ingunn on 23.02.2018.
//

#ifndef PROJECT_2_NSGA2_H
#define PROJECT_2_NSGA2_H

#pragma once
#include "genotype.h"
#include <vector>

using namespace std;


class Nsga2 {
    vector<Genotype> population;
    void fastNonDominatedSort();
    void crowdingDistanceAssignment();
    void crowdedComparison();
    void mainLoop();
};


#endif //PROJECT_2_NSGA_II_H
