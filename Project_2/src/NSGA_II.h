//
// Created by Ingunn on 23.02.2018.
//

#ifndef PROJECT_2_NSGA_II_H
#define PROJECT_2_NSGA_II_H

#include "Genotype.h"
#include <vector>

using namespace std;


class NSGA_II {
    vector<Genotype> population;
    void fastNonDominatedSort();
    void crowdingDistanceAssignment();
    void crowdedComparison();
    void mainLoop();
};


#endif //PROJECT_2_NSGA_II_H
