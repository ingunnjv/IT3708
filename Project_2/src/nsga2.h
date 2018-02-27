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
    void primMST(Eigen::MatrixXi red, Eigen::MatrixXi green, Eigen::MatrixXi blue);
    uint32_t minKey(double key[], bool mstSet[], uint32_t num_pixels);

    void fastNonDominatedSort();
    void crowdingDistanceAssignment();
    void crowdedComparison();
    void mainLoop();

};


#endif //PROJECT_2_NSGA_II_H
