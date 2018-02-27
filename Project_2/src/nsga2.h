//
// Created by Ingunn on 23.02.2018.
//

#ifndef PROJECT_2_NSGA2_H
#define PROJECT_2_NSGA2_H

#include "genotype.h"
#include <vector>
#include <Eigen/Dense>

class Nsga2 {
private:
    std::vector<Genotype> population;
    double mutation_rate;
    double crossover_rate;
    double tournament_size;
    double time_limit;
    uint16_t generation_limit;



public:
    // Constructors
    Nsga2();
    Nsga2(double mutation_rate, double crossover_rate, double tournament_size, double time_limit,
          uint16_t generation_limit);

    void primMST(Eigen::MatrixXi &red, Eigen::MatrixXi &green, Eigen::MatrixXi &blue);

    void fastNonDominatedSort();
    void crowdingDistanceAssignment();
    void crowdedComparison();
    void runMainLoop();
};


#endif //PROJECT_2_NSGA_II_H
