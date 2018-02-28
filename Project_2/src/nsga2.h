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
    double generation_limit;
    uint16_t population_size;



public:
    // Constructors
    Nsga2();
    Nsga2(double mutation_rate, double crossover_rate, double tournament_size, double time_limit,
          double generation_limit);

    void initializePopulation(const Eigen::MatrixXi &red, const Eigen::MatrixXi &green, const Eigen::MatrixXi &blue);
    std::vector<int> primMST(const Eigen::MatrixXi &red, const Eigen::MatrixXi &green, const Eigen::MatrixXi &blue);

    std::vector< std::vector<Genotype> > fastNonDominatedSort();
    void crowdingDistanceAssignment();
    void crowdedComparison();
    void runMainLoop();
};


#endif //PROJECT_2_NSGA_II_H
