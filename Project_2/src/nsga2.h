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
    double time_limit;
    uint16_t tournament_size;
    uint16_t generation_limit;
    uint16_t population_size;



public:
    // Constructors
    Nsga2();
    Nsga2(double mutation_rate, double crossover_rate, uint16_t tournament_size, double time_limit,
          uint16_t generation_limit, uint16_t population_size);

    void primMST(const Eigen::MatrixXi &red, const Eigen::MatrixXi &green, const Eigen::MatrixXi &blue);

    std::vector< std::vector<Genotype> > fastNonDominatedSort();
    std::tuple<double, double> objectiveValueSort(std::vector<Genotype> &genotypes, uint8_t objective_num);
    void crowdingDistanceAssignment(std::vector<Genotype> &front);
    Genotype crowdedComparison(const Genotype &gt1, const Genotype &gt2);
    void runMainLoop();
    std::vector<Genotype> makeNewPop(std::vector<Genotype> parent_pop);
};


#endif //PROJECT_2_NSGA_II_H
