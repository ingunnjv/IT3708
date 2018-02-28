#ifndef PROJECT_2_NSGA2_H
#define PROJECT_2_NSGA2_H

#include "genotype.h"
#include <vector>
#include <Eigen/Dense>

class Nsga2 {
private:
    std::vector<Genotype> population;

    /* Hyper parameters */
    double mutation_rate;
    double crossover_rate;
    double time_limit;
    uint16_t tournament_size;
    uint16_t generation_limit;
    uint16_t population_size;

public:
    /* Constructors */
    Nsga2();
    Nsga2(double mutation_rate, double crossover_rate, uint16_t tournament_size, double time_limit,
          uint16_t generation_limit, uint16_t population_size);

    /* Initializes a population using a MST as the basis for gene encoding*/
    void initializePopulation(const Eigen::MatrixXi &red, const Eigen::MatrixXi &green, const Eigen::MatrixXi &blue);

    /* Creates a MST based on the distance in RGB space of a picture as basis for edges */
    std::vector<int> primMST(const Eigen::MatrixXi &red, const Eigen::MatrixXi &green, const Eigen::MatrixXi &blue);

    /* Sorts the population based non-domination fronts*/
    std::vector< std::vector<Genotype> > fastNonDominatedSort();

    /* Sorts a vector of genotypes based on the objective number */
    std::tuple<double, double> objectiveValueSort(std::vector<Genotype> &genotypes, uint8_t objective_num);

    /* Sorts a front of genotypes based on their crowding distance (best to worst -> greatest to lowest) */
    void crowdedDistanceSort(std::vector<Genotype> front);

    /* Assigns every genotype in a front a crowding_distance */
    void crowdingDistanceAssignment(std::vector<Genotype> &front);

    /* Compares two genotypes returning the best based on rank and crowding measure*/
    Genotype crowdedComparison(const Genotype &gt1, const Genotype &gt2);

    /* Creates a new offspring population using crossover and mutation */
    std::vector<Genotype> makeNewPop(std::vector<Genotype> parent_pop);

    /* Run the main loop of the algorithm */
    void runMainLoop();
};


#endif //PROJECT_2_NSGA_II_H
