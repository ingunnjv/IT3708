//
// Created by Ingunn on 23.02.2018.
//

#ifndef PROJECT_2_NSGA2_H
#define PROJECT_2_NSGA2_H



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

    void primMST(const Eigen::MatrixXi &red, const Eigen::MatrixXi &green, const Eigen::MatrixXi &blue);
    uint32_t minKey(double key[], bool mstSet[], uint32_t num_pixels);

    std::vector< std::vector<Genotype> > fastNonDominatedSort();
    void crowdingDistanceAssignment();
    void crowdedComparison();
    void runMainLoop();
};


#endif //PROJECT_2_NSGA_II_H
