#ifndef PROJECT_3_ACO_H
#define PROJECT_3_ACO_H

#include "jssp.h"
#include "graph.h"

class ACO {
private:
    JSSP* jssp;
    int swarm_size; // number of ants
    double alpha; // influence weight of pheromone
    double beta; // influence weight of heuristic
    double rho; // evaporation rate of pheromone
    double cycles; // iterations of the algorithm
    double initial_pheromone; // initial pheromone for all edges
    std::vector<std::vector<double>> pheromone_trails; // pheromone on all edges

public:
    ACO(JSSP &jssp, int swarm_size, double initial_pheromone);
    void initializePheromoneTrails();
    void printPheromoneTrailsTable();
    std::vector<task*> getPossibleStates(std::vector<std::vector<int>> visited_states);
    std::vector<double> getStateTransitionProbs(std::vector<task*> possibleStates);
    void addAntPheromoneContribution(std::vector<std::vector<double>> pheromone_accumulator);
    void updatePheromoneTrails(std::vector<std::vector<double>> pheromone_accumulator);

    void runOptimization();

};


#endif //PROJECT_3_ACO_H
