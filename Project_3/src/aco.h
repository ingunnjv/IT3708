#ifndef PROJECT_3_ACO_H
#define PROJECT_3_ACO_H

#include "jssp.h"
#include "graph.h"

class ACO {
private:
    JSSP* jssp;

    int swarm_size;
    double initial_pheromone;
    std::vector<std::vector<double>> pheromone_trails;
public:
    ACO(JSSP &jssp, int swarm_size, double initial_pheromone);
    void initializePheromoneTrails();

};


#endif //PROJECT_3_ACO_H
