#ifndef PROJECT_3_ACO_H
#define PROJECT_3_ACO_H

#include "jssp.h"
#include "graph.h"

class ACO {
private:
    JSSP* jssp;
    Node source;
    std::vector<std::vector<Node>> nodes;
    std::vector<Edge> edges;

    int swarm_size;
    double initial_pheromone;
public:
    ACO(JSSP &jssp, int swarm_size, double initial_pheromone);
    void createJobGraph();

};


#endif //PROJECT_3_ACO_H
