#ifndef PROJECT_2_ACO_H
#define PROJECT_2_ACO_H

#include "jssp.h"
#include "graph.h"

class ACO {
private:
    JSSP* jssp;
    Node* source;
    std::vector<Node*> nodes;
    std::vector<Edge*> edges;

    int swarm_size;
public:
    ACO(JSSP &jssp, int swarm_size);
    void createJobGraph();

};


#endif //PROJECT_2_ACO_H
