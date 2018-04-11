#ifndef PROJECT_3_GRAPH_H
#define PROJECT_3_GRAPH_H

#include <vector>

class Node;
class Edge;

class Node{
private:
    int id;
    std::vector<Edge*> edges;
    int machine_no;
    int process_time;
public:
    Node();
    Node(int id, int machine_no, int process_time);
    int getMachineNo();
    int getProcessTime();
    int getId();
    void addEdge(Edge &edge);
};

class Edge{
private:
    int id;
    std::vector<Node*> nodes;
    double pheromone;
public:
    Edge();
    Edge(int id, double initial_pheromone);
    void setPheromone(double pheromone);
    double getPheromone();
    int getId();
    void addNode(Node &node);
    std::vector<Node*> getNodes();
};






#endif //PROJECT_3_GRAPH_H
