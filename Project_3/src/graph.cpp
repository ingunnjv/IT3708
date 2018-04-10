#include "graph.h"

using namespace std;

/* Node implementation */
Node::Node(int id, int machine_no, int process_time){
    this->id = id;
    this->machine_no = machine_no;
    this->process_time = process_time;
}
int Node::getMachineNo(){
    return this->machine_no;
}
int Node::getProcessTime(){
    return this->process_time;
}
void Node::addEdge(Edge &edge){
    this->edges.push_back(&edge);
}
int Node::getId(){
    return this->id;
}


/* Edge implementation */
Edge::Edge(int id, double initial_pheromone){
    this->id = id;
    this->pheromone = initial_pheromone;
}
void Edge::setPheromone(double pheromone){
    this->pheromone = pheromone;
}
double Edge::getPheromone(){
    return this->pheromone;
}
void Edge::addNode(Node &node){
    this->nodes.push_back(&node);
}
std::vector<Node*> Edge::getNodes(){
    return this->nodes;
}
int Edge::getId(){
    return this->id;
}