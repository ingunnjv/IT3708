#include "aco.h"

using namespace std;

ACO::ACO(JSSP &jssp, int swarm_size, double initial_pheromone){
    this->jssp = &jssp;
    this->swarm_size = swarm_size;
    this->initial_pheromone = initial_pheromone;
    this->source = Node(0, -1, -1);
}

void ACO::createJobGraph(){
    // Create nodes
    int id = 1;
    for (const auto &jobs: jssp->job_matrix){
        vector<Node> job_nodes;
        for (const auto &task: jobs){
            Node node = Node(id, task.first, task.second);
            job_nodes.push_back(node);
            id++;
        }
        nodes.push_back(job_nodes);
    }

    // Create machine tasks matrix
    vector<vector<Node*>> machine_tasks;
    machine_tasks.resize(jssp->getNumMachines());
    for (auto &jobs: nodes){
        for (auto &task: jobs){
            machine_tasks[task.getMachineNo()].push_back(&task);
        }
    }


    id = 0;
    for (int job = 0; job < jssp->getNumJobs(); job++){
        Edge edge = Edge(id, initial_pheromone);
        edge.addNode(nodes[job][0]);
        edges.push_back(edge);
        source.addEdge(edge);
        id++;
    }


    // Create unidirectional edges
    for (int job = 0; job < nodes.size(); job++){
        for (int task = 0; task < nodes[job].size(); task++){
            if (task < nodes[job].size()) {
                Edge edge = Edge(id, initial_pheromone);
                edge.addNode(nodes[job][task + 1]);
                edges.push_back(edge);
                nodes[job][task].addEdge(edges.back());
                id++;
            }
        }
    }
    // Create bidirectional edges
    for (int machine = 0; machine < machine_tasks.size(); machine++){
        for (int task_i = 0; task_i < machine_tasks[machine].size() - 1; task_i++){
            for (int task_j = task_i + 1; task_j < machine_tasks[machine].size(); task_j++){
                Edge edge = Edge(id, initial_pheromone);
                edges.push_back(edge);
                edges.back().addNode(*machine_tasks[machine][task_i]);
                edges.back().addNode(*machine_tasks[machine][task_j]);

                machine_tasks[machine][task_i]->addEdge(edges.back());
                machine_tasks[machine][task_j]->addEdge(edges.back());
                id++;
            }
        }
    }
}