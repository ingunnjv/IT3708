#include "aco.h"

using namespace std;

ACO::ACO(JSSP &jssp, int swarm_size, double initial_pheromone){
    this->jssp = &jssp;
    this->swarm_size = swarm_size;
    this->initial_pheromone = initial_pheromone;
    this->pheromone_trails.resize(jssp.getNumTasks(), vector<double>(jssp.getNumTasks()));
}

void ACO::initializePheromoneTrails(){
    // Set all elements to zero (including diagonals)
    for (int row = 0; row < jssp->getNumTasks(); row++) {
        for (int col = 0; col < jssp->getNumTasks(); col++) {
            pheromone_trails[row][col] = -1;
        }
    }
    // Set unidirectional edges
    for (int row = 0; row < jssp->job_tasks.size(); row++) {
        for (int i = 0; i < jssp->job_tasks[row].size() - 1; i++) {
            int task_id_row = jssp->job_tasks[row][i].task_id;
            int task_id_col = jssp->job_tasks[row][i+1].task_id;
            pheromone_trails[task_id_row][task_id_col] = initial_pheromone;
            pheromone_trails[task_id_col][task_id_row] = -1;
        }
    }
    // Set bidirectional edges
    for (int machine = 0; machine < jssp->machine_tasks.size(); machine++){
        for (int task_i = 0; task_i < jssp->machine_tasks[machine].size() - 1; task_i ++){
            for (int task_j = task_i + 1; task_j < jssp->machine_tasks[machine].size(); task_j ++){
                int task_id_i = jssp->machine_tasks[machine][task_i].task_id;
                int task_id_j = jssp->machine_tasks[machine][task_j].task_id;
                pheromone_trails[task_id_i][task_id_j] = initial_pheromone;
                pheromone_trails[task_id_j][task_id_i] = initial_pheromone;
            }
        }
    }
}
