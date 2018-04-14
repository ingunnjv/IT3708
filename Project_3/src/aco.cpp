#include "aco.h"
#include <chrono>
#include <random>

using namespace std;

ACO::ACO(JSSP &jssp, int swarm_size, int cycles, double alpha, double beta, double rho, double initial_pheromone){
    this->jssp = &jssp;
    this->swarm_size = swarm_size;
    this->cycles = cycles;
    this->alpha = alpha;
    this->beta = beta;
    this->rho = rho;
    this->initial_pheromone = initial_pheromone;
    this->pheromone_trails.resize(jssp.getNumTasks()+1, vector<double>(jssp.getNumTasks()+1));
}

void ACO::initializePheromoneTrails(){
    // Set all elements to zero (including diagonals)
    for (int row = 0; row < this->pheromone_trails.size(); row++) {
        for (int col = 0; col < this->pheromone_trails[row].size(); col++) {
            pheromone_trails[row][col] = -1;
        }
    }
    // Set edges from source to first tasks
    for (int i = 1; i < jssp->getNumJobs() + 1; i++){
        pheromone_trails[0][i] = initial_pheromone;
        pheromone_trails[i][0] = initial_pheromone;
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

void ACO::printPheromoneTrailsTable(){
    for(const auto &trails: pheromone_trails){
        for(const auto &trail: trails){
            printf("%d\t", int(trail));
        }
        printf("\n");
    }
    printf("\n");
}

void ACO::runOptimization() {
    vector<vector<int>> tabu(jssp->job_tasks.size(), vector<int>(jssp->job_tasks[0].size()));
    vector<vector<double>> pheromone_accumulator(jssp->getNumTasks()+1, vector<double>(jssp->getNumTasks()+1));
    vector<ant> ants(this->swarm_size, ant(0));
    task source = task(0, 0, 0, 0);
    initializePheromoneTrails();

    int cycle = 0;
    while(cycle < cycles){
        setMatrixToZero(pheromone_accumulator);

        for (int k = 0; k < this->swarm_size; k++){
            setMatrixToZero(tabu);
            /* Define decidability rule for each ant k */
            // Either eta = process_time(LPT) or eta = 1/process_time(SPT)
            auto decidability_rule = uint8_t(k % 2);
            ants[k] = ant(decidability_rule);

            /* Random assignment of the first operation */
//            int start_task = k % jssp->getNumJobs();
//            ants[k].path.push_back(make_pair(&source, &jssp->job_tasks[start_task][0]));
//            updateTabu(tabu, &jssp->job_tasks[start_task][0]);

            while(!isTabuFull(tabu)){

                /* Determine the set of operations achievable from the current state */
                vector<pair<task*, task*>> possible_edges = getStateTransitions(tabu);

                /* Select next state according to equation 1 */
                vector<double> state_transistion_probs = getStateTransitionProbs(possible_edges);
                int next_edge_index = chooseNextState(state_transistion_probs);
                pair<task*, task*> next_edge = possible_edges[next_edge_index];

                /* Move ant to selected state */
                ants[k].path.push_back(next_edge);

                /* Save selected state in tabu */
                updateTabu(tabu, next_edge.second);
            }
        }
        cycle++;
    }


/* AFTER ANTS HAVE ACQUIRED PATHS:
 *
 * current_cycles_best = INF
 * for every ant
 *  decode path to schedule
 *  find schedule length
 *  save ants schedule
 *  if schedule solution better than all_time_best
 *      save it as best
 *      save ant as elite
 *  if schedule better (or equal) than current_cycles_best
 *      save ant as elite
 * for every ant not elite
 *  delta t += ant pheromone contribution to edges w/ eq 3
 * for every ant elite
 *  delta t += ant pheromone contribution to edges w/ eq 4
 * for every edge
 *  eq 2
 * */

}

template<typename T>
void ACO::setMatrixToZero(vector<vector<T>> &matrix){
    for (int i = 0; i < matrix.size(); i++){
        for (int j = 0; j < matrix[i].size(); j++){
            matrix[i][j] = 0;
        }
    }
}

bool ACO::isTabuFull(vector<vector<int>> &tabu){
    for (int i = 0; i < tabu.size(); i++){
        for (int j = 0; j < tabu[i].size(); j++){
            if (tabu[i][j] == 0){
                return false;
            }
        }
    }
    return true;
}


int ACO::chooseNextState(vector<double> &state_transistion_probs){
    unsigned seed = (unsigned) chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    uniform_real_distribution<double> rand_distribution(0.0, 1.0);

    double random = rand_distribution(generator);
    double probability_accumulator = 0.0;
    for (int state_i = 0; state_i < state_transistion_probs.size(); state_i++){
        probability_accumulator += state_transistion_probs[state_i];
        if(random <= probability_accumulator){
            return state_i;
        }
    }
    return state_transistion_probs.size() - 1;
}

void ACO::updateTabu(vector<vector<int>> &tabu, task* next_task){
    uint16_t tabu_row = next_task->job_id;
    for (int i = 0; i < jssp->job_tasks[tabu_row].size(); i++){
        if (jssp->job_tasks[tabu_row][i].task_id == next_task->task_id){
            tabu[tabu_row][i] = 1;
            return;
        }
    }
}

vector<pair<task *, task *>> ACO::getStateTransitions(const vector<vector<int>> &tabu){
    vector<pair<task *, task *>> possible_edges;
    vector<task*> fronts;

    for(int i = 0; i < tabu.size(); i++){
        for (int j = 0; j < tabu[i].size(); j++){
            if (tabu[i][j] == 0){
                if (j == 0)
                    possible_edges.push_back(make_pair(&jssp->source_task, &jssp->job_tasks[i][j]));
                else {
                    possible_edges.push_back(make_pair(&jssp->job_tasks[i][j - 1], &jssp->job_tasks[i][j]));
                    fronts.push_back(&jssp->job_tasks[i][j - 1]);
                }
                break;
            }
        }
    }
    for(auto &front: fronts) {
        for (int i = 0; i < jssp->machine_tasks[front->machine_no].size(); i++) {
            if (front->task_id != jssp->machine_tasks[front->machine_no][i].task_id) {
                int col = jssp->machine_tasks[front->machine_no][i].task_id % jssp->getNumMachines();
                int row = jssp->machine_tasks[front->machine_no][i].job_id;
                if (tabu[row][col - 1] == 1) {
                    possible_edges.push_back(make_pair(front, &jssp->machine_tasks[front->machine_no][i]));
                }
            }
        }
    }
    return possible_edges;

}

std::vector<double> ACO::getStateTransitionProbs(vector<pair<task *, task *>> possible_edges){

}

void ACO::addAntPheromoneContribution(vector<vector<double>> pheromone_accumulator){

}

void ACO::updatePheromoneTrails(vector<vector<double>> pheromone_accumulator){

}