#include "aco.h"
#include <chrono>
#include <random>
#include <climits>
#include <cfloat>
#include <fstream>
#include "schedule_builder.h"
#include "utils.h"

using namespace std;

ACO::ACO(JSSP &jssp, int swarm_size, int cycles, double alpha, double beta, double rho, double initial_pheromone,
         double Q, double max_pheromone, double min_pheromone, double optimal_solution_val) {
    this->jssp = &jssp;
    this->swarm_size = swarm_size;
    this->cycles = cycles;
    this->alpha = alpha;
    this->beta = beta;
    this->rho = rho;
    this->Q = Q;
    this->initial_pheromone = initial_pheromone;
    this->max_pheromone_on_trails = max_pheromone;
    this->min_pheromone_on_trails = min_pheromone;
    this->optimal_solution_val = optimal_solution_val;
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
    for (int i = 1; i < this->pheromone_trails.size(); i += jssp->getNumJobs()){
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
    int i = 0;
    for(const auto &trails: pheromone_trails){
        printf("Task %d: ", i);
        for(const auto &trail: trails){
            printf("%.3f\t", trail);
        }
        printf("\n");
        i++;
    }
    printf("\n");
}

void ACO::runOptimization() {
    vector<vector<int>> tabu(jssp->job_tasks.size(), vector<int>(jssp->job_tasks[0].size()));
    vector<vector<double>> pheromone_accumulator(jssp->getNumTasks()+1, vector<double>(jssp->getNumTasks()+1));
    vector<ant> ants(this->swarm_size, ant(0));
    initializePheromoneTrails();

    schedule current_cycles_best_schedule;
    schedule all_time_best_schedule;
    double average_makespan = 0;
    all_time_best_schedule.makespan = DBL_MAX;

    int cycle = 0;
    while(cycle < cycles){
        if(cycle % 1 == 0 and cycle != 0){
            printf("Cycle: %d\n", cycle);
            printf("- Shortest makespan all time: %f\n", all_time_best_schedule.makespan);
            printf("- Average makespan size: %f\n", average_makespan);
        }
        setMatrixToZero(pheromone_accumulator);
        generatePaths(tabu, ants);

        current_cycles_best_schedule.makespan = DBL_MAX;
        average_makespan = 0;
        vector<int> elites;
        vector <schedule> schedules(this->swarm_size);
        for (int k = 0; k < this->swarm_size; k++){
            /* Build schedule */
            buildSchedule(schedules[k], ants[k].path, jssp);
            average_makespan += schedules[k].makespan;

            if (schedules[k].makespan < all_time_best_schedule.makespan){
                all_time_best_schedule = schedules[k];
            }
            if (schedules[k].makespan == current_cycles_best_schedule.makespan){
                elites.push_back(k);
            }
            else if (schedules[k].makespan < current_cycles_best_schedule.makespan){
                elites.clear();
                elites.push_back(k);
                current_cycles_best_schedule = schedules[k];
            }

        }

        for (int k = 0; k < this->swarm_size; k++){
            addAntPheromoneContribution(pheromone_accumulator, elites, ants[k], k, schedules[k].makespan);
        }
        updatePheromoneTrails(pheromone_accumulator);
        average_makespan = average_makespan / swarm_size;

        cycle++;
    }
    saveScheduleAsCSV(all_time_best_schedule, "Best-ant-schedule", jssp);
    printScoreToScreen(all_time_best_schedule.makespan, optimal_solution_val);
    callPythonGanttChartPlotter("Best-ant-schedule");
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
    vector<pair<task *, task *>> state_transitions;
    vector<task*> fronts;
    // Find edges to next available task within a job
    for(int i = 0; i < tabu.size(); i++){
        for (int j = 0; j < tabu[i].size(); j++){
            if (tabu[i][j] == 0){
                if (j == 0)
                    state_transitions.push_back(make_pair(&jssp->source_task, &jssp->job_tasks[i][j]));
                else {
                    state_transitions.push_back(make_pair(&jssp->job_tasks[i][j - 1], &jssp->job_tasks[i][j]));

                }
                break;
            }
            else if(tabu[i][j] == 1){
                fronts.push_back(&jssp->job_tasks[i][j]);
            }
        }
    }
    // Find edges to available tasks within machines
    for(auto &front: fronts) {
        for (int i = 0; i < jssp->machine_tasks[front->machine_no].size(); i++) {
            if (front->task_id != jssp->machine_tasks[front->machine_no][i].task_id) {
                int col = (jssp->machine_tasks[front->machine_no][i].task_id - 1) % jssp->getNumMachines();
                int row = jssp->machine_tasks[front->machine_no][i].job_id;
                if (tabu[row][col] == 0) {
                    if(col == 0){
                        state_transitions.push_back(make_pair(front, &jssp->machine_tasks[front->machine_no][i]));
                    }
                    else if(tabu[row][col - 1] == 1){
                        state_transitions.push_back(make_pair(front, &jssp->machine_tasks[front->machine_no][i]));
                    }
                }
            }
        }
    }
    return state_transitions;
}

vector<double> ACO::getStateTransitionProbs(vector<pair<task *, task *>> state_transitions, uint8_t decidability_rule,
                                            const vector<vector<int>> &tabu) {
    vector<double> state_transitions_probs;
    double pheromone_and_heuristic_sum = 0;

//    double max_edge_pheromone = 0;
//    double min_edge_pheromone = DBL_MAX;
//    double max_heuristic_val = 0;
//    double min_heuristic_val = DBL_MAX;

    vector<double> remaining_time_per_job;
    vector<int> remaining_tasks_per_job;

    srand(unsigned(time(0)));
    remaining_time_per_job.resize((unsigned int)jssp->getNumJobs());
    remaining_tasks_per_job.resize((unsigned int)jssp->getNumJobs());

    // Create counter to tasks added per jobb, remaining time per job, and remaining tasks per jobb
    for(int i = 0; i < jssp->getNumJobs(); i++){
        remaining_time_per_job[i] = 0;
        remaining_tasks_per_job[i] = jssp->getNumMachines();
        for(int j = 0; j < jssp->getNumMachines(); j++){
            if(tabu[i][j] != 1){
                remaining_time_per_job[i] += jssp->job_tasks[i][j].process_time;
            }
            if(tabu[i][j] == 1){
                remaining_tasks_per_job[i]--;
            }
        }
    }

//    /* Find minmax of pheromone edges and process time for normalization */
//    for(auto &transition: state_transitions){
//        double decidability_heuristic = 0;
//        task* task_i = transition.first;
//        task* task_j = transition.second;
//        double edge_pheromone = pheromone_trails[task_i->task_id][task_j->task_id];
//        if(decidability_rule == decidability_rules::MRTasks){
//            decidability_heuristic = remaining_tasks_per_job[task_j->job_id];
//        }
//        else if(decidability_rule == decidability_rules::MRTime){
//            decidability_heuristic = remaining_time_per_job[task_j->job_id];
//        }
//        max_edge_pheromone = max(max_edge_pheromone, edge_pheromone);
//        min_edge_pheromone = min(min_edge_pheromone, edge_pheromone);
//        max_heuristic_val = max(max_heuristic_val, decidability_heuristic);
//        min_heuristic_val = min(min_heuristic_val, decidability_heuristic);
//    }

    for(auto &transition: state_transitions){
        double decidability_heuristic = 0;
        task* task_i = transition.first;
        task* task_j = transition.second;
        double edge_pheromone = pheromone_trails[task_i->task_id][task_j->task_id];
        if(decidability_rule == decidability_rules::MRTasks){
            decidability_heuristic = remaining_tasks_per_job[task_j->job_id];
        }
        else if(decidability_rule == decidability_rules::MRTime){
            decidability_heuristic = remaining_time_per_job[task_j->job_id];
        }
        // Normalize
        //edge_pheromone = (edge_pheromone)/(max_edge_pheromone);
        //decidability_heuristic = (decidability_heuristic)/(max_heuristic_val);
        pheromone_and_heuristic_sum += pow(edge_pheromone, alpha)*pow(decidability_heuristic, beta);
    }
    for(auto &transition: state_transitions){
        double decidability_heuristic = 0;
        task* task_i = transition.first;
        task* task_j = transition.second;
        double edge_pheromone = pheromone_trails[task_i->task_id][task_j->task_id];
        if(decidability_rule == decidability_rules::MRTasks){
            decidability_heuristic = remaining_tasks_per_job[task_j->job_id];
        }
        else if(decidability_rule == decidability_rules::MRTime){
            decidability_heuristic = remaining_time_per_job[task_j->job_id];
        }

        //edge_pheromone = (edge_pheromone)/(max_edge_pheromone);
        //decidability_heuristic = (decidability_heuristic)/(max_heuristic_val);

        double state_transitions_prob = (pow(edge_pheromone, alpha)*pow(decidability_heuristic, beta))/pheromone_and_heuristic_sum;
        state_transitions_probs.push_back(state_transitions_prob);
    }
    return state_transitions_probs;
}

void ACO::addAntPheromoneContribution(vector<vector<double>> &pheromone_accumulator,
                                      vector<int> elites,
                                      const ant &colony_ant, const int ant_nr, double makespan) {
    if (isElementInVector(ant_nr, elites)){
        for (auto &edge: colony_ant.path){
            task* task_i= edge.first;
            task* task_j = edge.second;
            pheromone_accumulator[task_i->task_id][task_j->task_id] += (Q / makespan) * elites.size();
            pheromone_accumulator[task_j->task_id][task_i->task_id] += (Q / makespan) * elites.size();
        }
    }
    else{
        for (auto &edge: colony_ant.path){
            task* task_i= edge.first;
            task* task_j = edge.second;
            pheromone_accumulator[task_i->task_id][task_j->task_id] += Q / makespan;
            pheromone_accumulator[task_j->task_id][task_i->task_id] += Q / makespan;
        }
    }
}

void ACO::updatePheromoneTrails(const vector<vector<double>> &pheromone_accumulator){
    for (int row = 0; row < this->pheromone_trails.size(); row++) {
        for (int col = 0; col < this->pheromone_trails[row].size(); col++) {
            pheromone_trails[row][col] = rho * pheromone_trails[row][col] + pheromone_accumulator[row][col];
            if (pheromone_trails[row][col] > max_pheromone_on_trails){
                pheromone_trails[row][col] = max_pheromone_on_trails;
            }
            else if (pheromone_trails[row][col] < min_pheromone_on_trails){
                pheromone_trails[row][col] = min_pheromone_on_trails;
            }
        }
    }
}

void ACO::generatePaths(vector<vector<int>> &tabu, vector<ant> &ants) {
    for (int k = 0; k < this->swarm_size; k++){
        setMatrixToZero(tabu);
        /* Define decidability rule for each ant k */
        // Either eta = process_time(LPT) or eta = 1/process_time(SPT)
        auto decidability_rule = uint8_t(k % 2);
        ants[k] = ant(decidability_rule);

        while(!isTabuFull(tabu)){
            /* Determine the set of operations achievable from the current state */
            vector<pair<task*, task*>> state_transitions = getStateTransitions(tabu);

            /* Select next state according to equation 1 */
            vector<double> state_transistion_probs = getStateTransitionProbs(state_transitions, ants[k].decidability_rule, tabu);
            int next_edge_index = chooseNextState(state_transistion_probs);
            pair<task*, task*> next_edge = state_transitions[next_edge_index];

            /* Move ant to selected state */
            ants[k].path.push_back(next_edge);

            /* Save selected state in tabu */
            updateTabu(tabu, next_edge.second);
        }
    }
}




