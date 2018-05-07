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
    this->acceptable_solution_makespan = optimal_solution_val;
    this->pheromone_trails.resize(jssp.getNumTasks()+1, vector<double>(jssp.getNumTasks()+1));
}

void ACO::initializePheromoneTrails(){
    // Set all elements to zero (including diagonals)
    for (int row = 0; row < this->pheromone_trails.size(); row++) {
        for (int col = 0; col < this->pheromone_trails[row].size(); col++) {
            pheromone_trails[row][col] = initial_pheromone;
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
        // Print progress to screen
        if(cycle % 1 == 0 and cycle != 0){
            printf("Cycle: %d\n", cycle);
            printf("- Shortest makespan all time: %f\n", all_time_best_schedule.makespan);
            printf("- Average makespan size: %f\n", average_makespan);
        }

        // Early stopping
        if(100.0*(all_time_best_schedule.makespan/acceptable_solution_makespan - 1) <= 10.0){
            break;
        }

        // Reset pheromone accumulator
        setMatrixToZero(pheromone_accumulator);

        // Create an ant-walk (a path) for each ant
        for (int k = 0; k < this->swarm_size; k++){
            setMatrixToZero(tabu);
            /* Define decidability rule for each ant k */
            auto decidability_rule = uint8_t(k % NUMBER_OF_RULES);
            ants[k] = ant(decidability_rule);

            /* Set a random job as the start task for the ant */
            unsigned seed = (unsigned) chrono::system_clock::now().time_since_epoch().count();
            default_random_engine generator(seed);
            uniform_int_distribution<int> rand_distribution(0, jssp->getNumJobs() - 1);
            double start_job = rand_distribution(generator);
            pair<task*, task*> initial_edge = make_pair(&jssp->source_task, &jssp->job_tasks[start_job][0]);
            tabu[start_job][0] = 1;
            ants[k].path.push_back(initial_edge);

            /* Generate edges until all tasks are included */
            while(!isTabuFull(tabu)){
                /* Determine the set of operations achievable from the current state */
                vector<pair<task*, task*>> state_transitions = getStateTransitions2(tabu, ants[k]);

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

        current_cycles_best_schedule.makespan = DBL_MAX;
        average_makespan = 0;
        vector<int> elites;
        vector <schedule> schedules(this->swarm_size);
        for (int k = 0; k < this->swarm_size; k++){
            /* Build schedule */
            buildSchedule(schedules[k], ants[k].path, jssp);
            average_makespan += schedules[k].makespan;

            // Compare solution quality to current best and all time best
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

        // Add each ants walk as contribution to the pheromone accumulator
        for (int k = 0; k < this->swarm_size; k++){
            addAntPheromoneContribution(pheromone_accumulator, elites, ants[k], k, schedules[k].makespan);
        }

        // Apply evaporation of pheromone and add the pheromone accumulator to the pheromone trails
        updatePheromoneTrails(pheromone_accumulator);
        average_makespan = average_makespan / swarm_size;

        cycle++;
    }
    saveScheduleAsCSV(all_time_best_schedule, "Best-ant-schedule", jssp);
    printScoreToScreen(all_time_best_schedule.makespan, acceptable_solution_makespan);
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

vector<pair<task *, task *>> ACO::getStateTransitions(const vector<vector<int>> &tabu) {
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

vector<pair<task *, task *>> ACO::getStateTransitions2(const std::vector<std::vector<int>> &tabu, const ant &colony_ant){
    vector<pair<task *, task *>> state_transitions;
    task* previous_task;
    if(colony_ant.path.empty())
        previous_task = &jssp->source_task;
    else
        previous_task = colony_ant.path.back().second;

    // Add transitions to all available jobs
    for(int i = 0; i < tabu.size(); i++){
        for (int j = 0; j < tabu[i].size(); j++){
            if (tabu[i][j] == 0){
                state_transitions.push_back(make_pair(previous_task, &jssp->job_tasks[i][j]));
                break;
            }
        }
    }
    return state_transitions;
};


vector<double> ACO::getStateTransitionProbs(vector<pair<task *, task *>> state_transitions, uint8_t decidability_rule,
                                            const vector<vector<int>> &tabu) {
    vector<double> state_transitions_probs;
    double pheromone_and_heuristic_sum = 0;

    vector<double> remaining_time_per_job;
    vector<int> remaining_tasks_per_job;

    srand(unsigned(time(0)));
    remaining_time_per_job.resize((unsigned int)jssp->getNumJobs());
    remaining_tasks_per_job.resize((unsigned int)jssp->getNumJobs());

    // Create counter to remaining time per job and remaining tasks per jobb
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
    for(auto &transition: state_transitions){
        task* task_i = transition.first;
        task* task_j = transition.second;
        double edge_pheromone = pheromone_trails[task_i->task_id][task_j->task_id];
        double decidability_heuristic = decodeHeuristicToValue(remaining_time_per_job, remaining_tasks_per_job, decidability_rule, task_j);
        pheromone_and_heuristic_sum += pow(edge_pheromone, alpha)*pow(decidability_heuristic, beta);
    }
    for(auto &transition: state_transitions){
        task* task_i = transition.first;
        task* task_j = transition.second;
        double edge_pheromone = pheromone_trails[task_i->task_id][task_j->task_id];
        double decidability_heuristic = decodeHeuristicToValue(remaining_time_per_job, remaining_tasks_per_job, decidability_rule, task_j);
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

double ACO::decodeHeuristicToValue(vector<double> &remaining_time_per_job, vector<int> &remaining_tasks_per_job,
                                    int decidability_rule, task* next_task) {
    if(decidability_rule == decidability_rules::MRTasks){
        return remaining_tasks_per_job[next_task->job_id];
    }
    else if(decidability_rule == decidability_rules::MRTime){
        return remaining_time_per_job[next_task->job_id];
    }
    else if(decidability_rule == decidability_rules::LRTasks){
        return 1/remaining_tasks_per_job[next_task->job_id];
    }
    else if(decidability_rule == decidability_rules::LRTime){
        return 1/remaining_time_per_job[next_task->job_id];
    }
}

void ACO::simpleTabuViz(vector<vector<int>> matrix){
    printf("Current tabu matrix: \n");
    for(const auto &items: matrix){
        for(const auto &item: items){
            printf("%d ", item);
        }
        printf("\n");
    }
}

void ACO::printAvailableStateTransitions(vector<pair<task*, task*>> availableStateTransitions) {
    printf("Available state transitions: \n");
    for(auto &transition: availableStateTransitions){
        int first_id = transition.first->task_id;
        int second_id = transition.second->task_id;
        printf("Task %i -> Task %i \n", first_id, second_id);
    }
}
