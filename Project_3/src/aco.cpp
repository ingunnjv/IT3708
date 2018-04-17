#include "aco.h"
#include <chrono>
#include <random>
#include <climits>
#include <cfloat>
#include <fstream>

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
                vector<pair<task*, task*>> state_transitions = getStateTransitions(tabu);

                /* Select next state according to equation 1 */
                vector<double> state_transistion_probs = getStateTransitionProbs(state_transitions, ants[k].decidability_rule);
                int next_edge_index = chooseNextState(state_transistion_probs);
                pair<task*, task*> next_edge = state_transitions[next_edge_index];

                /* Move ant to selected state */
                ants[k].path.push_back(next_edge);

                /* Save selected state in tabu */
                updateTabu(tabu, next_edge.second);
                int dummy = 0;
            }
        }

        int current_cycles_best_makespan = INT_MAX;
        int all_time_best_makespan = INT_MAX;
        // all_time_best_schedule;
        // cycles_best_schedule;
        vector<ant> elites;
        vector<ant> not_elites;
        vector <schedule> schedules(this->swarm_size);
        for (int k = 0; k < this->swarm_size; k++){
            /* Build schedule */
            schedules[k].ant = &ants[k];
            buildSchedule(schedules[k], k);
            saveScheduleAsCSV(schedules[k]);


            if (schedules[k].makespan < all_time_best_makespan){
                continue;
            }
            if (schedules[k].makespan == current_cycles_best_makespan){
                elites.push_back(ants[k]);
            }
            else if (schedules[k].makespan < current_cycles_best_makespan){
                elites.clear();
                elites.push_back(ants[k]);
            }

        }

        for (auto & schedule: schedules){
            // if schedule->ant in elites
                /* for (auto &edge: schedule->ant.path){
                 *    task* task_i = edge.first;
                 *    task* task_j = edge.second;
                 *    pheromone_accumulator[task_i->task_id][task_j->task_id] = Q/current_cycles_best_makespan * elites.size()
                 *
                 *
                 *
                 * */
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


        cycle++;
    }

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

std::vector<double> ACO::getStateTransitionProbs(vector<pair<task *, task *>> state_transitions, uint8_t decidability_rule) {
    vector<double> state_transitions_probs;
    double pheromone_and_process_time_sum = 0;
    for(auto &transition: state_transitions){
        double edge_pheromone = 0;
        double task_process_time = 0;
        task* task_i = transition.first;
        task* task_j = transition.second;
        edge_pheromone = pheromone_trails[task_i->task_id][task_j->task_id];
        if(decidability_rule == decidability_rules::SPT){
            task_process_time = 1/task_j->process_time;
        }
        else if(decidability_rule == decidability_rules::LPT){
            task_process_time = task_j->process_time;
        }
        pheromone_and_process_time_sum += pow(edge_pheromone, alpha)*pow(task_process_time, beta);
    }
    for(auto &transition: state_transitions){
        double edge_pheromone = 0;
        double task_process_time = 0;
        task* task_i = transition.first;
        task* task_j = transition.second;
        edge_pheromone = pheromone_trails[task_i->task_id][task_j->task_id];
        if(decidability_rule == decidability_rules::SPT){
            task_process_time = 1/task_j->process_time;
        }
        else if(decidability_rule == decidability_rules::LPT){
            task_process_time = task_j->process_time;
        }
        double state_transitions_prob = (pow(edge_pheromone, alpha)*pow(task_process_time, beta))/pheromone_and_process_time_sum;
        state_transitions_probs.push_back(state_transitions_prob);
    }
    return state_transitions_probs;
}

void ACO::addAntPheromoneContribution(vector<vector<double>> pheromone_accumulator){

}

void ACO::updatePheromoneTrails(vector<vector<double>> pheromone_accumulator){

}

void ACO::buildSchedule(schedule &schedule, int k) {
    schedule.filename = "Schedule_" + to_string(k);

    // Add tasks to machines in correct order
    schedule.machine_schedules.resize(jssp->getNumMachines());
    int edge_count = 0;
    for (auto &edge: schedule.ant->path){
        task* task = edge.second;
        schedule_block block;
        block.task = task;
        if (schedule.machine_schedules[task->machine_no].empty()){
            block.start_time = 0;
        }
        else{
            block.start_time = schedule.machine_schedules[task->machine_no].back().start_time + schedule.machine_schedules[task->machine_no].back().task->process_time;
        }
        block.path_index = edge_count;
        schedule.machine_schedules[task->machine_no].push_back(block);
        edge_count++;
    }

    // Avoid collisions between jobs in different machines
    for (int machine_no_i = 0; machine_no_i < jssp->getNumMachines(); machine_no_i++){
        int block_no_i = 0;
        for (auto &machine_schedule_block_i: schedule.machine_schedules[machine_no_i]){

            for (int machine_no_j = 0; machine_no_j < jssp->getNumMachines(); machine_no_j++){
                if (machine_no_i != machine_no_j){
                    int block_no_j = 0;
                    for (auto &machine_schedule_block_j: schedule.machine_schedules[machine_no_j]){

                        if (machine_schedule_block_i.task->job_id == machine_schedule_block_j.task->job_id){
                            task* task_i = machine_schedule_block_i.task;
                            double start_i = machine_schedule_block_i.start_time;
                            task* task_j = machine_schedule_block_j.task;
                            double start_j = machine_schedule_block_j.start_time;
                            if  (!(start_i + task_i->process_time <= start_j or start_j + task_j->process_time <= start_i)){
                                // collision
                                if (machine_schedule_block_i.path_index < machine_schedule_block_j.path_index){
                                    // shift start time of block j and subsequent blocks in machine j
                                    for (int block_no = block_no_j; block_no < schedule.machine_schedules[machine_no_j].size(); block_no++){
                                        schedule.machine_schedules[machine_no_j][block_no].start_time += ((start_i + task_i->process_time) - start_j);
                                    }
                                }
                                else{
                                    // shift start time of block i and subsequent blocks in machine i
                                    for (int block_no = block_no_i; block_no < schedule.machine_schedules[machine_no_i].size(); block_no++){
                                        schedule.machine_schedules[machine_no_i][block_no].start_time += ((start_j + task_j->process_time) - start_i);
                                    }
                                }
                            }
                        }
                        block_no_j++;
                    }
                }
            }
            block_no_i++;
        }
    }

    // Find makespan of schedule
    double earliest_finish = DBL_MAX;
    for (auto &machine: schedule.machine_schedules){
        schedule_block last_task = machine.back();
        if (last_task.start_time + last_task.task->process_time < earliest_finish){
            earliest_finish = last_task.start_time + last_task.task->process_time;
        }
    }
    schedule.makespan = earliest_finish;
}

void ACO::saveScheduleAsCSV(schedule &schedule) {
    ofstream file;
    file.open ("../solutions/" + schedule.filename + ".csv");
    file << "Task,Start,Finish,Job,Description\n";

    for (auto &machine_schedule: schedule.machine_schedules){
        for (auto &schedule_block: machine_schedule){
            file << "Machine " << schedule_block.task->machine_no << ",";
            file << schedule_block.start_time << "," << schedule_block.start_time + schedule_block.task->process_time << ",";
            file << "Job " << schedule_block.task->job_id << ",J" << schedule_block.task->job_id << "\n";
        }
    }
    file.close();
}



