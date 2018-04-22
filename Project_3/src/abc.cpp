#include <cfloat>
#include "abc.h"

using namespace std;

ABC::ABC(JSSP &jssp, int food_sources, int abandonment_limit, int cycles){
    this->jssp = &jssp;
    this->num_food_sources = food_sources;
    this->abandonment_limit = abandonment_limit;
    this->cycles = cycles;
    this->employed_bees.resize(num_food_sources);
    initColony();
}

void ABC::initColony() {
    double worst_makespan = 0;
    auto best_makespan = DBL_MAX;
    int i = 0;
    for (auto &bee: employed_bees){
        bee.sequence_age = 0;
        initOperationSequence(bee);

        vector<pair<task *, task *>> path = decodeOperationsToPath(bee);
        buildSchedule(bee.schedule, path, jssp);
        //printf("Bee %d: \tMakespan %f\n",i, bee.schedule.makespan);

        if (bee.schedule.makespan > worst_makespan){
            worst_makespan = bee.schedule.makespan;
            idiet_loser_bee = &employed_bees[i];
        }
        if (bee.schedule.makespan < best_makespan){
            best_makespan = bee.schedule.makespan;
            super_amazing_bee = &employed_bees[i];
        }
        i++;
    }
    saveScheduleAsCSV(super_amazing_bee->schedule, "Best_bee", jssp);
}

void ABC::initOperationSequence(bee &colony_bee) {
    unsigned seed = (unsigned) chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    uniform_real_distribution<double> task_choice_rule_rand_distribution(0.0, 1.0);

    vector<double> remaining_time_per_job;
    vector<int> remaining_tasks_per_job;
    vector<int> added_tasks_per_job;

    int num_tasks_remaining = jssp->getNumTasks();

    srand(unsigned(time(0)));
    remaining_time_per_job.resize((unsigned int)jssp->getNumJobs());
    remaining_tasks_per_job.resize((unsigned int)jssp->getNumJobs());
    added_tasks_per_job.resize((unsigned int)jssp->getNumJobs());

    // Create counter to tasks added per jobb, remaining time per job, and remaining tasks per jobb
    for(int i = 0; i < jssp->getNumJobs(); i++){
        remaining_time_per_job[i] = 0;
        remaining_tasks_per_job[i] = jssp->getNumMachines();
        added_tasks_per_job[i] = 0;
        for(int j = 0; j < jssp->getNumMachines(); j++){
            remaining_time_per_job[i] += jssp->job_tasks[i][j].process_time;
        }
    }

    // Choose next operation until there are no more tasks left
    while(num_tasks_remaining > 0){
        double r = task_choice_rule_rand_distribution(generator);
        int index = 0;
        if(r < 0.2){
            vector<int> job_choices;
            for(int job = 0; job < jssp->getNumJobs(); job++){
                if(remaining_tasks_per_job[job] != 0){
                    job_choices.push_back(job);
                }
            }
            if(job_choices.size() == 1){
                index = job_choices[0];
            }
            else{
                uniform_int_distribution<int> task_choice_rand_distribution(0, job_choices.size() - 1);
                index = job_choices[task_choice_rand_distribution(generator)];
            }

        }
        else if(r < 0.6 and r >= 0.2){
            // Choose next task based on time remaining for each job
            vector<double>::iterator max_element_it;
            max_element_it = max_element(remaining_time_per_job.begin(), remaining_time_per_job.end());
            index = (int)distance(remaining_time_per_job.begin(), max_element_it);
        }
        else if(r >= 0.6){
            // Choose next task based on tasks remaining for each job
            vector<int> job_choices;
            int maximum_remaining_tasks = 0;
            for(int i = 0; i < remaining_tasks_per_job.size(); i++){
                if(remaining_tasks_per_job[i] == maximum_remaining_tasks){
                    maximum_remaining_tasks = remaining_tasks_per_job[i];
                    job_choices.push_back(i);
                }
                else if(remaining_tasks_per_job[i] > maximum_remaining_tasks){
                    maximum_remaining_tasks = remaining_tasks_per_job[i];
                    job_choices.clear();
                    job_choices.push_back(i);
                }
            }
            if(job_choices.size() == 1){
                index = job_choices[0];
            }
            else{
                uniform_int_distribution<int> task_choice_rand_distribution(0, job_choices.size() - 1);
                index = job_choices[task_choice_rand_distribution(generator)];
            }
        }
        colony_bee.operations_sequence.push_back(index);
        remaining_time_per_job[index] -= jssp->job_tasks[index][added_tasks_per_job[index]].process_time;
        added_tasks_per_job[index]++;
        remaining_tasks_per_job[index]--;
        num_tasks_remaining--;
    }
}

vector<pair<task *, task *>> ABC::decodeOperationsToPath(const bee &colony_bee) {
    vector<pair<task *, task *>> path;
    vector<int> added_tasks_per_job;
    added_tasks_per_job.resize((unsigned)jssp->getNumJobs());

    for (auto &operation: colony_bee.operations_sequence){
        int task_nr = added_tasks_per_job[operation];
        task* task = &jssp->job_tasks[operation][task_nr];
        path.push_back(make_pair(&jssp->source_task, task));

        added_tasks_per_job[operation]++;
    }
    return path;
}

void ABC::employedBeePhase() {
    for (int i =0; i < num_food_sources; i++){
        /* Produce new solution according to self-adaptive strategy */
        /* Perform local search on new solution if r < P_L */
        /* Replace bee with new solution if it is better */
    }
}

void ABC::onlookerBeePhase() {
    for (int i =0; i < num_food_sources; i++) {
        /* Select food source according tournament selection */
        /* Produce new solution according to self-adaptive strategy */
        /* Update population if the new solution is better than or equal to the selected one */
    }
}

void ABC::scoutBeePhase() {
    /* If a solution in the population has not been improved during the last limit
     * number of trials , abandon it and put a new solution generated by performing
     * several insert operators to the best solution into the population. If there are
     * several such solutions, randomly select one. */
}

void ABC::runOptimization() {
    int cycle = 0;
    while(cycle < cycles) {

        employedBeePhase();
        onlookerBeePhase();
        scoutBeePhase();

        /* Memorize the best solution achieved so far */

        cycle++;
    }

}









