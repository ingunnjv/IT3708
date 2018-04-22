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
    for (auto &bee: employed_bees){
        bee.sequence_age = 0;
        initOperationSequence(bee);
    }
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
            // Choose next task RANDOMLY
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
            // Choose next task based on TIME remaining for each job
            vector<double>::iterator max_element_it;
            max_element_it = max_element(remaining_time_per_job.begin(), remaining_time_per_job.end());
            index = (int)distance(remaining_time_per_job.begin(), max_element_it);
        }
        else if(r >= 0.6){
            // Choose next task based on TASKS remaining for each job
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









