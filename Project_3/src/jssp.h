#ifndef PROJECT_3_JSSP_H
#define PROJECT_3_JSSP_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <utility>

struct task{
    int task_id;
    int job_id;
    int machine_no;
    int process_time;
    task(int task_id, int job_id, int machine_no, int process_time){
        this->task_id = task_id;
        this->job_id = job_id;
        this->machine_no = machine_no;
        this->process_time = process_time;
    }
};
typedef std::vector<std::vector<task>> task_matrix;


class JSSP{
private:
    int num_jobs;
    int num_tasks;
    int num_machines;
public:
    JSSP();

    int readInputData(std::string problem_no);
    int getNumJobs() { return this->num_jobs; }
    int getNumTasks() { return this->num_tasks; }
    int getNumMachines() { return this->num_machines; }

    task_matrix jobs;
    task_matrix machine_tasks;
};


#endif //PROJECT_3_JSSP_H
