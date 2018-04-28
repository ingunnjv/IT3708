#ifndef PROJECT_3_JSSP_H
#define PROJECT_3_JSSP_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <utility>

struct task{
    uint16_t task_id;
    uint16_t job_id;
    uint16_t machine_no;
    double process_time;
    task(){
        this->task_id = 0;
        this->job_id = 0;
        this->machine_no = 0;
        this->process_time = 0;
    };
    task(uint16_t task_id, uint16_t job_id, uint16_t machine_no, double process_time){
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

    task_matrix job_tasks;
    task_matrix machine_tasks;
    task source_task;
};


#endif //PROJECT_3_JSSP_H
