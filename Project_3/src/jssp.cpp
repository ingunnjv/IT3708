#include "jssp.h"
using namespace std;

JSSP::JSSP(){

}

int JSSP::readInputData(std::string problem_no){
    string line;
    ifstream file ("../data/" + problem_no + ".txt");
    int machine_no, process_time;
    if (file.is_open()){
        getline(file, line);
        stringstream stream(line);
        stream >> this->num_jobs;
        stream >> this->num_machines;
        this->job_tasks.resize(this->num_jobs);
        int task_id = 0;
        int job_id = 0;
        while(getline(file,line)){
            stringstream stream(line);
            while(1){
                stream >> machine_no;
                stream >> process_time;
                if(!stream) {
                    job_id++;
                    break;
                }
                task task_desc = task(task_id, job_id, machine_no, process_time);
                this->job_tasks[job_id].push_back(task_desc);
                task_id++;
            }
        }
        this->num_tasks = task_id;
        file.close();

        // Create machine tasks matrix
        machine_tasks.resize(num_machines);
        for (auto &jobs: this->job_tasks){
            for (auto &task: jobs){
                machine_tasks[task.machine_no].push_back(task);
            }
        }

        return 0;
    }
    else cout << "Unable to open file\n";
    return 1;
}