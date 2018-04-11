#include "jssp.h"
using namespace std;

JSSP::JSSP(){

}

int JSSP::readInputData(std::string problem_no){
    string line;
    ifstream file ("../data/" + problem_no + ".txt");
    int machine_no, process_time;
    pair<int,int> task_desc;
    if (file.is_open()){
        getline(file, line);
        stringstream stream(line);
        stream >> this->num_jobs;
        stream >> this->num_machines;
        int job_i = 0;
        this->job_matrix.resize(this->num_jobs);
        while(getline(file,line)){
            stringstream stream(line);
            while(1){
                stream >> machine_no;
                stream >> process_time;
                if(!stream) {
                    job_i++;
                    break;
                }
                task_desc = make_pair(machine_no, process_time);
                this->job_matrix[job_i].push_back(task_desc);
            }
        }
        file.close();
        return 0;
    }
    else cout << "Unable to open file\n";
    return 1;
}