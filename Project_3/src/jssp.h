#ifndef PROJECT_3_JSSP_H
#define PROJECT_3_JSSP_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <utility>

typedef std::vector<std::vector<std::pair<int, int>>> int_pair_matrix;

class JSSP{
private:
    int num_jobs;
    int num_machines;


public:
    JSSP();

    int readInputData(std::string problem_no);
    int getNumJobs() { return this->num_jobs; }
    int getNumMachines() { return this->num_machines; }

    int_pair_matrix job_matrix;
};


#endif //PROJECT_3_JSSP_H
