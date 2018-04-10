#ifndef PROJECT_2_JSSP_H
#define PROJECT_2_JSSP_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <utility>

typedef std::vector<std::vector<std::pair<int, int>>> int_pair_matrix;

class JSSP{
private:
    int jobs;
    int machines;
    int_pair_matrix job_matrix;

public:
    JSSP();

    int readInputData(std::string problem_no);
};


#endif //PROJECT_2_JSSP_H
