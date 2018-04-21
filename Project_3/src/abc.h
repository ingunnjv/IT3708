#ifndef PROJECT_3_ABC_H
#define PROJECT_3_ABC_H

#pragma once
#include "jssp.h"
#include <chrono>
#include <random>
#include <algorithm>
#include <ctime>
#include <cstdlib>

enum initType{};

struct bee{
    std::vector<int> operations_sequence;
    double makespan;
    int sequence_age;
};

class ABC{
private:
    JSSP* jssp;
    std::vector<bee> employed_bees;
    int num_food_sources;
    int abandonment_limit;
    int cycles;

public:
    ABC(JSSP &jssp, int num_food_sources, int abandonment_limit, int cycles);
    void initColony();
    void initOperationSequence(bee &colony_bee);
};


#endif //PROJECT_3_ABC_H
