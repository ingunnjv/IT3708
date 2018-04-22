#ifndef PROJECT_3_ABC_H
#define PROJECT_3_ABC_H

#pragma once
#include "jssp.h"
#include <chrono>
#include <random>
#include <algorithm>
#include <ctime>
#include <cstdlib>
#include "schedule_builder.h"

struct bee{
    schedule schedule;
    std::vector<int> operations_sequence;
    int sequence_age;
};

class ABC{
private:
    JSSP* jssp;
    std::vector<bee> employed_bees;
    bee* idiet_loser_bee;
    bee* super_amazing_bee;
    int num_food_sources;
    int abandonment_limit;
    int cycles;
    // maybe save abandoned_bees?

public:
    ABC(JSSP &jssp, int num_food_sources, int abandonment_limit, int cycles);
    void initColony();
    void initOperationSequence(bee &colony_bee);
    std::vector<std::pair<task*, task*>> decodeOperationsToPath(const bee &colony_bee);
    void employedBeePhase();
    void onlookerBeePhase();
    void scoutBeePhase();

    void runOptimization();
};


#endif //PROJECT_3_ABC_H
