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
#include "utils.h"

enum neighbouring_approaches {ONE_SWAP = 0, ONE_INSERT, TWO_SWAP, TWO_INSERT};

struct bee{
    schedule schedule;
    std::vector<int> operations_sequence;
    int sequence_age;
};

class ABC{
private:
    JSSP* jssp;
    std::vector<bee> employed_bees;
    std::vector<int> old_bees_indices;
    bee* super_amazing_bee;
    schedule best_schedule;

    int num_food_sources;
    int abandonment_limit;
    int cycles;
    double p_local_search;

    std::vector<int> neighbour_list;
    std::vector<int> winning_neighbour_list;
    int nl_length;

public:
    ABC(JSSP &jssp, int num_food_sources, int abandonment_limit, int cycles, int nl_length, double p_local_search);
    void initColony();
    void initOperationSequence(bee &colony_bee);
    void initNeighbourList();
    void refillNeighbourList();

    std::vector<std::pair<task*, task*>> decodeOperationsToPath(const bee &colony_bee);
    void employedBeePhase();
    void onlookerBeePhase();
    void scoutBeePhase();
    std::pair<bee, int> selfAdaptiveStrategy(bee &colony_bee);
    void localSearch(bee &original_bee, int approach);
    std::pair<int, int> binaryTournamentSelection(int size);
    double computeAverageMakespan();

    void oneInsertion(bee &colony_bee);
    void oneSwap(bee &colony_bee);
    void twoInsertions(bee &colony_bee);
    void twoSwaps(bee &colony_bee);

    void runOptimization();
};


#endif //PROJECT_3_ABC_H
