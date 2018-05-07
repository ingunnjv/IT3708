#ifndef PROJECT_3_ACO_H
#define PROJECT_3_ACO_H

#pragma once
#include "jssp.h"
#include "schedule_builder.h"

const int NUMBER_OF_RULES = 4;

enum decidability_rules{MRTasks = 0, MRTime, LRTasks, LRTime};

struct ant{
    std::vector<std::pair<task*,task*>> path;
    uint8_t decidability_rule;
    ant(uint8_t rule){
        this->decidability_rule = rule;
    }
};



class ACO {
private:
    JSSP* jssp;
    int swarm_size; // number of ants
    int cycles; // iterations of the algorithm
    double alpha; // influence weight of pheromone
    double beta; // influence weight of heuristic
    double rho; // evaporation rate of pheromone
    double Q; // constant factor in ant pheromone contribution
    double initial_pheromone; // initial pheromone for all edges
    double max_pheromone_on_trails;
    double min_pheromone_on_trails;
    double acceptable_solution_makespan;
    std::vector<std::vector<double>> pheromone_trails; // pheromone on all edges

public:
    ACO(JSSP &jssp, int swarm_size, int cycles, double alpha, double beta, double rho, double initial_pheromone,
            double Q, double max_pheromone, double min_pheromone, double optimal_solution_val);
    void initializePheromoneTrails();
    void printPheromoneTrailsTable();
    std::vector<std::pair<task *, task *>> getStateTransitions(const std::vector<std::vector<int>> &tabu);
    std::vector<std::pair<task *, task *>> getStateTransitions2(const std::vector<std::vector<int>> &tabu, const ant &colony_ant);
    std::vector<double> getStateTransitionProbs(std::vector<std::pair<task *, task *>> state_transitions, uint8_t decidability_rule,
                                                    const std::vector<std::vector<int>> &tabu);
    void addAntPheromoneContribution(std::vector<std::vector<double>> &pheromone_accumulator,
                                     std::vector<int> elites,
                                     const ant &colony_ant, const int ant_nr, double makespan);
    void updatePheromoneTrails(const std::vector<std::vector<double>> &pheromone_accumulator);

    void runOptimization();

    template<typename T>
    void setMatrixToZero(std::vector<std::vector<T>> &matrix);
    bool isTabuFull(std::vector<std::vector<int>> &tabu);
    int chooseNextState(std::vector<double> &state_transistion_probs);
    void updateTabu(std::vector<std::vector<int>> &tabu, task* next_task);
    double decodeHeuristicToValue(std::vector<double> &remaining_time_per_job, std::vector<int> &remaining_tasks_per_job,
                                  int decidability_rule, task* next_task);
    void simpleTabuViz(std::vector<std::vector<int>> matrix);
    void printAvailableStateTransitions(std::vector<std::pair<task*, task*>> availableStateTransitions);

};


#endif //PROJECT_3_ACO_H
