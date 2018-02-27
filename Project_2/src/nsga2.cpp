//
// Created by Ingunn on 23.02.2018.
//

#include "nsga2.h"

/////////////////////////////////////////////////////////
Nsga2::Nsga2()
{
    this->mutation_rate = 0;
    this->crossover_rate = 0;
    this->tournament_size = 0;
    this->generation_limit = 0;
}

/////////////////////////////////////////////////////////
Nsga2::Nsga2(double mutation_rate, double crossover_rate, double tournament_size, double time_limit,
      uint16_t generation_limit)
{
    this->mutation_rate = mutation_rate;
    this->crossover_rate = crossover_rate;
    this->tournament_size = tournament_size;
    this->time_limit = time_limit;
    this->generation_limit = generation_limit;
}

/////////////////////////////////////////////////////////
void Nsga2::fastNonDominatedSort()
{

}

/////////////////////////////////////////////////////////
void Nsga2::crowdingDistanceAssignment()
{

}

/////////////////////////////////////////////////////////
void Nsga2::crowdedComparison()
{

}

/////////////////////////////////////////////////////////
void Nsga2::runMainLoop()
{
    // combine parent and offspring population Pt + Qt = Rt
    // fastNonDominatedSort() returns all nondominated fronts of Rt
    // Pt+1 = Ø and i = 1
    // until until the parent population Pt+1 is filled:
    // - - calculate crowding-distance in Fi
    // include ith nondominated front in the parent pop
    // check the next front for inclusion, i = i + 1

    // sort the next nondominated front in descending order <n
    // choose the first (N - |Pt+1|) elements of Fi
    // - - use selection, crossover and mutation to create a new population Qt+1
    // increment the generation counter, t = t + 1
}