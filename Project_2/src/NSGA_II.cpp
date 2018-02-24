//
// Created by Ingunn on 23.02.2018.
//

#include "NSGA_II.h"


void NSGA_II::mainLoop()
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