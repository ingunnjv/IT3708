//
// Created by Ingunn on 23.02.2018.
//
#pragma once
#ifndef PROJECT_2_GENOTYPE_H
#define PROJECT_2_GENOTYPE_H

#include <vector>
#include <set>

enum genValues {left, right, up, down, none}; // all possible values of a gene

class Genotype {
private:
    std::vector<uint8_t> chromosome;                    // storage the entire set of genes
    std::vector<double> objectiveValues;          // values of the two objectives that are optimized


public:
    std::vector<Genotype> dominates;                 // set of solutions that p dominates
    int domination_counter;                          // number of solution that dominates the solution p
    int rank;

    void setRank(int rank);
    void insertToDominationSet(Genotype &i);

    Genotype();
    Genotype(int num_pixels, int num_cols, std::vector<int> &parents);
    friend bool operator<(const Genotype &left, const Genotype &right);
    friend bool operator>(const Genotype &left, const Genotype &right);
};


#endif //PROJECT_2_GENOTYPE_H
