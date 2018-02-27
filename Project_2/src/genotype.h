//
// Created by Ingunn on 23.02.2018.
//
#pragma once
#ifndef PROJECT_2_GENOTYPE_H
#define PROJECT_2_GENOTYPE_H

#include <vector>
#include <set>

class Genotype {
private:
    enum genValues {left, right, up, down, none}; // all possible values of a gene
    std::vector<uint8_t> chromosome;                    // storage the entire set of genes
    std::vector<double> objectiveValues;          // values of the two objectives that are optimized


    //uint32_t minKey(double key[], bool mstSet[]);
public:
    std::vector<Genotype> dominates;                 // set of solutions that p dominates
    int domination_counter;                  // number of solution that dominates the solution p
    int rank;

    void setRank(int rank);
    void insertToDominationSet(Genotype &i);

    Genotype();
    Genotype(const Eigen::MatrixXi &red, const Eigen::MatrixXi &green, const Eigen::MatrixXi &blue);
    friend bool operator<(const Genotype &left, const Genotype &right);
    friend bool operator>(const Genotype &left, const Genotype &right);
};


#endif //PROJECT_2_GENOTYPE_H
