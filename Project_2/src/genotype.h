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
    std::vector<uint8_t> chromosome; // storage the entire set of genes

public:
    std::vector<Genotype> dominates; // set of solutions that this dominates
    int domination_counter; // number of solution that dominates this solution
    uint16_t rank; // ranges from 0 to maximum of the size of the population
    double crowding_distance; // measure of how far away the genotype is from others in the population
    std::vector<double> objective_values; // values of the two objectives that are optimized
    uint8_t num_objectives; // length of objective_values

    void setRank(int rank);
    void insertToDominationSet(Genotype &i);

    Genotype();
    Genotype(const Eigen::MatrixXi &red, const Eigen::MatrixXi &green, const Eigen::MatrixXi &blue);
    bool operator<(const Genotype &right) const;
    bool operator>(const Genotype &right) const;
};


#endif //PROJECT_2_GENOTYPE_H
