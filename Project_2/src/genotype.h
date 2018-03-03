#pragma once
#ifndef PROJECT_2_GENOTYPE_H
#define PROJECT_2_GENOTYPE_H

#include <vector>
#include <set>
#include <random>
#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include "utils.h"



enum genValues {left, right, up, down, none}; // all possible values of a gene
struct GeneNode {
    int segment;
    uint8_t value;
    GeneNode* child;
    std::vector<GeneNode*> parents;
};
typedef Eigen::Matrix<struct GeneNode, Eigen::Dynamic, Eigen::Dynamic> GeneMatrix;

class Genotype {
private:
    GeneMatrix chromosome; // storage the entire set of genes
    uint16_t num_cols;
    uint16_t num_rows;
    std::vector< std::vector<int> > segments; // segments of a solution

public:
    std::vector<Genotype> dominates; // set of solutions that this dominates
    int domination_counter; // number of solution that dominates this solution
    uint16_t rank; // ranges from 0 to maximum of the size of the population
    double crowding_distance; // measure of how far away the genotype is from others in the population
    std::vector<double> objective_values; // values of the two objectives that are optimized
    uint8_t num_objectives; // length of objective_values

    void setRank(int rank);
    void insertToDominationSet(Genotype &i);
    void calculateObjectives(const Eigen::MatrixXi &red, const Eigen::MatrixXi &green, const Eigen::MatrixXi &blue);
    //void genotypeToPhenotypeDecoding(uint16_t num_rows, uint16_t num_cols);
    //void visualize();
    double calcEuclideanRgbDiff(signed short dir_y, signed short dir_x, int this_col, int this_row, int this_segment,
                                const Eigen::MatrixXi &red, const Eigen::MatrixXi &green, const Eigen::MatrixXi &blue);

    Genotype();
    Genotype(uint16_t num_rows, uint16_t num_cols,  std::vector<int> &parents);


    static bool sortByObj1(const Genotype &lhs, const Genotype &rhs);
    static bool sortByObj2(const Genotype &lhs, const Genotype &rhs);
    static bool sortByCrowdedComparison(const Genotype &lhs, const Genotype &rhs);

    bool operator<(const Genotype &rhs) const;
    bool operator>(const Genotype &rhs) const;

};


#endif //PROJECT_2_GENOTYPE_H
