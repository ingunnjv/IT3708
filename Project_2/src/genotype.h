#pragma once
#ifndef PROJECT_2_GENOTYPE_H
#define PROJECT_2_GENOTYPE_H

#include <vector>
#include <set>
#include <random>
#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <opencv2/core.hpp>
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
    int num_cols;
    int num_rows;

public:
    std::vector<Genotype*> dominates; // set of solutions that this dominates
    int domination_counter; // number of solution that dominates this solution
    uint16_t rank; // ranges from 0 to maximum of the size of the population
    double crowding_distance; // measure of how far away the genotype is from others in the population
    std::vector<double> objective_values; // values of the two objectives that are optimized
    uint8_t num_objectives; // length of objective_values

    void setRank(int rank);
    void insertToDominationSet(Genotype &i);
    void genotypeToPhenotypeDecoding(int num_rows, int num_cols);
    void visualize(Eigen::MatrixXi &blue_ch, Eigen::MatrixXi &green_ch, Eigen::MatrixXi &red_ch, int num_rows, int num_cols);
    void calculateObjectives(const Eigen::MatrixXi &red, const Eigen::MatrixXi &green, const Eigen::MatrixXi &blue);
    double calcEuclideanRgbDiff(signed short dir_y, signed short dir_x, int this_col, int this_row, int this_segment,
                                const Eigen::MatrixXi &red, const Eigen::MatrixXi &green, const Eigen::MatrixXi &blue);


    Genotype();
    Genotype(uint16_t num_rows, uint16_t num_cols,  std::vector<int> &parents);
    //~Genotype() = default;


    static bool sortByObj1(const Genotype &lhs, const Genotype &rhs);
    static bool sortByObj2(const Genotype &lhs, const Genotype &rhs);
    static bool sortByCrowdedComparison(const Genotype &lhs, const Genotype &rhs);

    bool operator<(const Genotype &rhs) const;
    bool operator>(const Genotype &rhs) const;

};


#endif //PROJECT_2_GENOTYPE_H
