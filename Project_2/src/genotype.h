#ifndef PROJECT_2_GENOTYPE_H
#define PROJECT_2_GENOTYPE_H

#include <vector>
#include <set>
#include <random>
#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <opencv2/core.hpp>



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
    void genotypeToPhenotypeDecoding(int num_rows, int num_cols);
    void visualize(Eigen::MatrixXi &blue_ch, Eigen::MatrixXi &green_ch, Eigen::MatrixXi &red_ch, int num_rows, int num_cols);

    Genotype();
    Genotype(int num_rows, int num_cols,  std::vector<int> &parents);

    bool operator<(const Genotype &rhs) const;
    bool operator>(const Genotype &rhs) const;

};


#endif //PROJECT_2_GENOTYPE_H
