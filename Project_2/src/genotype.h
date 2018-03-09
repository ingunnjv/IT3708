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



enum genValues {left = 0, right, up, down, none}; // all possible values of a gene

struct GeneNode {
    uint16_t segment;
    uint8_t value;
};
typedef Eigen::Matrix<struct GeneNode, Eigen::Dynamic, Eigen::Dynamic> GeneMatrix;

class Genotype {
private:
    GeneMatrix chromosome; // storage the entire set of genes

public:
    uint16_t num_cols;
    uint16_t num_rows;
    uint16_t tot_segment_count;

    std::vector<Genotype*> dominates; // set of solutions that this dominates
    uint16_t domination_counter; // number of solution that dominates this solution

    uint16_t rank; // ranges from 0 to maximum of the size of the population
    double crowding_distance; // measure of how far away the genotype is from others in the population

    std::vector<double> objective_values; // values of the two objectives that are optimized
    uint8_t num_objectives; // length of objective_values

    void setRank(int rank);
    void insertToDominationSet(Genotype &i);
    void genotypeToPhenotypeDecoding();
    void visualizeSegments(const Eigen::MatrixXi &blue_ch, const Eigen::MatrixXi &green_ch,
                           const Eigen::MatrixXi &red_ch);
    /// Show two images for each solution in the Pareto-optimal set
    void visualizeEdges(cv::Mat test_image, std::string title);
    bool isEdgePixel(signed short dir_y, signed short dir_x, int this_col, int this_row, int this_segment);
    void calculateObjectives(const Eigen::MatrixXi &red, const Eigen::MatrixXi &green, const Eigen::MatrixXi &blue);
    double calcEuclideanRgbDiff(signed short dir_y, signed short dir_x, int this_col, int this_row, int this_segment,
                                const Eigen::MatrixXi &red, const Eigen::MatrixXi &green, const Eigen::MatrixXi &blue);


    Genotype();
    Genotype(uint16_t num_rows, uint16_t num_cols);
    Genotype(uint16_t num_rows, uint16_t num_cols,  std::vector<int> &parents);
    void setChromosomeValue(uint8_t value, int row, int col);
    void findAndAddParentNodesToList(std::vector<std::tuple<int,int>> &connected_nodes_pos, const std::tuple<int,int> &current_node_pos, Eigen::MatrixXi &discovered);
    void addChildNodeToList(std::vector<std::tuple<int,int>> &connected_nodes_pos, const std::tuple<int,int> &current_node_pos, Eigen::MatrixXi &discovered);
    void setChromosomeSegment(int segment, int row, int col);
    GeneNode * getChromosomeGeneNode(int row, int col);
    uint8_t getChromosomeValue(int row, int col);


    static bool sortByObj1(const Genotype* lhs, const Genotype* rhs);
    static bool sortByObj2(const Genotype* lhs, const Genotype* rhs);
    static bool sortByCrowdedComparison(const Genotype* lhs, const Genotype* rhs);

    bool operator<(const Genotype &rhs) const;
    bool operator>(const Genotype &rhs) const;

};


#endif //PROJECT_2_GENOTYPE_H
