#pragma once
#ifndef PROJECT_2_UTILS_H
#define PROJECT_2_UTILS_H

#include <Eigen/Dense>
#include <vector>
#include "genotype.h"

struct pixel_t {
    int row;
    int col;
};

double rgbDistance(pixel_t x, pixel_t y, const Eigen::MatrixXi &red, const Eigen::MatrixXi &green,
                   const Eigen::MatrixXi &blue);
void setUserArgs(int argc, char **argv, double &mutation_rate, double &crossover_rate,
                 uint16_t &tournament_size, double &time_limit, uint16_t &generation_limit,
                 uint16_t &population_size, int &problem_num);
void printMST(std::vector<int> parent, int num_pixels,
              const Eigen::MatrixXi &red, const Eigen::MatrixXi &green, const Eigen::MatrixXi &blue);
bool sortByObj1(const Genotype &lhs, const Genotype &rhs) ;
bool sortByObj2(const Genotype &lhs, const Genotype &rhs);
bool sortByCrowdedComparison(const Genotype &lhs, const Genotype &rhs);

#endif //PROJECT_2_UTILS_H
