//
// Created by Ingunn on 27.02.2018.
//
#pragma once
#ifndef PROJECT_2_UTILS_H
#define PROJECT_2_UTILS_H


#include <Eigen/Dense>

struct pixel_t {
    uint16_t row;
    uint16_t col;
};

double rgbDistance(pixel_t x, pixel_t y, const Eigen::MatrixXi &red, const Eigen::MatrixXi &green,
                   const Eigen::MatrixXi &blue);
void setUserArgs(const int argc, char **argv, double &mutation_rate, double &crossover_rate,
                 uint16_t &tournament_size,
                 double &time_limit, uint16_t &generation_limit, uint16_t &population_size);


#endif //PROJECT_2_UTILS_H
