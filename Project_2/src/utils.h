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

double rgbDistance(pixel_t x, pixel_t y, const Eigen::MatrixXi &red_channel, const Eigen::MatrixXi &green_channel,
                   const Eigen::MatrixXi &blue_channel);
void setUserArgs(int argc, char **argv, double &mutation_rate, double &crossover_rate, double &tournament_size,
                 double &time_limit, double &generation_limit, double &population_size);


#endif //PROJECT_2_UTILS_H
