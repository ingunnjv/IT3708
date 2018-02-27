//
// Created by Ingunn on 27.02.2018.
//

#ifndef PROJECT_2_UTILS_H
#define PROJECT_2_UTILS_H

#pragma once
#include <Eigen/Dense>

struct pixel_t {
    uint16_t row;
    uint16_t col;
};

double rgbDistance(pixel_t x, pixel_t y, Eigen::MatrixXi &red, Eigen::MatrixXi &green, Eigen::MatrixXi &blue);

#endif //PROJECT_2_UTILS_H
