//
// Created by Ingunn on 27.02.2018.
//

#ifndef PROJECT_2_UTILS_H
#define PROJECT_2_UTILS_H

#pragma once
//#include <cmath>
#include <Eigen/Dense>

struct pixel_t {
    uint16_t row;
    uint16_t col;
};

double rgbDistance(pixel_t x, pixel_t y, Eigen::MatrixXi red_channel, Eigen::MatrixXi green_channel, Eigen::MatrixXi blue_channel)
{
    double dist = sqrt( pow(red_channel(x.row, x.col) - red_channel(y.row, y.col), 2)
                        + pow(green_channel(x.row, x.col) - green_channel(y.row, y.col), 2)
                        + pow(blue_channel(x.row, x.col) - blue_channel(y.row, y.col), 2) );
    return dist;
}

#endif //PROJECT_2_UTILS_H
