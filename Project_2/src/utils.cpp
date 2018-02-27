//
// Created by Ingunn on 27.02.2018.
//

#include "utils.h"

double rgbDistance(pixel_t x, pixel_t y, Eigen::MatrixXi &red, Eigen::MatrixXi &green, Eigen::MatrixXi &blue)
{
    double dist = sqrt( pow(red(x.row, x.col) - red(y.row, y.col), 2)
                        + pow(green(x.row, x.col) - green(y.row, y.col), 2)
                        + pow(blue(x.row, x.col) - blue(y.row, y.col), 2) );
    return dist;
}