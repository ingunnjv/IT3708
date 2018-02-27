//
// Created by Ingunn on 23.02.2018.
//

#ifndef PROJECT_2_GENOTYPE_H
#define PROJECT_2_GENOTYPE_H

#pragma once
#include <Eigen/Dense>



class Genotype {
private:
    enum nodeConnections {left, right, up, down, none};
    uint16_t num_pixels;
    uint16_t num_rows;
    uint16_t num_cols;
   Eigen::ArrayXi chromosome;

public:
    Genotype();
    Genotype(Eigen::MatrixXi red, Eigen::MatrixXi green, Eigen::MatrixXi blue);

};


#endif //PROJECT_2_GENOTYPE_H
