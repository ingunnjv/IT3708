//
// Created by Ingunn on 23.02.2018.
//

#ifndef PROJECT_2_GENOTYPE_H
#define PROJECT_2_GENOTYPE_H

#pragma once
#include <random>
#include <Eigen/Dense>



class Genotype {
private:
    // Node connections types for genes
    enum nodeConnections {left, right, up, down, none};
    uint16_t num_pixels;
    uint16_t num_rows;
    uint16_t num_cols;
   Eigen::ArrayXi chromosome;

    uint32_t minKey(double key[], bool mstSet[]);
public:
    Genotype();
    Genotype(Eigen::MatrixXi red, Eigen::MatrixXi green, Eigen::MatrixXi blue);

};


#endif //PROJECT_2_GENOTYPE_H
