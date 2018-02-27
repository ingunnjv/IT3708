//
// Created by Ingunn on 23.02.2018.
//

#ifndef PROJECT_2_GENOTYPE_H
#define PROJECT_2_GENOTYPE_H

#pragma once
#include <random>
#include <Eigen/Dense>

struct pixel_t {
    uint16_t row;
    uint16_t col;
};

class Genotype {
private:
    // Node connections types for genes
    enum nodeConnections {left, right, up, down, none};
    uint16_t num_pixels;
    uint16_t num_rows;
    uint16_t num_cols;
    Eigen::MatrixXi red_channel;
    Eigen::MatrixXi green_channel;
    Eigen::MatrixXi blue_channel;
    Eigen::ArrayXi chromosome;

    uint32_t minKey(double key[], bool mstSet[]);
public:
    Genotype();
    Genotype(Eigen::MatrixXi red, Eigen::MatrixXi green, Eigen::MatrixXi blue);
    double rgbDistance(pixel_t x, pixel_t y);
    void primMST();
    int printMST(int parent[]);
};


#endif //PROJECT_2_GENOTYPE_H
