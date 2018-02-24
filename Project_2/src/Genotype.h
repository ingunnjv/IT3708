//
// Created by Ingunn on 23.02.2018.
//

#ifndef PROJECT_2_GENOTYPE_H
#define PROJECT_2_GENOTYPE_H

#include <Eigen/Dense>
struct pixel_t {
    int row;
    int col;
};

class Genotype {
private:
    enum nodeConnections {left, right, up, down, none};
    int num_pixels;
    int num_rows;
    int num_cols;
    Eigen::MatrixXd adjacencyMatrix;
    Eigen::MatrixXd redChannel;
    Eigen::MatrixXd greenChannel;
    Eigen::MatrixXd blueChannel;
    Eigen::ArrayXi chromosome;
public:
    Genotype();
    Genotype(Eigen::MatrixXd red, Eigen::MatrixXd green, Eigen::MatrixXd blue);
    double RGB_distance(pixel_t x, pixel_t y);
    void constructAdjacencyMatrix();
};


#endif //PROJECT_2_GENOTYPE_H
