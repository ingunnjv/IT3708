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
    Eigen::MatrixXd red_channel;
    Eigen::MatrixXd green_channel;
    Eigen::MatrixXd blue_channel;
    Eigen::ArrayXi chromosome;

    int minKey(double key[], bool mstSet[]);
public:
    Genotype();
    Genotype(Eigen::MatrixXd red, Eigen::MatrixXd green, Eigen::MatrixXd blue);
    double rgbDistance(pixel_t x, pixel_t y);
    void primMST();
    int printMST(int parent[]);
};


#endif //PROJECT_2_GENOTYPE_H
