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
    Eigen::MatrixXi red_channel;
    Eigen::MatrixXi green_channel;
    Eigen::MatrixXi blue_channel;
    Eigen::ArrayXi chromosome;

    int minKey(double key[], bool mstSet[]);
public:
    Genotype();
    Genotype(Eigen::MatrixXi red, Eigen::MatrixXi green, Eigen::MatrixXi blue);
    double rgbDistance(pixel_t x, pixel_t y);
    void primMST();
    int printMST(int parent[]);
};


#endif //PROJECT_2_GENOTYPE_H
