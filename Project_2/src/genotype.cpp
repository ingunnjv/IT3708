//
// Created by Ingunn on 23.02.2018.
//

#include "genotype.h"
using namespace std;
using Eigen::MatrixXd;

Genotype::Genotype()
{
    num_pixels = 10;
    adjacency_matrix = MatrixXd::Zero(num_pixels, num_pixels);
    chromosome.resize(num_pixels);
}

Genotype::Genotype(MatrixXd red, MatrixXd green, MatrixXd blue)
{
    num_rows = red.rows();
    num_cols = red.cols();
    num_pixels = num_rows * num_cols;
    adjacency_matrix = MatrixXd::Zero(num_pixels, num_pixels);
    chromosome.resize(num_pixels);
    red_channel = red;
    green_channel = green;
    blue_channel = blue;
}

double Genotype::rgbDistance(pixel_t x, pixel_t y)
{
    double dist = sqrt( pow(red_channel(x.row, x.col) - red_channel(y.row, y.col), 2)
                        + pow(green_channel(x.row, x.col) - green_channel(y.row, y.col), 2)
                        + pow(blue_channel(x.row, x.col) - blue_channel(y.row, y.col), 2) );
    return dist;
}

void Genotype::constructAdjacencyMatrix()
{
    pixel_t x, y;

    for (int i = 0; i < num_pixels; i++){
        x.row = i / num_rows; // integer division
        x.col = i % num_cols;
        if (i % num_cols != 0){
            y.row = (i - 1) / num_rows;
            y.col = (i - 1) % num_cols;
            adjacency_matrix(i, i - 1) = rgbDistance(x, y);
        }
        if ((i + 1) % num_cols != 0){
            y.row = (i + 1) / num_rows;
            y.col = (i + 1) % num_cols;
            adjacency_matrix(i, i + 1) = rgbDistance(x, y);
        }
        if (i >= num_cols){
            y.row = (i - num_cols) / num_rows;
            y.col = (i - num_cols) % num_cols;
            adjacency_matrix(i, i - num_cols) = rgbDistance(x, y);
        }
        if (i + num_cols < num_pixels){
            y.row = (i + num_cols) / num_rows;
            y.col = (i + num_cols) % num_cols;
            adjacency_matrix(i, i + num_cols) = rgbDistance(x, y);
        }
    }
}