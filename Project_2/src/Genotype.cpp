//
// Created by Ingunn on 23.02.2018.
//

#include "Genotype.h"
using namespace std;
using Eigen::MatrixXd;

Genotype::Genotype()
{
    num_pixels = 10;
    adjacencyMatrix = MatrixXd::Zero(num_pixels, num_pixels);
    chromosome.resize(num_pixels);
}

Genotype::Genotype(MatrixXd red, MatrixXd green, MatrixXd blue)
{
    num_rows = red.rows();
    num_cols = red.cols();
    num_pixels = num_rows * num_cols;
    adjacencyMatrix = MatrixXd::Zero(num_pixels, num_pixels);
    chromosome.resize(num_pixels);
    redChannel = red;
    greenChannel = green;
    blueChannel = blue;
}

double Genotype::RGB_distance(pixel_t x, pixel_t y)
{
    double dist = sqrt( pow(redChannel(x.row, x.col) - redChannel(y.row, y.col), 2)
                        + pow(greenChannel(x.row, x.col) - greenChannel(y.row, y.col), 2)
                        + pow(blueChannel(x.row, x.col) - blueChannel(y.row, y.col), 2) );
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
            adjacencyMatrix(i, i - 1) = RGB_distance(x, y);
        }
        if ((i + 1) % num_cols != 0){
            y.row = (i + 1) / num_rows;
            y.col = (i + 1) % num_cols;
            adjacencyMatrix(i, i + 1) = RGB_distance(x, y);
        }
        if (i >= num_cols){
            y.row = (i - num_cols) / num_rows;
            y.col = (i - num_cols) % num_cols;
            adjacencyMatrix(i, i - num_cols) = RGB_distance(x, y);
        }
        if (i + num_cols < num_pixels){
            y.row = (i + num_cols) / num_rows;
            y.col = (i + num_cols) % num_cols;
            adjacencyMatrix(i, i + num_cols) = RGB_distance(x, y);
        }
    }
}