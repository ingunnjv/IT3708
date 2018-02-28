//
// Created by LasseBot on 27-Feb-18.
//
#pragma once
#include <Eigen/Dense>
#include "utils.h"

/////////////////////////////////////////////////////////
double rgbDistance(pixel_t x, pixel_t y, const Eigen::MatrixXi &red, const Eigen::MatrixXi &green,
                   const Eigen::MatrixXi &blue)
{
    double dist = sqrt( pow(red(x.row, x.col) - red(y.row, y.col), 2)
                        + pow(green(x.row, x.col) - green(y.row, y.col), 2)
                        + pow(blue(x.row, x.col) - blue(y.row, y.col), 2) );
    return dist;
}

/////////////////////////////////////////////////////////
void setUserArgs(int argc, char **argv, double &mutation_rate, double &crossover_rate,
                 uint16_t &tournament_size, double &time_limit, uint16_t &generation_limit, uint16_t &population_size,
                 int &problem_num)
{
    mutation_rate = 0;
    crossover_rate = 0;
    tournament_size = 0;
    time_limit = 0;
    generation_limit = 0;
    population_size = 0;
    for (uint8_t i = 1; i < argc; i+=2)
    {
        char* arg_name = argv[i];
        double arg_val = atof(argv[i+1]);
        if (!strcmp("mutation_rate", arg_name))
        {
            mutation_rate = arg_val;
            continue;
        }
        else if(!strcmp("crossover_rate", arg_name))
        {
            crossover_rate = arg_val;
            continue;
        }
        else if(!strcmp("tournament_size", arg_name))
        {
            tournament_size = uint16_t(arg_val);
            continue;
        }
        else if(!strcmp("time_limit", arg_name))
        {
            time_limit = arg_val;
            continue;
        }
        else if(!strcmp("generation_limit", arg_name))
        {
            generation_limit = uint16_t(arg_val);
            continue;
        }
        else if(!strcmp("population_size", arg_name))
        {
            population_size = uint16_t(arg_val);
            continue;
        }
        else if(!strcmp("problem_num", arg_name))
        {
            problem_num = int(arg_val);
            continue;
        }
    }
}

/////////////////////////////////////////////////////////
void printMST(std::vector<int> parent, int num_pixels,
              const Eigen::MatrixXi &red, const Eigen::MatrixXi &green, const Eigen::MatrixXi &blue)
{
    pixel_t x, y;
    auto cols = uint16_t(red.cols());
    printf("Edge   Weight\n");
    for (int i = 1; i < num_pixels; i++) {
        x.row = i / cols;
        x.col = i % cols;
        y.row = parent[i] / cols;
        y.col = parent[i] % cols;
        printf("%d - %d    %f \n", parent[i], i, rgbDistance(y, x, red, green, blue));
    }
}

/////////////////////////////////////////////////////////
bool sortByObj1(const Genotype &lhs, const Genotype &rhs) { return lhs.objective_values[0] > rhs.objective_values[0]; }

/////////////////////////////////////////////////////////
bool sortByObj2(const Genotype &lhs, const Genotype &rhs) { return lhs.objective_values[1] > rhs.objective_values[1]; }

/////////////////////////////////////////////////////////
bool sortByCrowdedComparison(const Genotype &lhs, const Genotype &rhs) {
    if (lhs.rank != rhs.rank) {
        return lhs.rank > rhs.rank;
    }
    else if (lhs.rank == rhs.rank)
    {
        return lhs.crowding_distance > rhs.crowding_distance;
    }
}
