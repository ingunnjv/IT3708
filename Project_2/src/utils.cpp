//
// Created by LasseBot on 27-Feb-18.
//
#pragma once
#include <Eigen/Dense>
#include "utils.h"

/////////////////////////////////////////////////////////
double rgbDistance(pixel_t x, pixel_t y, const Eigen::MatrixXi &red_channel, const Eigen::MatrixXi &green_channel,
                   const Eigen::MatrixXi &blue_channel)
{
    double dist = sqrt( pow(red_channel(x.row, x.col) - red_channel(y.row, y.col), 2)
                        + pow(green_channel(x.row, x.col) - green_channel(y.row, y.col), 2)
                        + pow(blue_channel(x.row, x.col) - blue_channel(y.row, y.col), 2) );
    return dist;
}

/////////////////////////////////////////////////////////
void setUserArgs(int argc, char **argv, double &mutation_rate, double &crossover_rate, double &tournament_size,
                 double &time_limit, double &generation_limit, double &population_size)
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
            tournament_size = arg_val;
            continue;
        }
        else if(!strcmp("time_limit", arg_name))
        {
            time_limit = arg_val;
            continue;
        }
        else if(!strcmp("generation_limit", arg_name))
        {
            generation_limit = arg_val;
            continue;
        }
        else if(!strcmp("population_size", arg_name))
        {
            population_size = arg_val;
            continue;
        }
    }
}