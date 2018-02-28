//
// Created by Ingunn on 23.02.2018.
//


#include <random>
#include <Eigen/Dense>
#include <vector>

#include "genotype.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::MatrixXi;

/////////////////////////////////////////////////////////
void Genotype::setRank(int rank)
{
    this->rank = uint16_t(rank);
}

/////////////////////////////////////////////////////////
void Genotype::insertToDominationSet(Genotype &i)
{
    this->dominates.push_back(i);
}

/////////////////////////////////////////////////////////
Genotype::Genotype()
{
    int num_pixels = 10;
    this->chromosome.resize(num_pixels);
    this->num_objectives = 2;
    this->objective_values.resize(this->num_objectives);
    this->domination_counter = 0;
    this->rank = 0;
    this->crowding_distance = 0;

}

/////////////////////////////////////////////////////////
Genotype::Genotype(const MatrixXi &red, const MatrixXi &green, const MatrixXi &blue)
{
    int num_rows = uint16_t(red.rows());
    int num_cols = uint16_t(red.cols());
    int num_pixels = num_rows * num_cols;
    this->chromosome.resize(num_pixels);
    this->num_objectives = 2;
    this->objective_values.resize(this->num_objectives);
    this->domination_counter = 0;
    this->rank = 0;
    this->crowding_distance = 0;
}

/////////////////////////////////////////////////////////
bool Genotype::operator<(const Genotype &left, const Genotype &right) const
{
    for (vector<double>::size_type i = 0; i != left.objective_values.size(); i++)
    {
        if (left.objective_values[i] > right.objective_values[i])
        {
            return false;
        };
    }
    return true;
}

/////////////////////////////////////////////////////////
bool Genotype::operator>(const Genotype &left, const Genotype &right) const
{
    for (vector<double>::size_type i = 0; i != left.objective_values.size(); i++)
    {
        if (left.objective_values[i] < right.objective_values[i])
        {
            return false;
        };
    }
    return true;
}



