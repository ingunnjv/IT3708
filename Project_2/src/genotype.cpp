//
// Created by Ingunn on 23.02.2018.
//


#include <random>
#include <Eigen/Dense>
#include <vector>
#include <iostream>

#include "genotype.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::MatrixXi;

/////////////////////////////////////////////////////////
void Genotype::setRank(int rank)
{
    this->rank = rank;
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
    this->objectiveValues.resize(2);
    this->domination_counter = 0;
    this->rank = 0;
}

/////////////////////////////////////////////////////////
Genotype::Genotype(int num_pixels, int num_cols,  vector<int> &parents)
{
    this->chromosome.resize(num_pixels);
    this->objectiveValues.resize(2);
    this->domination_counter = 0;
    this->rank = 0;

    for (int i = 0; i < num_pixels; i++) {

        if (parents[i] == -1){
            this->chromosome[i] = genValues::none;
        }
        else{
            if (parents[i] - i == 1) //right
                this->chromosome[i] = genValues::right;
            else if (parents[i] - i == -1)//left
                this->chromosome[i] = genValues::left;
            else if (parents[i] - i == num_cols)
                this->chromosome[i] = genValues::down;
            else if (parents[i] - i == -num_cols)
                this->chromosome[i] = genValues::up;
            else{
                cout << "Error in chromosome initialization" << endl;
            }
        }
    }
}

/////////////////////////////////////////////////////////
bool operator<(const Genotype &left, const Genotype &right)
{
    for (vector<double>::size_type i = 0; i != left.objectiveValues.size(); i++)
    {
        if (left.objectiveValues[i] > right.objectiveValues[i])
        {
            return false;
        };
    }
    return true;
}

/////////////////////////////////////////////////////////
bool operator>(const Genotype &left, const Genotype &right)
{
    for (vector<double>::size_type i = 0; i != left.objectiveValues.size(); i++)
    {
        if (left.objectiveValues[i] < right.objectiveValues[i])
        {
            return false;
        };
    }
    return true;
}



