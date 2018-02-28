#pragma once
#include "genotype.h"

using namespace std;

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
Genotype::Genotype(int num_pixels, int num_cols,  vector<int> &parents)
{
    this->chromosome.resize(num_pixels);
    this->num_objectives = 2;
    this->objective_values.resize(this->num_objectives);
    this->domination_counter = 0;
    this->rank = 0;
    this->crowding_distance = 0;

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
bool Genotype::operator<(const Genotype &right) const
{
    for (vector<double>::size_type i = 0; i != right.objective_values.size(); i++)
    {
        if (this->objective_values[i]  > right.objective_values[i])
        {
            return false;
        };
    }
    return true;
}

/////////////////////////////////////////////////////////
bool Genotype::operator>(const Genotype &right) const
{
    for (vector<double>::size_type i = 0; i != right.objective_values.size(); i++)
    {
        if (this->objective_values[i] < right.objective_values[i])
        {
            return false;
        };
    }
    return true;
}



