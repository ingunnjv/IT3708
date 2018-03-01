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

void Genotype::genotypeToPhenotypeDecoding(int num_rows, int num_cols)
{
    vector< vector<int> > segments;
    Eigen::MatrixXi segmented_image(num_rows, num_cols);
    segmented_image = Eigen::MatrixXi::Ones(num_rows, num_cols)*(-1);

    int segment_number = 0;
    for (int i = 0; i < this->chromosome.size(); i++) {
        int row = i / num_cols, col = i % num_cols;
        if (segmented_image(row, col) == -1){
            // pixel is not yet assigned to a segment
            segmented_image(row, col) = segment_number;
        }
        else{
            // already assigned to a segment, skip
            continue;
        }
        int next = i;
        while(this->chromosome[next] != none){
            // find next gene
            if (this->chromosome[next] == genValues::right)
                next = next + 1;
            else if (this->chromosome[next] == genValues::left)
                next = next - 1;
            else if (this->chromosome[next] == genValues::down)
                next = next + num_cols;
            else if (this->chromosome[next] == genValues::up)
                next = next - num_cols;
            row = next / num_cols, col = next % num_cols;

            if (segmented_image(row, col) != -1 && segmented_image(row, col) < segment_number){
                // this is part of a previously defined segment. go back
                segment_number = segmented_image(row, col);
                next = i;
                row = next / num_cols, col = next % num_cols;
            }
            segmented_image(row, col) = segment_number;
        }
        segment_number++;
    }
    cout << segmented_image << endl;
}



