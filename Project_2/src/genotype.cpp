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
    this->chromosome.resize(10,10);
    this->num_objectives = 2;
    this->objective_values.resize(this->num_objectives);
    this->domination_counter = 0;
    this->rank = 0;
    this->crowding_distance = 0;
}

/////////////////////////////////////////////////////////
Genotype::Genotype(int num_rows, int num_cols,  vector<int> &parents)
{
    this->chromosome.resize(num_rows, num_cols);
    this->num_objectives = 2;
    this->objective_values.resize(this->num_objectives);
    this->domination_counter = 0;
    this->rank = 0;
    this->crowding_distance = 0;

    for (int row = 0; row < num_rows; row++){
        for (int col = 0; col < num_cols; col++){
            int i = row * num_cols + col;
            if (parents[i] == -1){
                chromosome(row, col)->value = genValues::none;
            }
            else{
                if (parents[i] - i == 1)
                    this->chromosome[row, col].value = genValues::right;
                else if (parents[i] - i == -1)
                    this->chromosome[row, col].value = genValues::left;
                else if (parents[i] - i == num_cols)
                    this->chromosome[row, col].value = genValues::down;
                else if (parents[i] - i == -num_cols)
                    this->chromosome[row, col].value = genValues::up;
                else{
                    cout << "Error in chromosome initialization" << endl;
                }
            }
        }
    }

    //for (int i = 0; i < num_pixels; i++) {
//
    //    if (parents[i] == -1){
    //        this->chromosome[i] = genValues::none;
    //    }
    //    else{
    //        if (parents[i] - i == 1)
    //            this->chromosome[i] = genValues::right;
    //        else if (parents[i] - i == -1)
    //            this->chromosome[i] = genValues::left;
    //        else if (parents[i] - i == num_cols)
    //            this->chromosome[i] = genValues::down;
    //        else if (parents[i] - i == -num_cols)
    //            this->chromosome[i] = genValues::up;
    //        else{
    //            cout << "Error in chromosome initialization" << endl;
    //        }
    //    }
    //}
}

/////////////////////////////////////////////////////////
bool Genotype::operator<(const Genotype &rhs) const
{
    for (vector<double>::size_type i = 0; i != rhs.objective_values.size(); i++)
    {
        if (this->objective_values[i]  > rhs.objective_values[i])
        {
            return false;
        };
    }
    return true;
}

/////////////////////////////////////////////////////////
bool Genotype::operator>(const Genotype &rhs) const
{
    for (vector<double>::size_type i = 0; i != rhs.objective_values.size(); i++)
    {
        if (this->objective_values[i] < rhs.objective_values[i])
        {
            return false;
        };
    }
    return true;
}
/*
void Genotype::genotypeToPhenotypeDecoding(int num_rows, int num_cols)
{
    Eigen::MatrixXi segmented_image(num_rows, num_cols);
    segmented_image = Eigen::MatrixXi::Ones(num_rows, num_cols)*(-1);

    int segment_number = 0;
    int total_number_of_segments = segment_number;
    for(vector<int>::size_type i = this->chromosome.size() - 1; i != (vector<int>::size_type) - 1; i--){
        int row = i / num_cols, col = i % num_cols;
        if (segmented_image(row, col) == -1){ // pixel is not yet assigned to a segment
            segment_number = total_number_of_segments;
            segmented_image(row, col) = segment_number;
            total_number_of_segments++;
        }
        else{ // already assigned to a segment, skip
            continue;
        }
        int next = i;
        while(this->chromosome[next] != none){ // find next gene
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
                total_number_of_segments--;
                segment_number = segmented_image(row, col);
                next = i;
                row = next / num_cols, col = next % num_cols;
            }
            segmented_image(row, col) = segment_number;
        }
        segment_number++;
    }
    //cout << segmented_image << endl;
    cout << "Total number of segments: " << total_number_of_segments << endl;
    // create vector of segments
    this->segments.resize(total_number_of_segments);
    for (int pixel = 0; pixel < this->chromosome.size(); pixel++){
        int row = pixel / num_cols, col = pixel % num_cols;
        this->segments[segmented_image(row, col)].push_back(pixel);
    }
}

void Genotype::visualize()
{

}

*/

