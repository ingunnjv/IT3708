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
Genotype::Genotype(const MatrixXi &red, const MatrixXi &green, const MatrixXi &blue)
{
    int num_rows = uint16_t(red.rows());
    int num_cols = uint16_t(red.cols());
    int num_pixels = num_rows * num_cols;
    this->chromosome;
    this->objectiveValues;
    this->domination_counter = 0;
    this->rank = 0;
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




//void Genotype::primMST() {
//    // Prim's algorithm
//    int parent[num_pixels];   // Array to store constructed MST
//    double key[num_pixels];   // Key values used to pick minimum weight edge in cut
//    bool mstSet[num_pixels];  // To represent set of vertices not yet included in MST
//
//    // Initialize all keys as INFINITE
//    for (int i = 0; i < num_pixels; i++)
//        key[i] = INT_MAX, mstSet[i] = false;
//
//    // Always include first one random vertex in MST.
//
//    // Seed with a real random value, if available
//    random_device r;
//    // Choose a random mean between 0 and num_pixels
//    default_random_engine e1(r());
//    uniform_int_distribution<int> uniform_dist(0, num_pixels);
//    int random_pixel = uniform_dist(e1);
//    key[random_pixel] = 0;     // Make key 0 so that this vertex is picked as first vertex
//    parent[0] = -1; // First node is always root of MST
//
//    // The MST will have num_pixels vertices
//    for (int count = 0; count < num_pixels - 1; count++) {
//        // Pick the minimum key vertex from the set of vertices not yet included in MST
//        int u = minKey(key, mstSet);
//
//        // Add the picked vertex to the MST Set
//        mstSet[u] = true;
//
//        pixel_t x, y;
//        x.row = u / num_cols;
//        x.col = u % num_cols;
//        int neighbors[4] = {false, false, false, false};
//        int neighbor_pos[4];
//        if (u % num_cols != 0) {
//            neighbors[left] = true;
//            neighbor_pos[left] = u - 1;
//        }
//        if ((u + 1) % num_cols != 0) {
//            neighbors[right] = true;
//            neighbor_pos[right] = u + 1;
//        }
//        if (u >= num_cols) {
//            neighbors[up] = true;
//            neighbor_pos[up] = u - num_cols;
//        }
//        if (u + num_cols < num_pixels) {
//            neighbors[down] = true;
//            neighbor_pos[down] = u + num_cols;
//        }
//
//        // Update key value and parent index of the adjacent vertices of the picked vertex.
//        for (int i = 0; i < 4; i++){
//            int v = neighbor_pos[i];
//
//
//            // calculate rgbDistance(u, v) only for adjacent vertices of u
//            // mstSet[v] is false for vertices not yet included in MST
//            // update the key only if rgbDistance(u, v) is smaller than key(v)
//            if (neighbors[i] && !mstSet[v]){
//                y.row = v / num_cols; // integer division
//                y.col = v % num_cols;
//                if (rgbDistance(x, y) < key[v]) {
//                    parent[v] = u;
//                    key[v] = rgbDistance(x, y);
//                }
//            }
//        }
//    }
//}
//


// A utility function to find the vertex with minimum key value, from
// the set of vertices not yet included in MST
//int32_t Genotype::minKey(double key[], bool mstSet[])
//
//   // Initialize min value
//   double min = INT_MAX;
//   uint32_t min_index;

//   for (uint32_t v = 0; v < num_pixels; v++)
//       if (! mstSet[v] && key[v] < min)
//           min = key[v], min_index = v;

//   return min_index;
//

