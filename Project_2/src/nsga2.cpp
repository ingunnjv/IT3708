//
// Created by Ingunn on 23.02.2018.
//
#pragma once
#include "utils.h"
#include "genotype.h"
#include <vector>
#include <Eigen/Dense>
#include <set>

#include "nsga2.h"
#include "utils.h"
#include <set>
#include <iostream>

using namespace std;

/////////////////////////////////////////////////////////
Nsga2::Nsga2()
{
    this->mutation_rate = 0;
    this->crossover_rate = 0;
    this->tournament_size = 0;
    this->generation_limit = 0;
    this->population_size = 0;
    this->population.resize(population_size);
}

/////////////////////////////////////////////////////////
Nsga2::Nsga2(double mutation_rate, double crossover_rate, double tournament_size, double time_limit,
      double generation_limit)
{
    this->mutation_rate = mutation_rate;
    this->crossover_rate = crossover_rate;
    this->tournament_size = tournament_size;
    this->time_limit = time_limit;
    this->generation_limit = generation_limit;
    this->population_size = 10;
    this->population.resize(population_size);
}

/////////////////////////////////////////////////////////
vector< vector<Genotype> > Nsga2::fastNonDominatedSort()
{
    // Create storage for the optimal solutions fronts (worst case if there is only 1 front with the entire pop)
    vector< vector<Genotype> > fronts;//(this->population_size, vector<Genotype>(this->population_size));
    vector<Genotype> first_front;

    // Iterate through the entire population to find who dominates who
    for (auto &p: this->population)
    {
        for (auto &q: this->population)
        {
            if (p < q) // check if p dominates q
            {
                p.insertToDominationSet(q);
            }
            else if (p > q) // check if q dominates p
            {
                p.domination_counter++;
            }
        }
        if (p.domination_counter == 0)
        {
            p.setRank(1);
            first_front.push_back(p); // non-dominated solutions belongs to the first front

        }
    }
    // Iterate through the population again to decide in which front to put each of them
    fronts[0] = first_front;
    int i = 0;  // initialize front counter
    while (!fronts[i].empty())
    {
        vector<Genotype> prev_front = fronts[i];
        vector<Genotype> new_front; // to store members of the next front
        for (auto &p: prev_front)
        {
            for (auto &q: p.dominates)
            {
                q.domination_counter--;
                if (q.domination_counter == 0)
                {
                    q.setRank(i+1);
                    new_front.push_back(q);
                }
            }
        }
        i++;
        fronts[i] = new_front;
    }
    return fronts;
}

/////////////////////////////////////////////////////////
void Nsga2::crowdingDistanceAssignment()
{

}

/////////////////////////////////////////////////////////
void Nsga2::crowdedComparison()
{

}

/////////////////////////////////////////////////////////
void Nsga2::runMainLoop()
{
    // combine parent and offspring population Pt + Qt = Rt
    // fastNonDominatedSort() returns all nondominated fronts of Rt
    // Pt+1 = Ø and i = 1
    // until until the parent population Pt+1 is filled:
    // - - calculate crowding-distance in Fi
    // include ith nondominated front in the parent pop
    // check the next front for inclusion, i = i + 1

    // sort the next nondominated front in descending order <n
    // choose the first (N - |Pt+1|) elements of Fi
    // - - use selection, crossover and mutation to create a new population Qt+1
    // increment the generation counter, t = t + 1
}



struct Cmp
{
    bool operator ()(const std::pair<uint32_t, double> &a, const std::pair<uint32_t, double> &b) {
        return a.second <= b.second;
    }
};
/////////////////////////////////////////////////////////
/// Prim's algorithm
void Nsga2::primMST(Eigen::MatrixXi &red, Eigen::MatrixXi &green, Eigen::MatrixXi &blue) {
    uint16_t num_rows = uint16_t(red.rows());
    uint16_t num_cols = uint16_t(red.cols());
    uint32_t num_pixels = num_rows * num_cols;

    int parent[num_pixels];   // Array to store constructed MST
    double key[num_pixels];   // Key values used to pick minimum weight edge in cut
    bool mstSet[num_pixels];  // To represent set of vertices not yet included in MST
    std::set <std::pair <uint32_t, double>, Cmp> vertices_considered;
    pixel_t x, y;

    // Initialize all keys as INFINITE
    for (int i = 0; i < num_pixels; i++)
        key[i] = INT_MAX, mstSet[i] = false;

    // Always include first vertex in MST.
    key[0] = 0;     // Make key 0 so that this vertex is picked as first vertex
    parent[0] = -1; // First node is always root of MST
    vertices_considered.insert(std::make_pair(0, key[0]));

    // The MST will have num_pixels vertices
    for (int count = 0; count < num_pixels - 1; count++) {
        // Pick the minimum key vertex from the set of vertices not yet included in MST
        auto it = vertices_considered.begin();
        while (mstSet[it->first]){
            // this pixel has already been added to the tree
            vertices_considered.erase(it);
            it = vertices_considered.begin();
        }
        uint32_t u = it->first;
        vertices_considered.erase(it);

        // Add the picked vertex to the MST Set
        mstSet[u] = true;

        x.row = u / num_cols;
        x.col = u % num_cols;
        bool neighbors[4] = {false, false, false, false};
        static uint32_t neighbor_pos[4];
        if (u % num_cols != 0) {
            neighbors[0] = true;
            neighbor_pos[0] = u - 1;
        }
        if ((u + 1) % num_cols != 0) {
            neighbors[1] = true;
            neighbor_pos[1] = u + 1;
        }
        if (u >= num_cols) {
            neighbors[2] = true;
            neighbor_pos[2] = u - num_cols;
        }
        if (u + num_cols < num_pixels) {
            neighbors[3] = true;
            neighbor_pos[3] = u + num_cols;
        }

        // Update key value and parent index of the adjacent vertices of the picked vertex.
        uint32_t v;
        for (uint8_t i = 0; i < 4; i++){
            v = neighbor_pos[i];
            // calculate rgbDistance(u, v) only for adjacent vertices of u
            // mstSet[v] is false for vertices not yet included in MST
            // update the key only if rgbDistance(u, v) is smaller than key(v)
            if (neighbors[i] && !mstSet[v]){

                y.row = v / num_cols; // integer division
                y.col = v % num_cols;
                double dist = rgbDistance(x, y, red, green, blue);
                if (dist < key[v]) {
                    parent[v] = u;
                    key[v] = dist;
                }
                vertices_considered.insert(std::make_pair(v, key[v]));
            }
        }
        //if (!vertices_considered.empty()){ std::cout << "set size: " << vertices_considered.size() << std::endl; }
        delete[] neighbors;
        delete[] neighbor_pos;
    }
}

