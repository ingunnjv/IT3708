//
// Created by Ingunn on 23.02.2018.
//

#include "nsga2.h"
#include "utils.h"


void Nsga2::mainLoop()
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

void Nsga2::primMST(const Eigen::MatrixXi red, const Eigen::MatrixXi green, const Eigen::MatrixXi blue) {
    // Prim's algorithm
    uint16_t num_rows = uint16_t(red.rows());
    uint16_t num_cols = uint16_t(red.cols());
    uint32_t num_pixels = num_rows * num_cols;

    int parent[num_pixels];   // Array to store constructed MST
    double key[num_pixels];   // Key values used to pick minimum weight edge in cut
    bool mstSet[num_pixels];  // To represent set of vertices not yet included in MST

    // Initialize all keys as INFINITE
    for (int i = 0; i < num_pixels; i++)
        key[i] = INT_MAX, mstSet[i] = false;

    // Always include first vertex in MST.
    key[0] = 0;     // Make key 0 so that this vertex is picked as first vertex
    parent[0] = -1; // First node is always root of MST

    // The MST will have num_pixels vertices
    for (int count = 0; count < num_pixels - 1; count++) {
        // Pick the minimum key vertex from the set of vertices not yet included in MST
        uint32_t u = minKey(key, mstSet, num_pixels);

        // Add the picked vertex to the MST Set
        mstSet[u] = true;

        pixel_t x, y;
        x.row = u / num_cols;
        x.col = u % num_cols;
        bool neighbors[4] = {false, false, false, false};
        uint32_t neighbor_pos[4];
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
        for (uint8_t i = 0; i < 4; i++){
            uint32_t v = neighbor_pos[i];
            // calculate rgbDistance(u, v) only for adjacent vertices of u
            // mstSet[v] is false for vertices not yet included in MST
            // update the key only if rgbDistance(u, v) is smaller than key(v)
            if (neighbors[i] && !mstSet[v]){
                y.row = v / num_cols; // integer division
                y.col = v % num_cols;
                if (rgbDistance(x, y, red, green, blue) < key[v]) {
                    parent[v] = u;
                    key[v] = rgbDistance(x, y, red, green, blue);
                }
            }
        }
    }
}

// A utility function to find the vertex with minimum key value, from the set of vertices not yet included in MST
uint32_t Nsga2::minKey(double key[], bool mstSet[], uint32_t num_pixels)
{
    // Initialize min value
    double min = INT_MAX;
    uint32_t min_index;

    for (uint32_t v = 0; v < num_pixels; v++)
        if (! mstSet[v] && key[v] < min)
            min = key[v], min_index = v;

    return min_index;
}