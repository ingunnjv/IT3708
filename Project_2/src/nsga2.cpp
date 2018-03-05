#include "nsga2.h"

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
Nsga2::Nsga2(double mutation_rate, double crossover_rate, uint16_t tournament_size, double time_limit,
             uint16_t generation_limit, uint16_t population_size)
{
    this->mutation_rate = mutation_rate;
    this->crossover_rate = crossover_rate;
    this->tournament_size = tournament_size;
    this->time_limit = time_limit;
    this->generation_limit = generation_limit;
    this->population_size = population_size;
    this->population.resize(population_size);
}

/////////////////////////////////////////////////////////
void Nsga2::initializePopulation(const Eigen::MatrixXi &red, const Eigen::MatrixXi &green, const Eigen::MatrixXi &blue)
{
    int num_rows = int(red.rows());
    int num_cols = int(red.cols());
    int num_pixels = num_rows * num_cols;
    ///In the initialization of
    ///the ith individual in the population, the (i−1) long links are
    ///removed from the MST individual.
    vector<int> parent_graph = primMST(red, green, blue);
    set < pair<uint32_t, double>, pairCmpGe > links;

    pixel_t x, y;
    for (int i = 1; i < num_pixels; i++) {
        x.row = i / num_cols, x.col = i % num_cols;
        y.row = parent_graph[i] / num_cols, y.col = parent_graph[i] % num_cols;
        links.insert(make_pair(i, rgbDistance(y, x, red, green, blue)));
    }
    for (int i = 0; i < population_size; i++){
        population[i] = Genotype(num_rows, num_cols, parent_graph);
        auto it = links.begin();
        parent_graph[it->first] = -1;
        if (!links.empty()) { links.erase(it); }
        population[i].genotypeToPhenotypeDecoding(num_rows, num_cols);
        population[i].calculateObjectives(red, green, blue);
        cout << "Created genotype with " << i << " segments." << endl;
    }
}

/////////////////////////////////////////////////////////
vector<vector<Genotype> > Nsga2::fastNonDominatedSort() {
    vector< vector<Genotype> > fronts;//(this->population_size, vector<Genotype>(this->population_size));
    fronts.reserve(this->population_size);
    vector<Genotype> first_front;
    first_front.reserve(this->population_size);

    // Iterate through the entire population to find who dominates who
     for (auto &p: this->population) {
        for (auto &q: this->population) {
            if (p < q){ // check if p dominates q
                p.insertToDominationSet(q);
            }
            else if (p > q){ // check if q dominates p
                p.domination_counter++;
            }
        }
        if (p.domination_counter == 0){ // if p is dominated by no one, it means it belongs to the first rank
            p.setRank(0);
            first_front.push_back(p); // non-dominated solutions belongs to the first front
        }
    }
    // Iterate through the population again to decide in which front to put each of them
    fronts.push_back(first_front);
    int i = 0; // initialize front counter
    while (!fronts[i].empty()){ // if last iteration gave no new front, the prev front was the last front
        vector<Genotype> new_front; // to store members of the next front
        new_front.reserve(this->population_size);
        for (auto &p: fronts[i]) {
            for (auto *q: p.dominates) {
                q->domination_counter--; // removing the counts given by the prev front
                if (q->domination_counter == 0){ // q must belong to the next front if its counter is now zero
                    q->setRank(i+1);
                    new_front.push_back(*q);
                }
            }
        }
        i++;
        fronts.push_back(new_front);
    }
    return fronts;
}

/////////////////////////////////////////////////////////
tuple<double, double> Nsga2::objectiveValueSort(std::vector<Genotype> &front, uint8_t obj_val_i)
{
    if (obj_val_i == 0){
        sort(front.begin(), front.end(), Genotype::sortByObj1);
        double fmin = front[0].objective_values[0];
        double fmax = front.back().objective_values[0];
        tuple<double, double> extreme_vals (make_tuple(fmin, fmax));
        return extreme_vals;
    }
    else if (obj_val_i == 1) {
        sort(front.begin(), front.end(), Genotype::sortByObj2);
        double fmin = front[0].objective_values[1];
        double fmax = front.back().objective_values[1];
        tuple<double, double> extreme_vals (make_tuple(fmin, fmax));
        return extreme_vals;
    }
}

/////////////////////////////////////////////////////////
void Nsga2::crowdingDistanceSort(std::vector<Genotype> &front)
{
    sort(front.begin(), front.end(), Genotype::sortByCrowdedComparison);
}

/////////////////////////////////////////////////////////
void Nsga2::crowdingDistanceAssignment(vector<Genotype> &front)
{
    vector<Genotype>::size_type front_size = front.size(); // number of solutions in the front
    uint8_t num_objectives = front[0].num_objectives; // number of objectives in problem
    for (auto &genotype: front){ // initialize all crowding distances to zero
        genotype.crowding_distance = 0;
    }
    for (uint8_t obj_val_num = 0; obj_val_num != num_objectives; obj_val_num++){
        // sort on objective value number and return min and max of the objective
        tuple<double, double> extreme_vals = objectiveValueSort(front, obj_val_num);
        double fmin = get<0>(extreme_vals);
        double fmax = get<1>(extreme_vals);
        // set first and last genotypes distance to inf such that boundary points are always selected
        front[0].crowding_distance = DBL_MAX;
        front.back().crowding_distance = DBL_MAX;
        for (vector<Genotype>::size_type i = 1; i != front_size - 1; i++)
        {
            if (fmin-fmax != 0) {
                front[i].crowding_distance +=
                        (front[i + 1].objective_values[obj_val_num] - front[i - 1].objective_values[obj_val_num]) /
                        (fmax - fmin);
            }
            else{ // avoids division by zero
                front[i].crowding_distance = DBL_MAX;
            }
        }
    }
}

/////////////////////////////////////////////////////////
void Nsga2::runMainLoop(const Eigen::MatrixXi &red, const Eigen::MatrixXi &green, const Eigen::MatrixXi &blue)
{
    vector< vector<Genotype> > fronts;
    fronts.reserve(this->population_size);
    vector<Genotype> next_gen;
    next_gen.reserve(this->population_size);
    vector<Genotype> offspring_pop;
    offspring_pop.reserve(this->population_size);

    int generation = 0;
    while (generation < generation_limit)
    {
        //TODO: implement some sort of early stopping by comparing solutions with PRI

        /* Create the next gen pop by the non-domination principle */
        fronts = fastNonDominatedSort();
        next_gen.clear();
        int i = 0;
        while (next_gen.size() + fronts[i].size() < population_size)
        {
            //crowdingDistanceAssignment(fronts[i]);
            next_gen.insert(next_gen.end(), fronts[i].begin(), fronts[i].end());
            i++;
        }
        //crowdingDistanceAssignment(fronts[i]);
        //crowdingDistanceSort(fronts[i]);
        next_gen.insert(next_gen.end(), fronts[i].begin(), fronts[i].begin() + (population_size - uint16_t(next_gen.size())));

        /* Create a new offspring population by crossover and mutation */
        offspring_pop = makeNewPop(next_gen); // TODO: Implement this function (+crossover & mutation)

        /* Combine the parent and offspring population for next iteration */
        population.clear();
        population.insert(population.end(), next_gen.begin(), next_gen.end());
        population.insert(population.end(), offspring_pop.begin(), offspring_pop.end());
        generation++;
    }

    // TODO: show solution
    // TODO: get and print PRI of the best solution
}

/////////////////////////////////////////////////////////
vector<Genotype> Nsga2::makeNewPop(std::vector<Genotype> &parent_pop) {
    return parent_pop;
}


/////////////////////////////////////////////////////////
vector<int> Nsga2::primMST(const Eigen::MatrixXi &red, const Eigen::MatrixXi &green, const Eigen::MatrixXi &blue) {
    auto num_rows = uint16_t(red.rows());
    auto num_cols = uint16_t(red.cols());
    uint32_t num_pixels = num_rows * num_cols;

    vector<int> parent(num_pixels);   // Array to store constructed MST
    double key[num_pixels];   // Key values used to pick minimum weight edge in cut
    bool mstSet[num_pixels];  // To represent set of vertices not yet included in MST
    set <pair <uint32_t, double>, pairCmpLe> vertices_considered;
    pixel_t x, y;

    // Initialize all keys as INFINITE
    for (int i = 0; i < num_pixels; i++)
        key[i] = INT_MAX, mstSet[i] = false;

    // Always include first vertex in MST.
    key[0] = 0;     // Make key 0 so that this vertex is picked as first vertex
    parent[0] = -1; // First node is always root of MST
    vertices_considered.insert(make_pair(0, key[0]));

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
                vertices_considered.insert(make_pair(v, key[v]));
            }
        }
        //if (!vertices_considered.empty()){ std::cout << "set size: " << vertices_considered.size() << std::endl; }
        //delete[] neighbors;
        //delete[] neighbor_pos;
    }
    return parent;
}
