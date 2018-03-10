#include "nsga2.h"
#include <stdlib.h>
#include <random>
#include <chrono>
#include <ctime>
#include <algorithm>
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
void Nsga2::initializePopulation(const Eigen::MatrixXi &red, const Eigen::MatrixXi &green, const Eigen::MatrixXi &blue,
                                 vector<Genotype> &initial_pop)
{
    uint16_t num_rows = uint16_t(red.rows());
    uint16_t num_cols = uint16_t(red.cols());
    int num_pixels = num_rows * num_cols;
    vector<int> parent_graph = primMST(red, green, blue);
    set < pair<uint32_t, double>, pairCmpGe > links;

    pixel_t x(0,0), y(0,0);
    for (int i = 0; i < num_pixels; i++) {
        if (parent_graph[i] == -1){
            continue;
        }
        x.row = i / num_cols, x.col = i % num_cols;
        y.row = parent_graph[i] / num_cols, y.col = parent_graph[i] % num_cols;
        links.insert(make_pair(i, rgbDistance(y, x, red, green, blue)));
    }
    for (int i = 0; i < population_size; i++){
        //this->population[i] = &initial_pop[i];
        // Remove largest edges
        //TODO: how many edges??
        if (i > 0){
            for(uint8_t removed_links = 0; removed_links < 10; removed_links++){
                if (!links.empty()){
                    auto it = links.begin();
                    parent_graph[it->first] = -1;
                    links.erase(it);
                }
            }
        }

        initial_pop[i] = Genotype(num_rows, num_cols, parent_graph);
        initialMutation(initial_pop[i]);
        initial_pop[i].decodeAndEvaluate(red, green, blue);
        printf("+ Created genotype with %d segments(s).\n", i+1);
    }
}

/////////////////////////////////////////////////////////
void Nsga2::fastNonDominatedSort(vector< vector<Genotype*> > &fronts) {
    vector<Genotype*> first_front;
    first_front.reserve(this->population_size);

    printf("\t+ Do fast non-domination sort of solutions...\n");
    // Iterate through the entire population to find who dominates who
    int p_i = 0, q_i = 0;
    for (auto &p: this->population) {
        p.dominates.clear();
        p.domination_counter = 0;
        for (auto &q: this->population) {
            if (p_i != q_i){
                if (p < q){ // check if p dominates q
                    p.insertToDominationSet(q);
                }
                else if (p > q){ // check if q dominates p
                    p.domination_counter++;
                }
            }
            q_i++;
        }
        if (p.domination_counter == 0){ // if p is dominated by no one, it means it belongs to the first rank
            p.setRank(0);
            first_front.push_back(&p); // non-dominated solutions belongs to the first front
        }
        p_i++;
        q_i = 0;
    }
    // Iterate through the population again to decide in which front to put each of them
    fronts.push_back(first_front);
    bool new_front_is_empty = true;
    int i = 0; // initialize front counter
    while (true){ // if last iteration gave no new front, the prev front was the last front
        std::vector<Genotype*> new_front; // to store members of the next front
        new_front.reserve(this->population_size);
        for (auto &p: fronts[i]) {
            for (auto &q: p->dominates) {
                q->domination_counter--; // removing the counts given by the prev front
                if (q->domination_counter == 0){ // q must belong to the next front if its counter is now zero
                    q->setRank(i+1);
                    new_front.push_back(q);
                    new_front_is_empty = false;
                }
            }
            p->dominates.clear();
            p->domination_counter = 0;
        }
        i++;
        if(!new_front_is_empty) {
            fronts.push_back(new_front);
            new_front_is_empty = true;
        }
        else{break;}
    }
}

/////////////////////////////////////////////////////////
tuple<double, double> Nsga2::objectiveValueSort(std::vector<Genotype*> &front, uint8_t obj_val_i)
{
    if (obj_val_i == 0){
        sort(front.begin(), front.end(), Genotype::sortByObj1);
        double fmin = front[0]->objective_values[0];
        double fmax = front.back()->objective_values[0];
        tuple<double, double> extreme_vals (make_tuple(fmin, fmax));
        return extreme_vals;
    }
    else if (obj_val_i == 1) {
        sort(front.begin(), front.end(), Genotype::sortByObj2);
        double fmin = front[0]->objective_values[1];
        double fmax = front.back()->objective_values[1];
        tuple<double, double> extreme_vals (make_tuple(fmin, fmax));
        return extreme_vals;
    }
}

/////////////////////////////////////////////////////////
void Nsga2::crowdingDistanceSort(std::vector<Genotype*> &front)
{
    sort(front.begin(), front.end(), Genotype::sortByCrowdedComparison);
}

/////////////////////////////////////////////////////////
void Nsga2::crowdingDistanceAssignment(vector<Genotype*> &front)
{
    vector<Genotype>::size_type front_size = front.size(); // number of solutions in the front
    uint8_t num_objectives = front[0]->num_objectives; // number of objectives in problem
    for (auto &genotype: front){ // initialize all crowding distances to zero
        genotype->crowding_distance = 0;
    }
    for (uint8_t obj_val_num = 0; obj_val_num != num_objectives; obj_val_num++){
        // sort on objective value number and return min and max of the objective
        tuple<double, double> extreme_vals = objectiveValueSort(front, obj_val_num);
        double fmin = get<0>(extreme_vals);
        double fmax = get<1>(extreme_vals);
        // set first and last genotypes distance to inf such that boundary points are always selected
        front[0]->crowding_distance = DBL_MAX;
        front.back()->crowding_distance = DBL_MAX;
        for (vector<Genotype>::size_type i = 1; i < front_size - 1; i++)
        {
            if (fmin-fmax != 0 && front[i]->crowding_distance != DBL_MAX) {
                front[i]->crowding_distance +=
                        abs((front[i + 1]->objective_values[obj_val_num] - front[i - 1]->objective_values[obj_val_num]) /
                        (fmax - fmin));
            }
        }
    }
}

/////////////////////////////////////////////////////////
void Nsga2::runMainLoop(const Eigen::MatrixXi &red, const Eigen::MatrixXi &green, const Eigen::MatrixXi &blue,
                        vector<Genotype> &initial_pop, cv::Mat &image) {
    /* Create and reserve storage space (big objects) */
    printf("+ Allocating memory for population buffers..\n");
    vector< vector<Genotype*> > fronts;//(this->population_size, vector<Genotype*>(this->population_size));
    vector<Genotype> parents_pop;
    vector<Genotype> offspring_pop;

    for (int i = 0; i < population_size; i++){
        population[i] = initial_pop[i];
        //population[i].visualizeEdges(image, "Initial segmentation");
    }

    initial_pop.erase(initial_pop.begin(), initial_pop.end());

    /* Run evolutionary process */
    int generation = 1;
    while (generation <= generation_limit)
    {
        /* Print process info every x'th generation */
        //TODO: implement some sort of early stopping by comparing solutions with PRI
        if (generation % 1 == 0){
            printf("+ Generation: %d (Best MPRI score: )\n", generation);
        }

        /* Create the next gen pop by the non-domination principle */
        fronts.clear();
        fastNonDominatedSort(fronts);
//        printf("\t+ Visualizing the final pareto optimal fronts..\n");
//        for(auto &genotype: fronts[0]){
//            genotype->visualizeSegments(red, green, blue);
//        }
        int front_index = 0;
        int parent_pop_index = 0;
        while (parent_pop_index + fronts[front_index].size() <= population_size && !fronts[front_index].empty())
        {
            crowdingDistanceAssignment(fronts[front_index]);
            for (const auto &genotype: fronts[front_index]){
                parents_pop.push_back(*genotype);
                parent_pop_index++;
            }
            front_index++;
        }
        /* Fill in the remaining population space with as many genotypes front front_index i there is space for,
         * prioritizing the ones with the best crowding_distance */
        if (parent_pop_index < population_size && !fronts[front_index].empty()){
            crowdingDistanceAssignment(fronts[front_index]);
            crowdingDistanceSort(fronts[front_index]);
        }
        int i = 0;
        while(parent_pop_index < population_size && !fronts[front_index].empty()){
            parents_pop.push_back(*fronts[front_index][i]);
            parent_pop_index++;
            i++;
        }

        /* Create a new offspring population by crossover and mutation */
        makeNewPop(red, green, blue, parents_pop, offspring_pop);

//        /* TEST: Visualize segments before termination */
//        if(generation == generation_limit - 1){
//            for(auto &genotype: fronts[0]){
//                genotype->visualizeSegments(red, green, blue);
//            }
//        }

        /* Combine the parent and offspring population */
        population.clear();
        for (i = 0; i < population_size; i++){
            population.push_back(parents_pop[i]);
            population.push_back(offspring_pop[i]);
        }
//        if (generation == 1){population.resize(population_size*2);}
//        int j = 0;
//        for (i = 0; i < population_size; i++){
//            population[i] = parents_pop[j];
//            j++;
//        }
//        j = 0;
//        for (i = population_size; i < 2*population_size; i++){
//            population[i] = offspring_pop[j];
//            j++;
//        }

        parents_pop.clear();
        offspring_pop.clear();
        generation++;


    }

    /* TEST: Visualize segments before termination */
    printf("\t+ Visualizing the final pareto optimal fronts..\n");
    fronts.clear();
    fastNonDominatedSort(fronts);
    int front_num = 0;
    for(auto &front: fronts){
        printf("\t+ Front: %d\n", front_num);
        for(auto &genotype: front){
            genotype->visualizeEdges(image, "...");
        }
        front_num++;
    }



//    fastNonDominatedSort(fronts);
//    int i_f = 0;
//    for (auto &front: fronts){
//        int i_im = 0;
//        for (auto &solution: front){
//            i_im++;
//            string title = "Front " + to_string(i_f)+ " image " + to_string(i_im);
//            solution->visualizeEdges(image, title);
//        }
//        i_f++;
//    }

    this->population.erase(population.begin(), population.end());
    fronts.erase(fronts.begin(), fronts.end());
    parents_pop.erase(parents_pop.begin(), parents_pop.end());
    offspring_pop.erase(offspring_pop.begin(), offspring_pop.end());
    initial_pop.erase(initial_pop.begin(), initial_pop.end());

    // TODO: show solution
    // TODO: get and print PRI of the best solution
}

/////////////////////////////////////////////////////////
void Nsga2::makeNewPop(const Eigen::MatrixXi &red, const Eigen::MatrixXi &green, const Eigen::MatrixXi &blue,
                       vector<Genotype> &parent_pop, vector<Genotype> &offspring_pop) {
    unsigned seed = (unsigned) chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    uniform_real_distribution<double> rand_distribution(0.0, 1.0);

//    vector<Genotype> offspring(2);
    vector<Genotype *> selected_parents(2);
    bool new_offspring_generated;
    printf("\t+ Make new offspring population...\n");

    for (int i = 0; i < population_size; i = i + 2) {
        new_offspring_generated = false;
        tournamentSelection(selected_parents, parent_pop);
        for (auto p: selected_parents) {
            offspring_pop.push_back(*p);
        }
        double crossover_event = rand_distribution(generator);
        if (crossover_event < this->crossover_rate) {
            uniformCrossover(offspring_pop[i], offspring_pop[i + 1]);
            new_offspring_generated = true;
        }
        for (int i_offspring = i; i_offspring < i + 2; i_offspring++){
            double mutation_event = rand_distribution(generator);
            if (mutation_event < this->mutation_rate || !new_offspring_generated) {
                mutation(offspring_pop[i_offspring], red, green, blue);
                new_offspring_generated = true;
            }

            // update segments & objectives
            offspring_pop[i_offspring].decodeAndEvaluate(red, green, blue);
        }
    }
}

void Nsga2::tournamentSelection(vector<Genotype *> &selected_parents, vector<Genotype> &parent_pop)
{
    unsigned seed = (unsigned)chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator (seed);
    // uniform distribution of integers in the range [0, this->population_size - 1]
    uniform_int_distribution<int> distribution(0, this->population_size - 1);

    vector<int> tournament_indices;
    vector<Genotype*> tournament_participants;
    int num_tournaments = 2;

    for (int i = 0; i < num_tournaments; i++){
        tournament_indices.clear();
        tournament_participants.clear();
        // Choose the participants
        while (tournament_indices.size() < this->tournament_size){
            int random = distribution(generator);
            if(find(tournament_indices.begin(), tournament_indices.end(), random) != tournament_indices.end()){
                // tournament_indices contains random already
                continue;
            }
            else{
                tournament_indices.push_back(random);
                // Add to participants list
                tournament_participants.push_back(&parent_pop[random]);
            }
        }

        // Sort and choose the best
        crowdingDistanceSort(tournament_participants);
        selected_parents[i] = tournament_participants[0];

    }
//    while (tournament_indices.size() < this->tournament_size){
//        int random = distribution(generator);
//        if(find(tournament_indices.begin(), tournament_indices.end(), random) != tournament_indices.end()){
//            // tournament_indices contains random already
//            continue;
//        }
//        else{
//            tournament_indices.push_back(random);
//        }
//    }

//    // clear parents from previous tournament
//    selected_parents.clear();
//    selected_parents.resize(this->tournament_size);
//    for (int i = 0; i < this->tournament_size; i++){
//        selected_parents[i] = &parent_pop[tournament_indices[i]];
//    }
//    crowdingDistanceSort(selected_parents);
//
//    // keep only the 2 best individuals
//    selected_parents.resize(2);
}

void Nsga2::uniformCrossover(Genotype &offspring1, Genotype &offspring2)
{
    unsigned seed = (unsigned)chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    uniform_real_distribution<double> distribution(0.0, 1.0);

    double random_variable;
    for (int i = 0; i < offspring1.num_cols * offspring1.num_rows; i++){
        int row = i / offspring1.num_cols, col = i % offspring1.num_cols;
        random_variable = distribution(generator);
        if(offspring1.getChromosomeValue(row, col) != offspring2.getChromosomeValue(row, col)){
            if (random_variable < 0.5){
                uint8_t value1 = offspring1.getChromosomeValue(row, col);
                uint8_t value2 = offspring2.getChromosomeValue(row, col);
                offspring1.setChromosomeValue(value2, row, col);
                offspring2.setChromosomeValue(value1, row, col);
            }
        }
    }
}

void Nsga2::initialMutation(Genotype &individual) {
    unsigned seed = (unsigned) chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);


    uniform_int_distribution<uint8_t> gene_value_distribution(0, 3);
    uniform_real_distribution<double> rand_distribution(0.0, 1.0);

    for (int row = 0; row < individual.num_rows; row++) {
        for (int col = 0; col < individual.num_cols; col++) {
            double mutate = rand_distribution(generator);
            if (mutate < 0.0001) {
                // Choose random gene value
                uint8_t gene_value = individual.getChromosomeValue(row, col);
                uint8_t random_gene_value = gene_value_distribution(generator);
                while (random_gene_value == gene_value) {
                    random_gene_value = gene_value_distribution(generator);
                }
                individual.setChromosomeValue(random_gene_value, row, col);
            }
        }
    }
}
void Nsga2::mutation(Genotype &individual, const Eigen::MatrixXi &red, const Eigen::MatrixXi &green, const Eigen::MatrixXi &blue)
{
    unsigned seed = (unsigned) chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    uniform_int_distribution<int> pixel_distribution(0, individual.num_rows * individual.num_cols - 1);
    uniform_int_distribution<uint8_t> gene_value_distribution(0, 3);
    uniform_real_distribution<double> rand_distribution(0.0, 1.0);

    vector<int> mutation_indices;
    while (mutation_indices.size() < (int)0.01*individual.num_rows * individual.num_cols){
        int random = pixel_distribution(generator);
        mutation_indices.push_back(random);
    }
    for (auto &gene: mutation_indices){
        // Choose random gene
        int row = gene / individual.num_cols, col = gene % individual.num_cols;

        double smart_mutation = rand_distribution(generator);
        if(smart_mutation < 0.75){
            // Choose value based on the neighbouring pixels
            tuple<uint16_t, uint16_t> most_similar_neighbour_pos = getMostSimilarNeighbourPixel(uint16_t(row), uint16_t(col), red, green, blue);
            uint16_t most_similar_row = get<0>(most_similar_neighbour_pos);
            uint16_t most_simiar_col = get<1>(most_similar_neighbour_pos);
            uint8_t most_similar_value = individual.getChromosomeValue(most_similar_row, most_simiar_col);
            individual.setChromosomeValue(most_similar_value, row, col);
        }
        else{
            // Set a random value
            uint8_t genValue = (uint8_t(3.0*rand_distribution(generator)));
            individual.setChromosomeValue(genValue, row, col);
        }
    }
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
    pixel_t x(0,0), y(0,0);

    // Initialize all keys as INFINITE
    for (int i = 0; i < num_pixels; i++)
        key[i] = DBL_MAX, mstSet[i] = false;

    // Always include initial random vertex in MST.
    unsigned seed = (unsigned)chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator (seed);
    uniform_int_distribution<int> vertex_distribution(0, num_rows * num_cols - 1);
    int random_vertex = vertex_distribution(generator);
    key[random_vertex] = 0;     // Make key 0 so that this vertex is picked as first vertex
    parent[random_vertex] = -1; // First node is always root of MST
    vertices_considered.insert(make_pair(random_vertex, key[random_vertex]));

    // The MST will have num_pixels vertices
    for (int count = 0; count < num_pixels - 1; count++) {
        // Pick the minimum key vertex from the set of vertices not yet included in MST
        auto it = vertices_considered.begin();
        while (mstSet[it->first]){  // this pixel has already been added to the tree
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
