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
Genotype::Genotype(uint16_t num_rows, uint16_t num_cols,  vector<int> &parents)
{
    this->chromosome.resize(num_rows, num_cols);
    this->num_rows = num_rows;
    this->num_cols = num_cols;
    this->num_objectives = 2;
    this->objective_values.resize(this->num_objectives);
    this->domination_counter = 0;
    this->rank = 0;
    this->crowding_distance = 0;

    for (int row = 0; row < num_rows; row++){
        for (int col = 0; col < num_cols; col++){
            int i = row * num_cols + col;
            if (parents[i] == -1){
                this->chromosome(row, col).value = genValues::none;
                this->chromosome(row, col).child = NULL;
            }
            else{
                if (parents[i] - i == 1) {
                    this->chromosome(row, col).value = genValues::right;
                    this->chromosome(row, col).child = &chromosome(row, col + 1);
                }
                else if (parents[i] - i == -1) {
                    this->chromosome(row, col).value = genValues::left;
                    this->chromosome(row, col).child = &chromosome(row, col - 1);
                }
                else if (parents[i] - i == num_cols) {
                    this->chromosome(row, col).value = genValues::down;
                    this->chromosome(row, col).child = &chromosome(row + 1, col);
                }
                else if (parents[i] - i == -num_cols){
                    this->chromosome(row, col).value = genValues::up;
                    this->chromosome(row, col).child = &chromosome(row - 1, col);
                }
                else{
                    cout << "Error in chromosome initialization" << endl;
                }
                this->chromosome(row, col).child->parents.push_back(&chromosome(row, col));
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

//Genotype::~Genotype(){
//    int rows = this->chromosome.rows();
//    int cols = this->chromosome.cols();
//    for (int row = 0; row < rows; row++) {
//        for (int col = 0; col < cols; col++) {
//            free(this->chromosome(row, col));
//        }
//    }
//}

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


/////////////////////////////////////////////////////////
/*
void Genotype::genotypeToPhenotypeDecoding(uint16_t num_rows, uint16_t num_cols)
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
 */
/*
/////////////////////////////////////////////////////////
void Genotype::visualize() {

}*/

/////////////////////////////////////////////////////////
void Genotype::calculateObjectives(const Eigen::MatrixXi &red, const Eigen::MatrixXi &green, const Eigen::MatrixXi &blue) {
    // find number of segments
    vector<int> segment_nums_found;
    for (int row = 0; row < num_rows; row++) {
        for (int col = 0; col < num_cols; col++) {
            int segment_num = chromosome(row,col).segment;
            if (find(segment_nums_found.begin(), segment_nums_found.end(), segment_num) != segment_nums_found.end()){
                segment_nums_found.push_back(segment_num);
            }
        }
    }
    // create grouping of pixels in corresponding segments
    vector<vector<pixel_t>> pixels_segment_affiliation (segment_nums_found.size());
    for (int row = 0; row < num_rows; row++) {
        for (int col = 0; col < num_cols; col++) {
            int segment_num = chromosome(row,col).segment;
            pixel_t pixel = {row, col};
            pixels_segment_affiliation[segment_num].push_back(pixel);
        }
    }
    // calculate the rgb centroid for each segment
    vector<rgb_centroid_t> centroids (segment_nums_found.size());
    int i = 0;
    for (const auto s: pixels_segment_affiliation){
        rgb_centroid_t centroid = {0, 0, 0};
        double intensity_sum_r = 0;
        double intensity_sum_g = 0;
        double intensity_sum_b = 0;
        int num_pixels = 0;
        for (auto p: s){
            intensity_sum_r += red(p.row, p.col);
            intensity_sum_g += green(p.row, p.col);
            intensity_sum_b += blue(p.row, p.col);
            num_pixels++;
        }
        centroid.r = intensity_sum_r / num_pixels;
        centroid.g = intensity_sum_g / num_pixels;
        centroid.b = intensity_sum_b / num_pixels;
        centroids[i] = centroid;
        i++;
    }
    // calculate objective 1 and 2
    objective_values[0] = 0;
    objective_values[1] = 0;
    i = 0;
    for (const auto s: pixels_segment_affiliation) {
        for (const auto p: s){
            // objective 1
            double pixel_centroid_deviation = sqrt( pow(red(p.row, p.col) - centroids[i].r, 2)
                                                    + pow(green(p.row, p.col) - centroids[i].g, 2)
                                                    + pow(blue(p.row, p.col) - centroids[i].b, 2));
            objective_values[0] += pixel_centroid_deviation;
            // objective 2
            int this_row = p.row;
            int this_col = p.col;
            int this_segment = chromosome(this_row, this_col).segment;
            double segment_boundary_diff = 0;
            segment_boundary_diff += calcEuclideanRgbDiff(-1, 0, this_col, this_row, this_segment, red,green,blue);
            segment_boundary_diff += calcEuclideanRgbDiff(1, 0, this_col, this_row, this_segment, red,green,blue);
            segment_boundary_diff += calcEuclideanRgbDiff(0, -1, this_col, this_row, this_segment, red,green,blue);
            segment_boundary_diff += calcEuclideanRgbDiff(0, 1, this_col, this_row, this_segment, red,green,blue);
            objective_values[1] += segment_boundary_diff;
            }
        i++;
    }
}

/////////////////////////////////////////////////////////
bool Genotype::sortByObj1(const Genotype &lhs, const Genotype &rhs) { return lhs.objective_values[0] < rhs.objective_values[0]; }

/////////////////////////////////////////////////////////
bool Genotype::sortByObj2(const Genotype &lhs, const Genotype &rhs) { return lhs.objective_values[1] < rhs.objective_values[1]; }

/////////////////////////////////////////////////////////
bool Genotype::sortByCrowdedComparison(const Genotype &lhs, const Genotype &rhs) {
    if (lhs.rank != rhs.rank) {
        return lhs.rank > rhs.rank;
    }
    else if (lhs.rank == rhs.rank)
    {
        return lhs.crowding_distance > rhs.crowding_distance;
    }
}

/////////////////////////////////////////////////////////
double Genotype::calcEuclideanRgbDiff(signed short dir_y, signed short dir_x, int this_col, int this_row, int this_segment,
                                      const Eigen::MatrixXi &red, const Eigen::MatrixXi &green, const Eigen::MatrixXi &blue) {
    double diff = 0;
    if (this_col + dir_x > 0 && this_col + dir_x < num_cols && this_row + dir_y > 0 && this_row + dir_y < num_rows) {
        if (chromosome(this_row + dir_y, this_col + dir_x).segment != this_segment) {
            diff = sqrt(pow(red(this_row, this_col) - red(this_row + dir_y, this_col + dir_x), 2)
                        + pow(green(this_row, this_col) - green(this_row + dir_y, this_col + dir_x), 2)
                        + pow(blue(this_row, this_col) - blue(this_row + dir_y, this_col + dir_x), 2));

        }
    }
    return diff;
}



