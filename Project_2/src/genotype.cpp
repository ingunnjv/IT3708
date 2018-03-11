#include "genotype.h"

#include <opencv2/imgproc.hpp>
#include <opencv2/core/eigen.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgcodecs.hpp>
#include <chrono>

using namespace std;

/////////////////////////////////////////////////////////
void Genotype::setRank(int rank)
{
    this->rank = uint16_t(rank);
}

/////////////////////////////////////////////////////////
void Genotype::insertToDominationSet(Genotype &i)
{
    Genotype* g_p = &i;
    this->dominates.push_back(g_p);
}

/////////////////////////////////////////////////////////
Genotype::Genotype()
{
    this->chromosome.resize(1,1);
    this->num_rows = 1;
    this->num_cols = 1;
    this->num_objectives = 2;
    this->objective_values.resize(this->num_objectives);
    this->domination_counter = 0;
    this->rank = 0;
    this->crowding_distance = 0;
    this->tot_segment_count = 0;
}

/////////////////////////////////////////////////////////
Genotype::Genotype(uint16_t num_rows, uint16_t num_cols)
{
    this->chromosome.resize(num_rows, num_cols);
    this->num_rows = num_rows;
    this->num_cols = num_cols;
    this->num_objectives = 2;
    this->objective_values.resize(this->num_objectives);
    this->domination_counter = 0;
    this->rank = 0;
    this->crowding_distance = 0;
    this->tot_segment_count = 0;
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
                //this->chromosome(row, col).value = genValues::none;

                /* Check if neighbouring node points at current node, if it does, point back to neighbour */
                if(row + 1 < num_rows) {
                    if (chromosome(row + 1, col).value == genValues::up) {
                        this->chromosome(row, col).value = genValues::down;
                    }
                }
                else if(row - 1 >= 0) {
                    if (chromosome(row - 1, col).value == genValues::down) {
                        this->chromosome(row, col).value = genValues::up;
                    }
                }
                else if(col + 1 < num_cols) {
                    if (chromosome(row, col + 1).value == genValues::left) {
                        this->chromosome(row, col).value = genValues::right;
                    }
                }
                else if(col - 1 >= 0) {
                    if (chromosome(row, col - 1).value == genValues::right) {
                        this->chromosome(row, col).value = genValues::left;
                    }
                }
                else{
                    unsigned seed = (unsigned) chrono::system_clock::now().time_since_epoch().count();
                    default_random_engine generator(seed);
                    uniform_int_distribution<uint8_t> gene_value_distribution(0, 3);
                    this->chromosome(row, col).value = (uint8_t) gene_value_distribution(generator);
                }
            }
            else{
                if (parents[i] - i == 1) {
                    this->chromosome(row, col).value = genValues::right;
                }
                else if (parents[i] - i == -1) {
                    this->chromosome(row, col).value = genValues::left;
                }
                else if (parents[i] - i == num_cols) {
                    this->chromosome(row, col).value = genValues::down;
                }
                else if (parents[i] - i == -num_cols){
                    this->chromosome(row, col).value = genValues::up;
                }
                else{
                    printf("Error in chromosome initialization\n");
                }
            }
        }
    }
}

/////////////////////////////////////////////////////////
void Genotype::setChromosomeValue(uint8_t value, int row, int col)
{
    if (value > -1 && value < 5)
        this->chromosome(row, col).value = value;
}

/////////////////////////////////////////////////////////
uint8_t Genotype::getChromosomeValue(int row, int col)
{
    return this->chromosome(row, col).value;
}

/////////////////////////////////////////////////////////
void Genotype::setChromosomeSegment(int segment, int row, int col)
{
    if (segment != -1)
        this->chromosome(row, col).segment = segment;
}

/////////////////////////////////////////////////////////
GeneNode* Genotype::getChromosomeGeneNode(int row, int col){
    return &this->chromosome(row, col);
}

/////////////////////////////////////////////////////////
bool Genotype::operator<(const Genotype &rhs) const
{
    if (this->objective_values[0] < rhs.objective_values[0] and
        this->objective_values[1] < rhs.objective_values[1]){
        return true;
    }
    else{
        return false;
    }
}

/////////////////////////////////////////////////////////
bool Genotype::operator>(const Genotype &rhs) const
{
    if (this->objective_values[0] > rhs.objective_values[0] and
        this->objective_values[1] > rhs.objective_values[1]){
        return true;
    }
    else{
        return false;
    }
}

/////////////////////////////////////////////////////////
void Genotype::genotypeToPhenotypeDecoding() {
    uint16_t segment_number = 0;
    int total_number_of_segments = segment_number;
    vector<tuple<int,int>> connected_nodes_pos;
    tuple<int, int> current_node_pos;
    Eigen::MatrixXi discovered = Eigen::MatrixXi::Zero(num_rows, num_cols);
    int current_node_row = -1, current_node_col = -1;
    GeneNode* current_node;

    for (int row = 0; row < num_rows; row++){
        for (int col = 0; col < num_cols; col++){
            if (!discovered(row,col)){
                /* Mark as discovered and assign segment number */
                current_node = &chromosome(row, col);
                current_node->segment = segment_number;
                discovered(row, col) = 1;

                /* Add it's connected neighbours to the list of nodes to do assigning */
                current_node_pos = make_tuple(row, col);
                findAndAddParentNodesToList(connected_nodes_pos, current_node_pos, discovered);
                addChildNodeToList(connected_nodes_pos, current_node_pos, discovered);

                /* Recursively do the same procedure until all members of the segment is found */
                while(!connected_nodes_pos.empty()){
                    current_node_row = get<0>(connected_nodes_pos.back());
                    current_node_col = get<1>(connected_nodes_pos.back());
                    current_node = &chromosome(current_node_row, current_node_col);
                    current_node->segment = segment_number;
                    connected_nodes_pos.pop_back();

                    current_node_pos = make_tuple(current_node_row, current_node_col);
                    findAndAddParentNodesToList(connected_nodes_pos, current_node_pos, discovered);
                    addChildNodeToList(connected_nodes_pos, current_node_pos, discovered);
                }
                segment_number++;
            }
        }
    }
    tot_segment_count = segment_number;
}

/////////////////////////////////////////////////////////
void Genotype::findAndAddParentNodesToList(vector<tuple<int,int>> &connected_nodes_pos, const tuple<int,int> &current_node_pos, Eigen::MatrixXi &discovered){
    int row = get<0>(current_node_pos);
    int col = get<1>(current_node_pos);

    /* Check if neighbouring node points at current node, if it does, its connected to the current node*/
    if(row + 1 < num_rows) {
        if (chromosome(row + 1, col).value == genValues::up && !discovered(row + 1, col)) {
            discovered(row + 1, col) = 1;
            tuple<int, int> new_connected_node_pos = make_tuple(row + 1, col);
            connected_nodes_pos.push_back(new_connected_node_pos);
        }
    }
    if(row - 1 >= 0) {
        if (chromosome(row - 1, col).value == genValues::down && !discovered(row - 1, col)) {
            discovered(row - 1, col) = 1;
            tuple<int, int> new_connected_node_pos = make_tuple(row - 1, col);
            connected_nodes_pos.push_back(new_connected_node_pos);
        }
    }
    if(col + 1 < num_cols) {
        if (chromosome(row, col + 1).value == genValues::left && !discovered(row, col + 1)) {
            discovered(row, col + 1) = 1;
            tuple<int, int> new_connected_node_pos = make_tuple(row, col + 1);
            connected_nodes_pos.push_back(new_connected_node_pos);
        }
    }
    if(col - 1 >= 0) {
        if (chromosome(row, col - 1).value == genValues::right && !discovered(row, col - 1)) {
            discovered(row, col - 1) = 1;
            tuple<int, int> new_connected_node_pos = make_tuple(row, col - 1);
            connected_nodes_pos.push_back(new_connected_node_pos);
        }
    }
}

/////////////////////////////////////////////////////////
void Genotype::addChildNodeToList(vector<tuple<int,int>> &connected_nodes_pos, const tuple<int,int> &current_node_pos, Eigen::MatrixXi &discovered){
    int row = get<0>(current_node_pos);
    int col = get<1>(current_node_pos);
    uint8_t nodeValue = chromosome(row,col).value;

    /* Check if child node (current node is poiting to) is discovered, if not, add to connected_nodes_pos */
    if(row + 1 < num_rows) {
        if (nodeValue == genValues::down && !discovered(row + 1, col)) {
            discovered(row + 1, col) = 1;
            tuple<int, int> new_connected_node_pos = make_tuple(row + 1, col);
            connected_nodes_pos.push_back(new_connected_node_pos);
            return;
        }
    }
    if(row - 1 >= 0) {
        if (nodeValue == genValues::up && !discovered(row - 1, col)) {
            discovered(row - 1, col) = 1;
            tuple<int, int> new_connected_node_pos = make_tuple(row - 1, col);
            connected_nodes_pos.push_back(new_connected_node_pos);
            return;
        }
    }
    if(col + 1 < num_cols) {
        if (nodeValue == genValues::right && !discovered(row, col + 1)) {
            discovered(row, col + 1) = 1;
            tuple<int, int> new_connected_node_pos = make_tuple(row, col + 1);
            connected_nodes_pos.push_back(new_connected_node_pos);
            return;
        }
    }
    if(col - 1 >= 0) {
        if (nodeValue == genValues::left && !discovered(row, col - 1)) {
            discovered(row, col - 1) = 1;
            tuple<int, int> new_connected_node_pos = make_tuple(row, col - 1);
            connected_nodes_pos.push_back(new_connected_node_pos);
            return;
        }
    }
}



/////////////////////////////////////////////////////////
void Genotype::visualizeSegments(const Eigen::MatrixXi &blue_ch, const Eigen::MatrixXi &green_ch,
                                 const Eigen::MatrixXi &red_ch) {
    Eigen::MatrixXi segment_eigen_matrix(num_rows, num_cols);
    cv::Mat segment_cv_image(num_rows, num_cols, CV_8UC3, cv::Scalar(0, 0, 0));

    vector<int> segments;
    vector<cv::Vec3b> segment_colors;
    cv::Vec3b current_segment_color;
    int current_segment_index = 0;
    int current_segment = 0;

    for (int i = 0; i < num_rows; i++) {
        for (int j = 0; j < num_cols; j++) {
            current_segment = this->chromosome(i, j).segment;
            auto iterator = find(segments.begin(), segments.end(), current_segment);
            if (iterator != segments.end()) {
                // current segment has been visited before
                current_segment_index = distance(segments.begin(), iterator);
                current_segment_color = segment_colors[current_segment_index];
            } else {
                // current segment is visited for the first time
                current_segment_color = cv::Vec3b(uint8_t(rand()*255), uint8_t(rand()*255),
                                                  uint8_t(rand()*255));
                segment_colors.push_back(current_segment_color);

                current_segment_index = segments.size();
                segments.push_back(current_segment);
            }

            segment_cv_image.at<cv::Vec3b>(i, j) = current_segment_color;
            segment_eigen_matrix(i, j) = current_segment;
        }
    }
    cv::namedWindow("Segments", cv::WINDOW_AUTOSIZE);
    cv::imshow("Segments", segment_cv_image);
    cv::waitKey(0);

}

/////////////////////////////////////////////////////////
void Genotype::decodeAndEvaluate(const Eigen::MatrixXi &red, const Eigen::MatrixXi &green, const Eigen::MatrixXi &blue){
    uint16_t segment_number = 0;
    int total_number_of_segments = segment_number;
    vector<tuple<int,int>> connected_nodes_pos;
    tuple<int, int> current_node_pos;
    Eigen::MatrixXi discovered = Eigen::MatrixXi::Zero(num_rows, num_cols);
    int current_node_row = -1, current_node_col = -1;
    GeneNode* current_node;
    int total_pixels_found = 0;

    vector<pixel_t> pixels_in_segment;
    objective_values[0] = 0;
    objective_values[1] = 0;
    int total_edge_pixels = 0;
    for (int row = 0; row < num_rows; row++){
        for (int col = 0; col < num_cols; col++){
            if (discovered(row,col) == 0){
                /* Mark as discovered and assign segment number */
                current_node = &chromosome(row, col);
                current_node->segment = segment_number;
                discovered(row, col) = 1;

                /*  Calculate the rgb centroid of the segment */
                rgb_centroid_t centroid;
                int num_pixels = 0;
                centroid.r += red(row, col);
                centroid.g += green(row, col);
                centroid.b += blue(row, col);
                num_pixels++;
                total_pixels_found++;

                /* Keep track of pixels in current segment */
                pixel_t pixel(row, col);
                pixels_in_segment.push_back(pixel);

                /* Add its connected neighbours to the list of nodes to do assigning */
                current_node_pos = make_tuple(row, col);
                findAndAddParentNodesToList(connected_nodes_pos, current_node_pos, discovered);
                addChildNodeToList(connected_nodes_pos, current_node_pos, discovered);

                /* Recursively do the same procedure until all members of the segment is found */
                while(!connected_nodes_pos.empty()){
                    current_node_row = get<0>(connected_nodes_pos.back());
                    current_node_col = get<1>(connected_nodes_pos.back());
                    current_node = &chromosome(current_node_row, current_node_col);
                    current_node->segment = segment_number;
                    connected_nodes_pos.pop_back();

                    centroid.r += red(current_node_row, current_node_col);
                    centroid.g += green(current_node_row, current_node_col);
                    centroid.b += blue(current_node_row, current_node_col);
                    num_pixels++;
                    total_pixels_found++;

                    pixel_t pixel(current_node_row, current_node_col);
                    pixels_in_segment.push_back(pixel);

                    current_node_pos = make_tuple(current_node_row, current_node_col);
                    findAndAddParentNodesToList(connected_nodes_pos, current_node_pos, discovered);
                    addChildNodeToList(connected_nodes_pos, current_node_pos, discovered);
                }
                centroid.r = centroid.r / num_pixels;
                centroid.g = centroid.g / num_pixels;
                centroid.b = centroid.b / num_pixels;

                /* Calculate objective values */
                int segment_pixel_centroid_deviation = 0;
                for(auto &p: pixels_in_segment){
                    // objective 1
                    double pixel_centroid_deviation = sqrt( pow(red(p.row, p.col) - centroid.r, 2)
                                                            + pow(green(p.row, p.col) - centroid.g, 2)
                                                            + pow(blue(p.row, p.col) - centroid.b, 2));
                    segment_pixel_centroid_deviation += pixel_centroid_deviation;
                    // objective 2
                    int segment_num = chromosome(p.row, p.col).segment;
                    double segment_boundary_diff = 0;
                    segment_boundary_diff -= calcEuclideanRgbDiff(-1, 0, p.col, p.row, segment_num, red,green,blue);
                    segment_boundary_diff -= calcEuclideanRgbDiff(1, 0, p.col, p.row, segment_num, red,green,blue);
                    segment_boundary_diff -= calcEuclideanRgbDiff(0, -1, p.col, p.row, segment_num, red,green,blue);
                    segment_boundary_diff -= calcEuclideanRgbDiff(0, 1, p.col, p.row, segment_num, red,green,blue);
                    objective_values[1] += segment_boundary_diff;
                    if (segment_boundary_diff < 0){ total_edge_pixels++; }
                }
                objective_values[0] += (segment_pixel_centroid_deviation);// / num_pixels);
                segment_number++;
                pixels_in_segment.clear();
                if (total_pixels_found == num_cols*num_rows){
                    tot_segment_count = segment_number;
                    //objective_values[1] = objective_values[1]/total_edge_pixels;
                    return;
                }
            }
        }
    }
    //objective_values[1] = objective_values[1]/total_edge_pixels;
    tot_segment_count = segment_number;
}

/////////////////////////////////////////////////////////
void Genotype::visualizeEdges(cv::Mat test_image, string title)
{
    // create white image
    cv::Mat segment_cv_image(num_rows, num_cols, CV_8UC1, cv::Scalar(255));

    // create grouping of pixels in corresponding segments
    vector<vector<pixel_t>> pixels_segment_affiliation (tot_segment_count);
    for (int row = 0; row < num_rows; row++) {
        for (int col = 0; col < num_cols; col++) {
            int segment_num = chromosome(row,col).segment;
            pixel_t pixel = {row, col};
            pixels_segment_affiliation[segment_num].push_back(pixel);
        }
    }
    for (const auto &s: pixels_segment_affiliation) {
        for (const auto &p: s){
            // edge pixel?
            int this_row = p.row;
            int this_col = p.col;
            int this_segment = chromosome(this_row, this_col).segment;

            if (isEdgePixel(1, 0, this_col, this_row, this_segment)){
                segment_cv_image.at<uchar>(this_row, this_col) = 0;
            }
            else if (isEdgePixel(-1, 0, this_col, this_row, this_segment)){
                segment_cv_image.at<uchar>(this_row, this_col) = 0;
            }
            else if (isEdgePixel(0, 1, this_col, this_row, this_segment)){
                segment_cv_image.at<uchar>(this_row, this_col) = 0;
            }
            else if (isEdgePixel(0, -1, this_col, this_row, this_segment)){
                segment_cv_image.at<uchar>(this_row, this_col) = 0;
            }
        }
    }

    // invert the edge image
    cv::Mat inverted_image(num_rows, num_cols, CV_8UC1, cv::Scalar(0));
    cv::bitwise_not(segment_cv_image, inverted_image);

    // thinner edges
    thinning(inverted_image, inverted_image);
    cv::bitwise_not(inverted_image, segment_cv_image);
//    cv::namedWindow("Segments 2.0", cv::WINDOW_AUTOSIZE);
//    cv::imshow("Segments 2.0", segment_cv_image);
//    cv::waitKey(0);

    // create copy and draw green edges in original image
    cv::Mat copy = test_image.clone();
    cv::Mat green_image(num_rows, num_cols, CV_8UC3, cv::Scalar(0, 255, 0));
    green_image.copyTo(copy, inverted_image);
//    cv::namedWindow("Green", cv::WINDOW_AUTOSIZE);
//    cv::imshow("Green", copy);
//    cv::waitKey(0);

    // Create mat for window
    int width = copy.cols;
    int height = copy.rows;
    cv::Mat win_mat(cv::Size(width*3, height), CV_8UC3);
    cv::Mat segment_cv_image_color(num_rows, num_cols, CV_8UC3);
    cv::cvtColor(segment_cv_image, segment_cv_image_color, cv::COLOR_GRAY2BGR);

    // Copy small images into big mat
    test_image.copyTo(win_mat(cv::Rect(  0, 0, width, height)));
    copy.copyTo(win_mat(cv::Rect(width, 0, width, height)));
    segment_cv_image_color.copyTo(win_mat(cv::Rect(width*2, 0, width, height)));

    // Display big mat
    cv::imshow(title, win_mat);
    cv::waitKey(0);
}



/////////////////////////////////////////////////////////
void Genotype::calculateObjectives(const Eigen::MatrixXi &red, const Eigen::MatrixXi &green, const Eigen::MatrixXi &blue) {
    // create grouping of pixels in corresponding segments
    vector<vector<pixel_t>> pixels_segment_affiliation (tot_segment_count);
    for (int row = 0; row < num_rows; row++) {
        for (int col = 0; col < num_cols; col++) {
            int segment_num = chromosome(row,col).segment;
            pixel_t pixel = {row, col};
            pixels_segment_affiliation[segment_num].push_back(pixel);
        }
    }
    // calculate the rgb centroid for each segment
    vector<rgb_centroid_t> centroids (tot_segment_count);
    objective_values[0] = 0;
    objective_values[1] = 0;
    int i = 0;
    for (const auto &s: pixels_segment_affiliation){
        rgb_centroid_t centroid;
        double intensity_sum_r = 0;
        double intensity_sum_g = 0;
        double intensity_sum_b = 0;
        int num_pixels = 0;
        for (const auto &p: s){
            intensity_sum_r += red(p.row, p.col);
            intensity_sum_g += green(p.row, p.col);
            intensity_sum_b += blue(p.row, p.col);
            num_pixels++;
        }
        centroid.r = intensity_sum_r / num_pixels;
        centroid.g = intensity_sum_g / num_pixels;
        centroid.b = intensity_sum_b / num_pixels;
        centroids[i] = centroid;
        for (const auto &p: s){
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
            segment_boundary_diff -= calcEuclideanRgbDiff(-1, 0, this_col, this_row, this_segment, red,green,blue);
            segment_boundary_diff -= calcEuclideanRgbDiff(1, 0, this_col, this_row, this_segment, red,green,blue);
            segment_boundary_diff -= calcEuclideanRgbDiff(0, -1, this_col, this_row, this_segment, red,green,blue);
            segment_boundary_diff -= calcEuclideanRgbDiff(0, 1, this_col, this_row, this_segment, red,green,blue);
            objective_values[1] += segment_boundary_diff;
        }
        i++;
    }
    // calculate objective 1 and 2

//    i = 0;
//    for (const auto &s: pixels_segment_affiliation) {
//        for (const auto &p: s){
//            // objective 1
//            double pixel_centroid_deviation = sqrt( pow(red(p.row, p.col) - centroids[i].r, 2)
//                                                    + pow(green(p.row, p.col) - centroids[i].g, 2)
//                                                    + pow(blue(p.row, p.col) - centroids[i].b, 2));
//            objective_values[0] += pixel_centroid_deviation;
//            // objective 2
//            int this_row = p.row;
//            int this_col = p.col;
//            int this_segment = chromosome(this_row, this_col).segment;
//            double segment_boundary_diff = 0;
//            segment_boundary_diff -= calcEuclideanRgbDiff(-1, 0, this_col, this_row, this_segment, red,green,blue);
//            segment_boundary_diff -= calcEuclideanRgbDiff(1, 0, this_col, this_row, this_segment, red,green,blue);
//            segment_boundary_diff -= calcEuclideanRgbDiff(0, -1, this_col, this_row, this_segment, red,green,blue);
//            segment_boundary_diff -= calcEuclideanRgbDiff(0, 1, this_col, this_row, this_segment, red,green,blue);
//            objective_values[1] += segment_boundary_diff;
//            }
//        i++;
//    }
}

/////////////////////////////////////////////////////////
bool Genotype::sortByObj1(const Genotype* lhs, const Genotype* rhs) { return lhs->objective_values[0] < rhs->objective_values[0]; }

/////////////////////////////////////////////////////////
bool Genotype::sortByObj2(const Genotype* lhs, const Genotype* rhs) { return lhs->objective_values[1] < rhs->objective_values[1]; }

/////////////////////////////////////////////////////////
bool Genotype::sortByCrowdedComparison(const Genotype* lhs, const Genotype* rhs) {
    if (lhs->rank != rhs->rank) {
        return lhs->rank < rhs->rank;
    }
    else if (lhs->rank == rhs->rank)
    {
        return lhs->crowding_distance > rhs->crowding_distance;
    }
}

/////////////////////////////////////////////////////////
double Genotype::calcEuclideanRgbDiff(signed short dir_y, signed short dir_x, int this_col, int this_row, int this_segment,
                                      const Eigen::MatrixXi &red, const Eigen::MatrixXi &green, const Eigen::MatrixXi &blue) {
    double diff = 0;
    if (this_col + dir_x >= 0 && this_col + dir_x < num_cols && this_row + dir_y >= 0 && this_row + dir_y < num_rows) {
        if (chromosome(this_row + dir_y, this_col + dir_x).segment != this_segment) {
            diff = sqrt(pow(red(this_row, this_col) - red(this_row + dir_y, this_col + dir_x), 2)
                        + pow(green(this_row, this_col) - green(this_row + dir_y, this_col + dir_x), 2)
                        + pow(blue(this_row, this_col) - blue(this_row + dir_y, this_col + dir_x), 2));

        }
    }
    return diff;
}

/////////////////////////////////////////////////////////
bool Genotype::isEdgePixel(signed short dir_y, signed short dir_x, int this_col, int this_row, int this_segment)
{
    if (this_col + dir_x >= 0 && this_col + dir_x < num_cols && this_row + dir_y >= 0 && this_row + dir_y < num_rows) {
        int segment = chromosome(this_row + dir_y, this_col + dir_x).segment;
        return (chromosome(this_row + dir_y, this_col + dir_x).segment != this_segment);
    }
    else{
        return false;
    }
}


void Genotype::mergeSegments(const Eigen::MatrixXi &red, const Eigen::MatrixXi &green, const Eigen::MatrixXi &blue)
{
    // assume decoding has been done
    // create grouping of pixels in corresponding segments
    vector<vector<pixel_t>> pixels_segment_affiliation (tot_segment_count);
    for (int row = 0; row < num_rows; row++) {
        for (int col = 0; col < num_cols; col++) {
            int segment_num = chromosome(row,col).segment;
            pixel_t pixel = {row, col};
            pixels_segment_affiliation[segment_num].push_back(pixel);
        }
    }
    int smallest_segment_size = INT_MAX;
    int smallest_segment = 0;
    int i_segment = 0;
    vector<rgb_centroid_t> centroids(tot_segment_count);
    for (auto &s: pixels_segment_affiliation){
        if (s.size() < smallest_segment_size){
            smallest_segment_size = s.size();
            smallest_segment =  i_segment;
        }
        /*  Calculate the rgb centroid of the segment */
        int num_pixels = 0;
        for (const auto &p: s){
            centroids[i_segment].r += red(p.row, p.col);
            centroids[i_segment].g += green(p.row, p.col);
            centroids[i_segment].b += blue(p.row, p.col);
            num_pixels++;
        }
        centroids[i_segment].r = centroids[i_segment].r / num_pixels;
        centroids[i_segment].g = centroids[i_segment].g / num_pixels;
        centroids[i_segment].b = centroids[i_segment].b / num_pixels;
        i_segment++;
    }

    int MIN_SEGMENT_SIZE = 700; //test
    int MAX_SEGMENT_NUM = 15;
    while (smallest_segment_size < MIN_SEGMENT_SIZE || tot_segment_count > MAX_SEGMENT_NUM){
        // choose smallest segment
        int i_current_segment = smallest_segment;
        rgb_centroid_t current_centroid = centroids[i_current_segment];
        vector<pixel_t> current_segment = pixels_segment_affiliation[i_current_segment];

        // Find neighbouring segments
        set<int16_t> neighbour_segments;
        int smallest_neighbour_size = INT_MAX;
        vector<pixel_t> smallest_neighbour_segment;
        int smallest_neighbour_i = 0;

        double smallest_centroid_diff = DBL_MAX;
        vector <pixel_t> smallest_centroid_diff_segment;
        int smallest_centroid_diff_i = 0;
        int smallest_centroid_diff_num_pixels = 0;

        for (const auto &p: current_segment){
            int this_row = p.row;
            int this_col = p.col;
            int16_t this_segment = chromosome(this_row, this_col).segment;

            if (isEdgePixel(1, 0, this_col, this_row, this_segment)){
                if (neighbour_segments.find(chromosome(this_row + 1, this_col).segment) == neighbour_segments.end()){
                    neighbour_segments.insert(chromosome(this_row + 1, this_col).segment);
                    // find size of neighbouring segment, keep track of smallest neighbour
                    // find centroid of neighbouring segment, keep track of most similar centroid
                    int i_s = -1;
                    for (const auto &s: pixels_segment_affiliation){
                        i_s++;
                        if (chromosome(this_row + 1, this_col).segment == chromosome(s[0].row, s[0].col).segment) {
//                            if (s.size() < smallest_neighbour_size) {
//                                smallest_neighbour_size = s.size();
//                                smallest_neighbour_segment = s;
//                                smallest_neighbour_i = i_s;
//                            }
                            double diff = sqrt(pow(current_centroid.r - centroids[i_s].r, 2)
                                  + pow(current_centroid.g - centroids[i_s].g, 2)
                                  + pow(current_centroid.b - centroids[i_s].b, 2));
                            if (diff < smallest_centroid_diff){
                                smallest_centroid_diff = diff;
                                smallest_centroid_diff_segment = s;
                                smallest_centroid_diff_i = i_s;
                            }
                        }
                    }
                }
            }
            else if (isEdgePixel(-1, 0, this_col, this_row, this_segment)){
                if (neighbour_segments.find(chromosome(this_row - 1, this_col).segment) == neighbour_segments.end()) {
                    neighbour_segments.insert(chromosome(this_row - 1, this_col).segment);

                    // find size of neighbouring segment, keep track of smallest neighbour
                    int i_s = -1;
                    for (const auto &s: pixels_segment_affiliation){
                        i_s++;
                        if (chromosome(this_row - 1, this_col).segment == chromosome(s[0].row, s[0].col).segment) {
//                            if (s.size() < smallest_neighbour_size) {
//                                smallest_neighbour_size = s.size();
//                                smallest_neighbour_segment = s;
//                                smallest_neighbour_i = i_s;
//                            }
                            double diff = sqrt(pow(current_centroid.r - centroids[i_s].r, 2)
                                               + pow(current_centroid.g - centroids[i_s].g, 2)
                                               + pow(current_centroid.b - centroids[i_s].b, 2));
                            if (diff < smallest_centroid_diff){
                                smallest_centroid_diff = diff;
                                smallest_centroid_diff_segment = s;
                                smallest_centroid_diff_i = i_s;
                            }
                        }
                    }
                }
            }
            else if (isEdgePixel(0, 1, this_col, this_row, this_segment)){
                if (neighbour_segments.find(chromosome(this_row, this_col + 1).segment) == neighbour_segments.end()) {
                    neighbour_segments.insert(chromosome(this_row, this_col + 1).segment);
                    // find size of neighbouring segment, keep track of smallest neighbour
                    int i_s = -1;
                    for (const auto &s: pixels_segment_affiliation){
                        i_s++;
                        if (chromosome(this_row, this_col + 1).segment == chromosome(s[0].row, s[0].col).segment) {
//                            if (s.size() < smallest_neighbour_size) {
//                                smallest_neighbour_size = s.size();
//                                smallest_neighbour_segment = s;
//                                smallest_neighbour_i = i_s;
//                            }
                            double diff = sqrt(pow(current_centroid.r - centroids[i_s].r, 2)
                                               + pow(current_centroid.g - centroids[i_s].g, 2)
                                               + pow(current_centroid.b - centroids[i_s].b, 2));
                            if (diff < smallest_centroid_diff){
                                smallest_centroid_diff = diff;
                                smallest_centroid_diff_segment = s;
                                smallest_centroid_diff_i = i_s;
                            }
                        }
                    }
                }
            }
            else if (isEdgePixel(0, -1, this_col, this_row, this_segment)){
                if (neighbour_segments.find(chromosome(this_row, this_col - 1).segment) == neighbour_segments.end()) {
                    neighbour_segments.insert(chromosome(this_row, this_col - 1).segment);
                    // find size of neighbouring segment, keep track of smallest neighbour
                    int i_s = -1;
                    for (const auto &s: pixels_segment_affiliation){
                        i_s++;
                        if (chromosome(this_row, this_col - 1).segment == chromosome(s[0].row, s[0].col).segment) {
//                            if (s.size() < smallest_neighbour_size) {
//                                smallest_neighbour_size = s.size();
//                                smallest_neighbour_segment = s;
//                                smallest_neighbour_i = i_s;
//                            }
                            double diff = sqrt(pow(current_centroid.r - centroids[i_s].r, 2)
                                               + pow(current_centroid.g - centroids[i_s].g, 2)
                                               + pow(current_centroid.b - centroids[i_s].b, 2));
                            if (diff < smallest_centroid_diff){
                                smallest_centroid_diff = diff;
                                smallest_centroid_diff_segment = s;
                                smallest_centroid_diff_i = i_s;
                            }
                        }
                    }
                }
            }
        }

//        // Smallest neighbour segment found: merge
//        if (!neighbour_segments.empty()){
//            int16_t current_segment_num = chromosome(current_segment[0].row, current_segment[0].col).segment;
//            for (auto &p: smallest_neighbour_segment){
//                chromosome(p.row, p.col).segment = current_segment_num;
//                pixels_segment_affiliation[i_current_segment].push_back(p);
//            }
//            pixels_segment_affiliation.erase(pixels_segment_affiliation.begin() + smallest_neighbour_i);
//            tot_segment_count--;
//        }

        // Smallest centroid diff segment found: merge
        if (!neighbour_segments.empty()){
            int16_t current_segment_num = chromosome(current_segment[0].row, current_segment[0].col).segment;
            centroids[i_current_segment].r = centroids[i_current_segment].r * current_segment.size();
            centroids[i_current_segment].g = centroids[i_current_segment].g * current_segment.size();
            centroids[i_current_segment].b = centroids[i_current_segment].b * current_segment.size();
            for (auto &p: smallest_centroid_diff_segment){
                chromosome(p.row, p.col).segment = current_segment_num;
                pixels_segment_affiliation[i_current_segment].push_back(p);
                centroids[i_current_segment].r += red(p.row, p.col);
                centroids[i_current_segment].g += green(p.row, p.col);
                centroids[i_current_segment].b += blue(p.row, p.col);
            }
            centroids[i_current_segment].r = centroids[i_current_segment].r / (smallest_centroid_diff_segment.size() + current_segment.size());
            centroids[i_current_segment].g = centroids[i_current_segment].g / (smallest_centroid_diff_segment.size() + current_segment.size());
            centroids[i_current_segment].b = centroids[i_current_segment].b / (smallest_centroid_diff_segment.size() + current_segment.size());
            pixels_segment_affiliation.erase(pixels_segment_affiliation.begin() + smallest_centroid_diff_i);
            centroids.erase(centroids.begin() + smallest_centroid_diff_i);
            tot_segment_count--;
        }

        // update centroids
        // find new smallest segment
        smallest_segment_size = INT_MAX;
        smallest_segment = 0;
        i_segment = 0;
        //centroids.resize(tot_segment_count);
        for (const auto &s: pixels_segment_affiliation){
            if (s.size() < smallest_segment_size){
                smallest_segment_size = s.size();
                smallest_segment =  i_segment;//chromosome(s[0].row, s[0].col).segment;
            }
//            /*  Calculate the rgb centroid of the segment */
//            centroids[i_segment].r = 0, centroids[i_segment].g = 0, centroids[i_segment].b = 0;
//            int num_pixels = 0;
//            for (const auto &p: s){
//                centroids[i_segment].r += red(p.row, p.col);
//                centroids[i_segment].g += green(p.row, p.col);
//                centroids[i_segment].b += blue(p.row, p.col);
//                num_pixels++;
//            }
//            centroids[i_segment].r = centroids[i_segment].r / num_pixels;
//            centroids[i_segment].g = centroids[i_segment].g / num_pixels;
//            centroids[i_segment].b = centroids[i_segment].b / num_pixels;
            i_segment++;
        }
    }

    // Update segment numbers in the range [0, tot_num_segments)
    i_segment = 0;
    for (const auto &s: pixels_segment_affiliation){
        for (const auto &p: s){
            chromosome(p.row, p.col).segment = i_segment;
        }
        i_segment++;
    }
}

/*
Genotype::~Genotype(){
    cout << "deleted genotype!" << endl;
}*/
