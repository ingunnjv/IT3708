#pragma once
#include "genotype.h"

#include <opencv2/imgproc.hpp>
#include <opencv2/core/eigen.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgcodecs.hpp>

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
            this->chromosome(row, col).segment = -1;
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

void Genotype::genotypeToPhenotypeDecoding(int num_rows, int num_cols)
{
    int segment_number = 0;
    int total_number_of_segments = segment_number;

    for(int i = 0; i < this->chromosome.size(); i++) {
        int row = i / num_cols, col = i % num_cols;

        if (this->chromosome(row, col).segment == -1) { // pixel is not yet assigned to a segment
            vector<GeneNode*> discovery_list;
            segment_number = total_number_of_segments;
            this->chromosome(row, col).segment = segment_number;
            total_number_of_segments++;

            if (this->chromosome(row, col).child != NULL && this->chromosome(row, col).child->segment == -1)
                discovery_list.push_back(this->chromosome(row, col).child);
            if (!this->chromosome(row, col).parents.empty()) {
                for (auto p: this->chromosome(row, col).parents){
                    if (p->segment == -1)
                        discovery_list.push_back(p);
                }
            }

            while(!discovery_list.empty()){
                GeneNode* current_gene = discovery_list.front();
                discovery_list.erase(discovery_list.begin());

                current_gene->segment = segment_number;
                if (current_gene->child != NULL && current_gene->child->segment == -1)
                    discovery_list.push_back(current_gene->child);
                if (!current_gene->parents.empty()) {
                    for (auto p: current_gene->parents){
                        if (p->segment == -1)
                            discovery_list.push_back(p);
                    }
                }
            }
        }
    }

    //// print image segments
    //for (int i = 0; i < num_rows; i++){
    //    for (int j = 0; j < num_cols; j++){
    //        cout << this->chromosome(i, j).segment;
    //    }
    //    cout << endl;
    //}

    cout << "Total number of segments: " << total_number_of_segments << endl;
    /*
    // create vector of segments
    this->segments.resize(total_number_of_segments);
    for (int pixel = 0; pixel < this->chromosome.size(); pixel++){
        int row = pixel / num_cols, col = pixel % num_cols;
        this->segments[segmented_image(row, col)].push_back(pixel);
    }
                 */
}

void Genotype::visualize(cv::Mat &test_image, int num_rows, int num_cols)
{
    Eigen::MatrixXi segment_eigen_image(num_rows, num_cols);
    cv::Mat segment_cv_image(num_rows, num_cols, CV_8UC1, cv::Scalar(0));
    for (int i = 0; i < num_rows; i++){
        for (int j = 0; j < num_cols; j++){
            segment_cv_image.at<char>(i, j) = uint8_t(this->chromosome(i, j).segment);
            segment_eigen_image(i, j) = this->chromosome(i, j).segment;

        }
    }
    cout << segment_eigen_image << endl;


    cv::Mat canny_output;
    vector<vector<cv::Point> > contours;
    vector<cv::Vec4i> hierarchy;
    cv::Canny( segment_cv_image, canny_output, 0, 1, 3 );
    cv::findContours( canny_output, contours, hierarchy, cv::RETR_TREE, cv::CHAIN_APPROX_SIMPLE, cv::Point(0, 0) );

    cv::Mat white_drawing = cv::Mat::ones( canny_output.size(), CV_8UC1 ) * 255;
    cv::Mat test_image_drawing(test_image);
    for( size_t i = 0; i< contours.size(); i++ )
    {
        cv::Scalar color = cv::Scalar( 0 );
        cv::drawContours( white_drawing, contours, (int)i, color, 1, 1, hierarchy, 0, cv::Point() );

        color = cv::Scalar( 0, 200, 0 );
        cv::drawContours( test_image_drawing, contours, (int)i, color, 1, 1, hierarchy, 0, cv::Point() );
    }
    //cv::namedWindow( "Solution type 2", cv::WINDOW_AUTOSIZE );
    //cv::imshow( "Solution type 2", white_drawing );
    //cv::namedWindow( "Solution type 1", cv::WINDOW_AUTOSIZE );
    //cv::imshow( "Solution type 1", test_image_drawing );
    //cv::waitKey(0); // Wait for a keystroke in the window
}



