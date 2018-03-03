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

void Genotype::visualize(Eigen::MatrixXi &blue_ch, Eigen::MatrixXi &green_ch, Eigen::MatrixXi &red_ch, int num_rows, int num_cols)
{
    Eigen::MatrixXi segment_eigen_image(num_rows, num_cols);
    cv::Mat segment_cv_image(num_rows, num_cols, CV_8UC3, cv::Scalar(0, 0, 0));
    int current_segment = 0, prev_segment = 0;
    cv::Vec3b segment_color = cv::Vec3b(uint8_t(blue_ch(0, 0)), uint8_t(green_ch(0, 0)), uint8_t(red_ch(0, 0)));
    for (int i = 0; i < num_rows; i++){
        for (int j = 0; j < num_cols; j++){
            current_segment = this->chromosome(i, j).segment;
            if (current_segment != prev_segment){
                prev_segment = current_segment;
                segment_color = cv::Vec3b(uint8_t(blue_ch(i, j)), uint8_t(green_ch(i, j)), uint8_t(red_ch(i, j)));
            }

            segment_cv_image.at<cv::Vec3b>(i, j) = segment_color;
            segment_eigen_image(i, j) = this->chromosome(i, j).segment;
        }
    }
    cout << segment_eigen_image << endl;
    cv::namedWindow( "Segments", cv::WINDOW_AUTOSIZE );
    cv::imshow( "Segments", segment_cv_image );
    cv::waitKey(0);
    cv::Mat canny_output;
    vector<vector<cv::Point> > contours;
    vector<cv::Vec4i> hierarchy;
    cv::Canny( segment_cv_image, canny_output, 0, 0, 3 );
    cv::findContours( canny_output, contours, hierarchy, cv::RETR_TREE, cv::CHAIN_APPROX_SIMPLE, cv::Point(0, 0) );

    cv::Mat white_drawing = cv::Mat::ones( canny_output.size(), CV_8UC1 ) * 255;
    //cv::Mat test_image_drawing(test_image);
    cv::Mat segments = cv::Mat::ones( canny_output.size(), CV_8UC3 ) * 100;
    int red = 200, green = 100, blue = 0;
    int offset = 50;
    for( size_t i = 0; i< contours.size(); i++ )
    {
        if (i % 2 == 0 && i > 0){
            if (red + offset < 255)
                red = red + offset;
            else
                red = 0;
        }
        if (i % 3 == 0 && i > 0){
            if (green + offset < 255)
                green = green + offset;
            else
                green = 0;
        }
        if (i % 2 == 1 && i > 0){
            if (blue + offset/2 < 255)
                blue = blue + offset/2;
            else
                blue = 0;
        }

        cv::Scalar segment_color = cv::Scalar(blue, green, red);
        cv::drawContours( segments, contours, (int)i, segment_color, -1, 8, hierarchy, 0, cv::Point() );


        cv::Scalar color = cv::Scalar( 0 );
        cv::drawContours( white_drawing, contours, (int)i, color, 1, 1, hierarchy, 0, cv::Point() );

        color = cv::Scalar( 0, 200, 0 );
        //cv::drawContours( test_image_drawing, contours, (int)i, color, 1, 1, hierarchy, 0, cv::Point() );
    }
    cv::namedWindow( "Segments", cv::WINDOW_AUTOSIZE );
    cv::imshow( "Segments", segments );
    cv::namedWindow( "Solution type 2", cv::WINDOW_AUTOSIZE );
    cv::imshow( "Solution type 2", white_drawing );
    //cv::namedWindow( "Solution type 1", cv::WINDOW_AUTOSIZE );
    //cv::imshow( "Solution type 1", test_image_drawing );
    cv::waitKey(0); // Wait for a keystroke in the window
}



