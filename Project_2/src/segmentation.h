//
// Created by Ingunn on 11.03.2018.
//

#ifndef PROJECT_2_SEGMENTATION_H
#define PROJECT_2_SEGMENTATION_H
#include <vector>
#include <opencv2/core.hpp>

float weight(int p1[],int p2[]);
int root (std::vector <int> &Arr ,int i);
int segmentation(cv::Mat &image);

#endif //PROJECT_2_SEGMENTATION_H
