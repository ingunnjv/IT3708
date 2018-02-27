//
// Created by LasseBot on 21-Feb-18.
//

#pragma once
#ifndef PROJECT_2_IMAGE_LOADER_H
#define PROJECT_2_IMAGE_LOADER_H

#endif //PROJECT_2_IMAGE_LOADER_H

#include <opencv2/core.hpp>

class ImageLoader
{
private:
    std::vector<cv::Mat> gt_images;
    std::vector<cv::Mat> segment_images;
    cv::Mat test_image;

public:
    Eigen::MatrixXi b_channel;
    Eigen::MatrixXi g_channel;
    Eigen::MatrixXi r_channel;
    ImageLoader();
    void LoadImagesFromFolder(std::string imagefolder);
    void ExtractRGBChannels();

};
