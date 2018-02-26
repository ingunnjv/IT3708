//
// Created by LasseBot on 21-Feb-18.
//

#ifndef PROJECT_2_IMAGE_LOADER_H
#define PROJECT_2_IMAGE_LOADER_H

#endif //PROJECT_2_IMAGE_LOADER_H
#pragma once
#include <opencv2/core.hpp>
#include <vector>


class ImageLoader
{
private:
    std::vector<cv::Mat> gt_images;
    std::vector<cv::Mat> segment_images;
    cv::Mat Test_image;

public:
    cv::Mat b_image;
    cv::Mat g_image;
    cv::Mat r_image;
    ImageLoader();
    void LoadImagesFromFolder(std::string imagefolder);
    void ExtractRGBChannels();

};
