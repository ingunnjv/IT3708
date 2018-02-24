//
// Created by LasseBot on 21-Feb-18.
//

#ifndef PROJECT_2_IMAGE_LOADER_H
#define PROJECT_2_IMAGE_LOADER_H

#endif //PROJECT_2_IMAGE_LOADER_H
#include <opencv2/core.hpp>
#include <vector>

class ImageLoader
{
private:
    double A;
    std::vector<cv::Mat> GT_images;
    std::vector<cv::Mat> Segment_images;
    cv::Mat Test_image;



public:
    cv::Mat B_image;
    cv::Mat G_image;
    cv::Mat R_image;
    ImageLoader();
    void LoadImagesFromFolder(std::string imagefolder);
    void ExtractRGBChannels();
    double B;



};
