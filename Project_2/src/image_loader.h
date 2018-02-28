#ifndef PROJECT_2_IMAGE_LOADER_H
#define PROJECT_2_IMAGE_LOADER_H

#include <opencv2/core.hpp>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <opencv/cxeigen.hpp>

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

#endif //PROJECT_2_IMAGE_LOADER_H
