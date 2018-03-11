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

public:
    cv::Mat test_image;
    cv::Mat nonmodified_image;
    Eigen::MatrixXi b_channel;
    Eigen::MatrixXi g_channel;
    Eigen::MatrixXi r_channel;
    ImageLoader();
    void loadImagesFromFolder(std::string imagefolder);
    void extractRGBChannels();

    int segmentation(std::string imagefolder);
    void kMeansClustering(std::string image_folder);

};

#endif //PROJECT_2_IMAGE_LOADER_H
