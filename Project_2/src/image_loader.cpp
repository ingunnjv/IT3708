//
// Created by LasseBot on 21-Feb-18.
//
#include "image_loader.h"
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <iostream>

using namespace std;


ImageLoader::ImageLoader()
{

}

void ImageLoader::LoadImagesFromFolder(string imagefolder)
{
    string basefolder = "../Test Images/";

    string filetypes[] = {"/GT*.jpg", "/s*.jpg", "/T*.jpg"};
    for (int i = 0; i < 3; i++){
        string pathname = basefolder + imagefolder + filetypes[i];
        cv::String path(pathname);
        vector<cv::String> fn;
        cv::glob(path, fn, true); // recurse
        for (size_t k=0; k<fn.size(); ++k)
        {
            cv::Mat im = cv::imread(fn[k], cv::IMREAD_COLOR);
            if (im.empty()) continue; //only proceed if successful

            // you probably want to do some preprocessing
            if (i == 0){ gt_images.push_back(im); }
            else if (i == 1){ segment_images.push_back(im); }
            else {Test_image = im;}
        }
    }
}

void ImageLoader::ExtractRGBChannels()
{
    cv::Mat planes[3];
    cv::split(Test_image, planes);  // BGR: planes[2] is the red channel
    b_image = planes[0];
    g_image = planes[1];
    r_image = planes[2];
}
