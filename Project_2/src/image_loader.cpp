#pragma once
#include "image_loader.h"

#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

ImageLoader::ImageLoader()
{

}

void ImageLoader::loadImagesFromFolder(string imagefolder, uint8_t data_aug)
{
    string basefolder = "../Test Images/";

    string filetypes[] = {"/GT*.jpg", "/s*.jpg", "/T*.jpg"};
    for (int i = 0; i < 3; i++){
        string pathname = basefolder + imagefolder + filetypes[i];
        cv::String path(pathname);
        vector<cv::String> fn;
        cv::glob(path, fn, true); // recurse
        for (size_t k=0; k<fn.size(); ++k){
            cv::Mat im = cv::imread(fn[k], cv::IMREAD_COLOR);
            if (im.empty()) continue; //only proceed if successful

            // you probably want to do some preprocessing
            if (i == 0){ gt_images.push_back(im); }
            else if (i == 1){ segment_images.push_back(im); }
            else {test_image = im; nonmodified_image = im;}
        }
    }
    if(data_aug){
        kMeansClustering(2, imagefolder);
    }
}

void ImageLoader::extractRGBChannels()
{
    cv::Mat planes[3];
    cv::Mat normalized_planes[3];
    cv::split(test_image, planes);  // BGR: planes[2] is the red channel

    for(int i = 0; i < 3; i++){
        cv::normalize(planes[i], normalized_planes[i], 0, 255, cv::NORM_MINMAX);
    }
    cv::cv2eigen(normalized_planes[0], b_channel);
    cv::cv2eigen(normalized_planes[1], g_channel);
    cv::cv2eigen(normalized_planes[2], r_channel);
}

int ImageLoader::segmentation(string image_folder){
    /**
    * @function Watershed_and_Distance_Transform.cpp
    * @brief Sample code showing how to segment overlapping objects using Laplacian filtering, in addition to Watershed and Distance Transformation
    * @author OpenCV Team
    */


    //! [load_image]
    // Load the image
    cv::Mat src = cv::imread("../Test Images/" + image_folder + "/Test image.jpg");

    // Check if everything was fine
    if (!src.data)
        return -1;

    // Show source image
//    cv::imshow("Source Image", src);
    //! [load_image]

    //! [black_bg]
    // Change the background from white to black, since that will help later to extract
    // better results during the use of Distance Transform
    for( int x = 0; x < src.rows; x++ ) {
        for( int y = 0; y < src.cols; y++ ) {
            if ( src.at<Vec3b>(x, y) == Vec3b(255,255,255) ) {
                src.at<cv::Vec3b>(x, y)[0] = 0;
                src.at<cv::Vec3b>(x, y)[1] = 0;
                src.at<cv::Vec3b>(x, y)[2] = 0;
            }
        }
    }

    // Show output image
//    cv::imshow("Black Background Image", src);
    //! [black_bg]

    //! [sharp]
    // Create a kernel that we will use for accuting/sharpening our image
    cv::Mat kernel = (Mat_<float>(3,3) <<
            1, 1, 1,
            1, -8, 1,
            1, 1, 1); // an approximation of second derivative, a quite strong kernel

    // do the laplacian filtering as it is
    // well, we need to convert everything in something more deeper then CV_8U
    // because the kernel has some negative values,
    // and we can expect in general to have a Laplacian image with negative values
    // BUT a 8bits unsigned int (the one we are working with) can contain values from 0 to 255
    // so the possible negative number will be truncated
    cv::Mat imgLaplacian;
    cv::Mat sharp = src; // copy source image to another temporary one
    cv::filter2D(sharp, imgLaplacian, CV_32F, kernel);
    src.convertTo(sharp, CV_32F);
    cv::Mat imgResult = sharp - imgLaplacian;

    // convert back to 8bits gray scale
    imgResult.convertTo(imgResult, CV_8UC3);
    imgLaplacian.convertTo(imgLaplacian, CV_8UC3);

    // imshow( "Laplace Filtered Image", imgLaplacian );
    //cv::imshow( "New Sharped Image", imgResult );
    //! [sharp]

    src = imgResult; // copy back

    //! [bin]
    // Create binary image from source image
    cv::Mat bw;
    cv::cvtColor(src, bw, cv::COLOR_BGR2GRAY);
    cv::threshold(bw, bw, 40, 255, cv::THRESH_BINARY | cv::THRESH_OTSU);
//    imshow("Binary Image", bw);
    //! [bin]

    /// dilation:
    cv::Mat dilated;
    cv::Mat kernel2 = cv::getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(3, 3));//, cv::Point(0,0));
    cv::morphologyEx(bw, dilated, MORPH_CLOSE, kernel2);
    //cv::imshow("Closed", dilated);
    //cv::waitKey(0);

    //! [dist]
    // Perform the distance transform algorithm
    cv::Mat dist;
    cv::distanceTransform(bw, dist, cv::DIST_LABEL_PIXEL, 3);

    // Normalize the distance image for range = {0.0, 1.0}
    // so we can visualize and threshold it
    cv::normalize(dist, dist, 0, 1., cv::NORM_MINMAX);
//    cv::imshow("Distance Transform Image", dist);
    //! [dist]

    //! [peaks]
    // Threshold to obtain the peaks
    // This will be the markers for the foreground objects
    cv::threshold(dist, dist, .3, 1., cv::THRESH_BINARY);

    // Dilate a bit the dist image
    cv::Mat kernel1 = cv::Mat::ones(2, 2, CV_8UC1);
    cv::dilate(dist, dist, kernel1);
//    imshow("Peaks", dist);
    //! [peaks]

    //! [seeds]
    // Create the CV_8U version of the distance image
    // It is needed for findContours()
    cv::Mat dist_8u;
    dist.convertTo(dist_8u, CV_8U);

    // Find total markers
    vector<vector<cv::Point> > contours;
    cv::findContours(dist_8u, contours, cv::RETR_EXTERNAL, cv::CHAIN_APPROX_SIMPLE);

    // Create the marker image for the watershed algorithm
    cv::Mat markers = cv::Mat::zeros(dist.size(), CV_32SC1);

    // Draw the foreground markers
    for (size_t i = 0; i < contours.size(); i++)
        cv::drawContours(markers, contours, static_cast<int>(i), cv::Scalar::all(static_cast<int>(i)+1), -1);

    // Draw the background marker
    cv::circle(markers, cv::Point(5,5), 3, CV_RGB(255,255,255), -1);
//    cv::imshow("Markers", markers*10000);
    //! [seeds]

    //! [watershed]
    // Perform the watershed algorithm
    cv::watershed(src, markers);

    cv::Mat mark = cv::Mat::zeros(markers.size(), CV_8UC1);
    markers.convertTo(mark, CV_8UC1);
    cv::bitwise_not(mark, mark);
//        imshow("Markers_v2", mark); // uncomment this if you want to see how the mark
    // image looks like at that point

    // Generate random colors
    vector<cv::Vec3b> colors;
    for (size_t i = 0; i < contours.size(); i++)
    {
        int b = cv::theRNG().uniform(0, 255);
        int g = cv::theRNG().uniform(0, 255);
        int r = cv::theRNG().uniform(0, 255);

        colors.push_back(cv::Vec3b((uchar)b, (uchar)g, (uchar)r));
    }

    // Create the result image
    cv::Mat dst = cv::Mat::zeros(markers.size(), CV_8UC3);

    // Fill labeled objects with random colors
    for (int i = 0; i < markers.rows; i++)
    {
        for (int j = 0; j < markers.cols; j++)
        {
            int index = markers.at<int>(i,j);
            if (index > 0 && index <= static_cast<int>(contours.size()))
                dst.at<cv::Vec3b>(i,j) = colors[index-1];
            else
                dst.at<cv::Vec3b>(i,j) = cv::Vec3b(0,0,0);
        }
    }

    // Visualize the final image
    cv::imshow("Final Result", dst);
    //! [watershed]

    cv::waitKey(0);


    test_image = dst;
    return 0;

}


void ImageLoader::kMeansClustering(int num_clusters, string image_folder){
    cv::Mat src = cv::imread("../Test Images/" + image_folder + "/Test image.jpg");
    Mat samples(src.rows * src.cols, 3, CV_32F);
    for( int y = 0; y < src.rows; y++ )
        for( int x = 0; x < src.cols; x++ )
            for( int z = 0; z < 3; z++)
                samples.at<float>(y + x*src.rows, z) = src.at<Vec3b>(y,x)[z];


    int clusterCount = num_clusters;
    Mat labels;
    int attempts = 20;
    Mat centers;
    kmeans(samples, clusterCount, labels, TermCriteria(CV_TERMCRIT_ITER|CV_TERMCRIT_EPS, 10, 1), attempts, KMEANS_PP_CENTERS, centers );


    Mat new_image( src.size(), src.type() );
    for( int y = 0; y < src.rows; y++ )
        for( int x = 0; x < src.cols; x++ )
        {
            int cluster_idx = labels.at<int>(y + x*src.rows,0);
            new_image.at<Vec3b>(y,x)[0] = centers.at<float>(cluster_idx, 0);
            new_image.at<Vec3b>(y,x)[1] = centers.at<float>(cluster_idx, 1);
            new_image.at<Vec3b>(y,x)[2] = centers.at<float>(cluster_idx, 2);
        }
    test_image = new_image;
    nonmodified_image = src;
    string title = to_string(num_clusters) + " clusters";
    imshow(title, new_image );
    waitKey( 0 );
}

