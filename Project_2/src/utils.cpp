#include <opencv2/core/mat.hpp>
#include "utils.h"

/////////////////////////////////////////////////////////
double rgbDistance(pixel_t x, pixel_t y, const Eigen::MatrixXi &red, const Eigen::MatrixXi &green,
                   const Eigen::MatrixXi &blue)
{
    double dist = sqrt( pow(red(x.row, x.col) - red(y.row, y.col), 2)
                        + pow(green(x.row, x.col) - green(y.row, y.col), 2)
                        + pow(blue(x.row, x.col) - blue(y.row, y.col), 2) );
    return dist;
}

/////////////////////////////////////////////////////////
void setUserArgs(int argc, char **argv, double &mutation_rate, double &crossover_rate,
                 uint16_t &tournament_size, double &time_limit, uint16_t &generation_limit, uint16_t &population_size,
                 int &problem_num)
{
    mutation_rate = 0;
    crossover_rate = 0;
    tournament_size = 0;
    time_limit = 0;
    generation_limit = 0;
    population_size = 0;
    for (uint8_t i = 1; i < argc; i+=2)
    {
        char* arg_name = argv[i];
        double arg_val = atof(argv[i+1]);
        if (!strcmp("mutation_rate", arg_name))
        {
            mutation_rate = arg_val;
            continue;
        }
        else if(!strcmp("crossover_rate", arg_name))
        {
            crossover_rate = arg_val;
            continue;
        }
        else if(!strcmp("tournament_size", arg_name))
        {
            tournament_size = uint16_t(arg_val);
            continue;
        }
        else if(!strcmp("time_limit", arg_name))
        {
            time_limit = arg_val;
            continue;
        }
        else if(!strcmp("generation_limit", arg_name))
        {
            generation_limit = uint16_t(arg_val);
            continue;
        }
        else if(!strcmp("population_size", arg_name))
        {
            population_size = uint16_t(arg_val);
            continue;
        }
        else if(!strcmp("problem_num", arg_name))
        {
            problem_num = int(arg_val);
            continue;
        }
    }
}

/////////////////////////////////////////////////////////
void printMST(std::vector<int> &parent, int num_pixels,
              const Eigen::MatrixXi &red, const Eigen::MatrixXi &green, const Eigen::MatrixXi &blue)
{
    pixel_t x(0,0), y(0,0);
    auto cols = uint16_t(red.cols());
    printf("Edge   Weight\n");
    for (int i = 1; i < num_pixels; i++) {
        x.row = i / cols;
        x.col = i % cols;
        y.row = parent[i] / cols;
        y.col = parent[i] % cols;
        printf("%d - %d    %f \n", parent[i], i, rgbDistance(x, y, red, green, blue));
    }
}


/**
 * Code for thinning a binary image using Zhang-Suen algorithm.
 *
 * Author:  Nash (nash [at] opencv-code [dot] com)
 * Website: http://opencv-code.com
 */

/**
 * Perform one thinning iteration.
 * Normally you wouldn't call this function directly from your code.
 *
 * Parameters:
 * 		im    Binary image with range = [0,1]
 * 		iter  0=even, 1=odd
 */
void thinningIteration(cv::Mat& img, int iter)
{
    CV_Assert(img.channels() == 1);
    CV_Assert(img.depth() != sizeof(uchar));
    CV_Assert(img.rows > 3 && img.cols > 3);

    cv::Mat marker = cv::Mat::zeros(img.size(), CV_8UC1);

    int nRows = img.rows;
    int nCols = img.cols;

    if (img.isContinuous()) {
        nCols *= nRows;
        nRows = 1;
    }

    int x, y;
    uchar *pAbove;
    uchar *pCurr;
    uchar *pBelow;
    uchar *nw, *no, *ne;    // north (pAbove)
    uchar *we, *me, *ea;
    uchar *sw, *so, *se;    // south (pBelow)

    uchar *pDst;

    // initialize row pointers
    pAbove = NULL;
    pCurr  = img.ptr<uchar>(0);
    pBelow = img.ptr<uchar>(1);

    for (y = 1; y < img.rows-1; ++y) {
        // shift the rows up by one
        pAbove = pCurr;
        pCurr  = pBelow;
        pBelow = img.ptr<uchar>(y+1);

        pDst = marker.ptr<uchar>(y);

        // initialize col pointers
        no = &(pAbove[0]);
        ne = &(pAbove[1]);
        me = &(pCurr[0]);
        ea = &(pCurr[1]);
        so = &(pBelow[0]);
        se = &(pBelow[1]);

        for (x = 1; x < img.cols-1; ++x) {
            // shift col pointers left by one (scan left to right)
            nw = no;
            no = ne;
            ne = &(pAbove[x+1]);
            we = me;
            me = ea;
            ea = &(pCurr[x+1]);
            sw = so;
            so = se;
            se = &(pBelow[x+1]);

            int A  = (*no == 0 && *ne == 1) + (*ne == 0 && *ea == 1) +
                     (*ea == 0 && *se == 1) + (*se == 0 && *so == 1) +
                     (*so == 0 && *sw == 1) + (*sw == 0 && *we == 1) +
                     (*we == 0 && *nw == 1) + (*nw == 0 && *no == 1);
            int B  = *no + *ne + *ea + *se + *so + *sw + *we + *nw;
            int m1 = iter == 0 ? (*no * *ea * *so) : (*no * *ea * *we);
            int m2 = iter == 0 ? (*ea * *so * *we) : (*no * *so * *we);

            if (A == 1 && (B >= 2 && B <= 6) && m1 == 0 && m2 == 0)
                pDst[x] = 1;
        }
    }

    img &= ~marker;
}

/**
 * Function for thinning the given binary image
 *
 * Parameters:
 * 		src  The source image, binary with range = [0,255]
 * 		dst  The destination image
 */
void thinning(const cv::Mat& src, cv::Mat& dst)
{
    dst = src.clone();
    dst /= 255;         // convert to binary image

    cv::Mat prev = cv::Mat::zeros(dst.size(), CV_8UC1);
    cv::Mat diff;

    do {
        thinningIteration(dst, 0);
        thinningIteration(dst, 1);
        cv::absdiff(dst, prev, diff);
        dst.copyTo(prev);
    }
    while (cv::countNonZero(diff) > 0);

    dst *= 255;
}


