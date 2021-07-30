#pragma once
#include <iostream>
#include <opencv2/opencv.hpp>
#include <cassert>  
#include <vector>  
#include <cmath>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#define num_img 2 /*6*/  // 更新2
//#define gr = 2

using namespace cv;
using namespace std;

struct GOMY
{
	vector<float> Bi;
	vector<float> Sicb;
	vector<float> Sicr;
	//待定值
	int w = 0;
	int h = 0;
};

struct minhang
{
	float mindiff = 0;
	int seamindex = 0;

};


GOMY get_overlap_mean_yuv(Mat img1, Mat img2, int shift_x, int shift_y);

vector<float> getCorrParameters(vector<vector<float> > B);

minhang minhh(vector<vector<float> > a, int m);

Mat RGB2YUV0(Mat src);
Mat YUV2RGB0(Mat src);

void Color_correction(Mat *im1, Mat *im2);


