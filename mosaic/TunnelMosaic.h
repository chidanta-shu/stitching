#pragma once
#include <vector>
#include "mosaic.h"
#include "stdafx.h"
#include <direct.h>
#include <stdlib.h>
#include <stdio.h>
#include <direct.h> //_mkdir������ͷ�ļ�
#include <io.h>     //_access������ͷ�ļ�
#include <iostream> //��ʵ���˼�����İ�����mathͷ
#include <algorithm>
#include <string>
//#include <math.h>
#include <opencv2/opencv.hpp>
#include "TunnelMosaic.h"
#include <ColorCorrection.h>
#include <cmath>
#include <Shlwapi.h>

//#include "Bitmap.h"
#define	MOSAIC_API __declspec(dllexport)

using namespace std;
using namespace cv;

//class ImageMap;
//class cv::Mat;
//class ImageOut;

class MOSAIC_API TunnelMosaic
{
public:
	TunnelMosaic() { ; };
	~TunnelMosaic();

	void setconfigpath(std::string cfgpath) {
		configpath = cfgpath;
	};
	void settemppath(std::string tpath) {
		temppath = tpath;
	};
	void setoutputpath(std::string opath) {
		outputpath = opath;
	};
	//void start(string rawPath, string resultPath);
	//bool loadconfig();
	//void runAxial();
	//void runAxial(TunnelMosaic& preRadia);
	//void runRadial(std::vector<std::string>  imgpathlist, int ipos = 0);
	//void runRadial2(std::vector<std::string>  imgpathlist, int ipos = 0);
	//void runRadial23(std::vector<std::string>  imgpathlist, int ipos = 0); //���� 
	void OptimizeSeam();
	void my_stitching(std::vector<std::string> imgpathlist, Mat* dst);

	void start_stitching(string rawPath, string resultPath);
	void startMosaicAlgorithm(string rawPath, string resultPath);
	int getCurrImageNum();
	
private:

	//��Ÿ�������Ĳ����Լ����ò���
	std::string configpath;
	//��ʱ·��
	std::string temppath;
	//���Ŀ¼
	std::string outputpath;
	//�������
	int sumRadial = 6;
	int sumAxial;
public:	
	int		medpos;
	long	lwidth;
	long	lheight;
	long	lmedleft;
	int num = 0;
	string srcPath;
	string dstPath;
	//ͼ���б�
	//std::vector<ImageMap*>  imagemaps;
	//���
	//ImageOut*    Radialout = NULL;
};

//static void startMosaicAlgorithm(string rawPath, string resultPath);

//Mat convertTo3Channels(const Mat& binImg)
//{
//	Mat three_channel = Mat::zeros(binImg.rows, binImg.cols, CV_8UC3);
//	vector<Mat> channels;
//	for (int i = 0; i < 3; i++)
//	{
//		channels.push_back(binImg);
//	}
//	merge(channels, three_channel);
//	return three_channel;
//}

//void startMosaicAlgorithm(string rawPath, string resultPath)
//{
//	TunnelMosaic   Mosaic;
//
//	Mosaic.setconfigpath(rawPath);
//	Mosaic.settemppath(rawPath);
//	Mosaic.setoutputpath(rawPath);
//	Mosaic.loadconfig();
//
//	Mosaic.start_stitching(rawPath, resultPath);
//
//
//
//
//
//}



//static void startMosaicAlgorithm(string rawPath, string resultPath)
//{
//	startMosaicAlgorithm(rawPath, resultPath);
//
//}

