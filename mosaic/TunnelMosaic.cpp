#include "stdafx.h"
#include <direct.h>
#include <stdlib.h>
#include <stdio.h>
#include <direct.h> //_mkdir函数的头文件
#include <io.h>     //_access函数的头文件
#include <iostream> //其实拐了几个弯的包含了math头
#include <algorithm>
#include <string>
//#include <math.h>
#include <opencv2/opencv.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/features2d.hpp>
#include <opencv2/flann.hpp>
//#include <opencv2/xfeatures2d.hpp>
#include "TunnelMosaic.h"
#include "ColorCorrection.h"
#include <cmath>

//#include "Bitmap.h"

//#ifdef LOG_OUTPUT
//#include <glog/logging.h>
//#endif
using namespace std;
using namespace cv;
//using namespace cv::xfeatures2d;

typedef struct
{
	Point2f left_top;
	Point2f left_bottom;
	Point2f right_top;
	Point2f right_bottom;
}four_corners_t;

static void CreateDir(const char *dir)
{
	int m = 0, n;
	std::string str1, str2;

	str1 = dir;
	str2 = str1.substr(0, 2);
	str1 = str1.substr(3, str1.size());

	while (m >= 0)
	{
		m = str1.find('\\');

		str2 += '\\' + str1.substr(0, m);
		//判断该目录是否存在
		n = _access(str2.c_str(), 0);
		if (n == -1)
		{	//创建目录
			_mkdir(str2.c_str());
		}
		str1 = str1.substr(m + 1, str1.size());
	}
}

static void CalcCorners(const cv::Mat& H, Size size, four_corners_t &corners)
{
	double v2[] = { 0, 0, 1 };//左上角
	double v1[3];//变换后的坐标值
	Mat V2 = Mat(3, 1, CV_64FC1, v2);  //列向量
	Mat V1 = Mat(3, 1, CV_64FC1, v1);  //列向量
	V1 = H * V2;

	//左上角(0,0,1)
	corners.left_top.x = v1[0] / v1[2];
	corners.left_top.y = v1[1] / v1[2];

	//左下角(0,src.rows,1)
	v2[0] = 0;
	v2[1] = size.height;//src.rows;
	v2[2] = 1;
	V2 = Mat(3, 1, CV_64FC1, v2);  //列向量
	V1 = Mat(3, 1, CV_64FC1, v1);  //列向量
	V1 = H * V2;
	corners.left_bottom.x = v1[0] / v1[2];
	corners.left_bottom.y = v1[1] / v1[2];

	//右上角(src.cols,0,1)
	v2[0] = size.width;// src.cols;
	v2[1] = 0;
	v2[2] = 1;
	V2 = Mat(3, 1, CV_64FC1, v2);  //列向量
	V1 = Mat(3, 1, CV_64FC1, v1);  //列向量
	V1 = H * V2;
	corners.right_top.x = v1[0] / v1[2];
	corners.right_top.y = v1[1] / v1[2];

	//右下角(src.cols,src.rows,1)
	v2[0] = size.width;// src.cols;
	v2[1] = size.height;//src.rows;
	v2[2] = 1;
	V2 = Mat(3, 1, CV_64FC1, v2);  //列向量
	V1 = Mat(3, 1, CV_64FC1, v1);  //列向量
	V1 = H * V2;
	corners.right_bottom.x = v1[0] / v1[2];
	corners.right_bottom.y = v1[1] / v1[2];
}

static void cylindrical_projection(cv::Mat *src) 
{
	//Mat srcImage = imread("E:\\Temp\\files\\1\\1.jpg");
	//Mat srcImage1 = imread("E:\\Temp\\files\\1\\2.jpg");

	Mat srcImage = *src;

	int height = srcImage.rows; //原图像的高（即原图像矩阵行数）
	int width = srcImage.cols; //原图像的宽（即原图像矩阵列数）
	int centerX = width / 2; //图像中心横坐标
	int centerY = height / 2; //图像中心纵坐标
	double alpha = CV_PI / 4; //相机视角角度
	double f = width / (2 * tan(alpha / 2)); //焦距（圆的半径）8

	//求左右黑色缝隙宽度
	int len = (width / 2 - f * alpha / 2); //cvRound:取整

	//定义合适的目标图
	Mat dstImage = Mat::zeros(srcImage.rows, width - 2 * len, CV_8UC3); //注意尺寸

	//循环遍历
	for (int i = 0; i < srcImage.rows; i++)
	{
		for (int j = 0; j < srcImage.cols; j++)
		{
			//注意图像坐标与像素矩阵坐标的区别
			float theta = atan((j - centerX) / f);
			int pointX = cvRound(width / 2 + f * theta); //注意这里用width / 2，用f * alpha / 2,两边缝隙会不均匀（只有右边有黑缝隙）
			int pointY = cvRound(f * (i - centerY) / sqrt((j - centerX) * (j - centerX) + f * f) + centerY);

			//像素赋值,此时要将列数（图像横坐标）往左移，与初始的黑色图像边缘对其，即pointX - len
			dstImage.at<Vec3b>(pointY, pointX - len)[0] = srcImage.at<Vec3b>(i, j)[0];
			dstImage.at<Vec3b>(pointY, pointX - len)[1] = srcImage.at<Vec3b>(i, j)[1];
			dstImage.at<Vec3b>(pointY, pointX - len)[2] = srcImage.at<Vec3b>(i, j)[2];
		}
	}




	*src = dstImage.clone();
}

static void my_stitching_method(Mat* im1, Mat* im2, Mat* dst, int shiftx, int shifty, int _gap)
{
	
		Mat left = *im1;
		Mat right = *im2;
		int shiftx12 = shiftx;
		int shifty12 = shifty + 50;
		int gap = _gap;
		int r = left.rows;

		if (left.rows == right.rows)
			r = r * 1.05;

		Mat temp(r, left.cols + right.cols + shiftx12 - gap, left.type(), Scalar(0));

		//imshow("1.jpg", dstImage);
		//waitKey();
		if (left.rows == right.rows) {
			Mat half(temp, cv::Rect(0, 50, left.cols, left.rows));
			Mat halfm(temp, cv::Rect(left.cols + shiftx12, shifty12, right.cols - gap, right.rows));//
			Mat halfo(right, cv::Rect(gap, 0, right.cols - gap, right.rows));
			left.copyTo(half);
			halfo.copyTo(halfm);
		}
		else {
			Mat half(temp, cv::Rect(0, 0, left.cols, left.rows));
			Mat halfm(temp, cv::Rect(left.cols + shiftx12, shifty12, right.cols - gap, right.rows));//
			Mat halfo(right, cv::Rect(gap, 0, right.cols - gap, right.rows));
			left.copyTo(half);
			halfo.copyTo(halfm);

		}
		//imwrite("dst.jpg", temp);
		*dst = temp.clone();
	
}

static bool cmp(string a[], string b[]) {
	// 如果身高相同
	if (a[0] == b[0])
		// 按照名字升序排序
		return a[1] <= b[1];
	else
		// 按照身高升序排序
		return a[0] <= b[0];
}

static bool comparator(const string& src1, const string& src2) {

	//string find_src = "_";//

	size_t num1 = stoi(src1.substr(0, src1.size()));//关于字符串的使用可以查看末尾的链接
	size_t num2 = stoi(src2.substr(0, src2.size()));//stoi是将字符转数字，
	if (num1 > num2) {
		return 0;
	}
	else
	{
		return 1;//return 1 实现从小到大 
	}
}


//static bool computePairNum(std::pair<double, std::string> pair1, std::pair<double, std::string> pair2)
//{
//	return pair1.first < pair2.first;
//}

static void sort_filelists(std::vector<std::string>& filists, std::vector<std::string>& newfilists)
{
	if (filists.empty())return;
	//std::vector<std::pair<double, std::string> > filelists_pair;
	//for (int i = 0; i < filists.size(); ++i) {
	//	std::string tmp_string = filists[i];
	//	int npos = tmp_string.find_first_of("_");
	//	std::string tmp_num_string = tmp_string.substr(tmp_string.at(0), npos );
	////	double tmp_num = atof(tmp_num_string.c_str());
	////	std::pair<double, std::string> tmp_pair;
	////	tmp_pair.first = tmp_num;
	////	tmp_pair.second = tmp_string;
	////	filelists_pair.push_back(tmp_pair);
	////}
	////std::sort(filelists_pair.begin(), filelists_pair.end(), computePairNum);
	//////filists.clear();
	//for (int i = 0; i < filelists_pair.size(); ++i) {
	//	newfilists.push_back(filelists_pair[i].second);
	//}
	//string ** info = new string *[filists.size()];

	vector<size_t> num;

	for (size_t i = 0; i < filists.size(); i++)
	{
		int npos = filists[i].find_first_of("_");
		num.push_back(stoi(filists[i].substr(0, npos)));
		
	}

	for (size_t i = 0; i < filists.size(); i++)
	{
		for (size_t j = 0; j < filists.size() - i - 1; j++)
		{
			if (num[j] > num[j + 1])//如果这个数比前面的大
			{
				swap(num[j], num[j + 1]);
				swap(filists[j], filists[j + 1]);
			}
		}
	}

	//for (int i = 0; i < filists.size(); ++i)
	//{
	//	info[i] = new string[2];
	//	info[i][0] = num[i];
	//	info[i][1] = filists[i];
	//	//cout << info[i][0] << endl;
	//	//cout << info[i][1] << endl;
	//}
	//// 排序
	////StrCmpLogicalW(,);
	//int num1 = stoi(num[0]);
	//cout << num1 << endl;
	//std::sort(info, info + filists.size(), cmp);//
	//for (int i = 0; i < filists.size(); ++i)
	//{
	//	cout << info[i][0] << endl;
	//	cout << info[i][1] << endl;
	//}
	//// 输出
	////filists.clear();
	
	//for (int i = 0; i < filists.size(); i++)//数组输出
	//{
	//	cout << filists[i] << "\t";
	//}
	newfilists.clear();

	newfilists.insert(newfilists.end(), filists.begin(), filists.end());
	
	//for (size_t i = 0; i < filists.size(); i++) {
	//	newfilists.push_back(info[i][1]);
	//	cout << newfilists[i] << endl;
	//}
	//cout << "sort:" << filists.size() << endl;
}


static void getFiles(string path, string exd, vector<string>& files)
{
	vector<string> tempfileslist;
	//文件句柄
	size_t hFile = 0;
	//文件信息
	struct __finddata64_t fileinfo;
	string pathName, exdName;

	if (0 != strcmp(exd.c_str(), "")) exdName = "\\*" + exd;
	else exdName = "\\*";

	if ((hFile = _findfirst64(pathName.assign(path).append(exdName).c_str(), &fileinfo)) != -1)
	{
		do
		{
			if (strcmp(fileinfo.name, ".") != 0 && strcmp(fileinfo.name, "..") != 0)
				tempfileslist.push_back(fileinfo.name);       // tempfileslist

		} while (_findnext64(hFile, &fileinfo) == 0);
		_findclose(hFile);
	}


	sort_filelists(tempfileslist, files);

	//files.clear();
	//files.insert(files.end(), zhongjian.begin(),zhongjian.end());

	//cout << files[0] << endl;

}

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
//static void lightchange()

//using namespace  cv::detail;

//class ImageOut
//{
//public:
//	ImageOut() { ; };
//	~ImageOut() { ; };
//
//public:
//	cv::Mat outimage;
//
//};
//
//class ImageMap
//{
//public:
//
//	ImageMap(std::string cfgpath) { 
//		configpath = cfgpath; 
//		//加载相机的配置
//		//读的形式打开yml。当然也可以打开xml，主要看后缀
//		FileStorage readfs(configpath, FileStorage::READ); 
//
//		if (readfs.isOpened())
//		{
//
//			readfs["camera_matrix"] >> cameraMatrix;
//			readfs["distortion_coefficients"] >> distCoeffs;
//			readfs["image_width"] >> image_width;
//			readfs["image_height"] >> image_height;
//			//std::cout<<"Camera:"<< cfgpath <<std::endl;
//			if (!cameraMatrix.empty() && !distCoeffs.empty()) 
//			{
//				//Formatter::FMT_PYTHON
//				//Formatter::FMT_NUMPY
//				//Formatter::FMT_CSV
//				//Formatter::FMT_C
//				//std::cout << "cameraMatrix:" << std::endl << format(cameraMatrix, Formatter::FMT_C) << std::endl << std::endl;
//				//std::cout << "distCoeffs:" << std::endl << format(distCoeffs, Formatter::FMT_C) << std::endl << std::endl;
//				bcfgValid = true;
//			}
//			
//		}
//		readfs.release();	
//
//	};
//	~ImageMap() { ; };
//	//是否数据合法
//	bool isvalid() { return bcfgValid; };
//	void setHoriangle(double angle) { Hangle = angle; };
//	Mat deDistort(Mat frame)
//	{
//		if (!bcfgValid)  return frame;
//
//		if (bfirst)
//		{
//			Size imageSize;
//			imageSize = frame.size();
//			initUndistortRectifyMap(cameraMatrix, distCoeffs, Mat(),
//				getOptimalNewCameraMatrix(cameraMatrix, distCoeffs, imageSize, 1, imageSize, 0),
//				imageSize, CV_16SC2, map1, map2);
//		}
//
//		Mat frameCalibration;
//		remap(frame, frameCalibration, map1, map2, INTER_LINEAR);
//		//return tocylinder(frameCalibration);
//		return frameCalibration;
//	}
//
//	Mat tocylinder(Mat src)
//	{
//
//		//是否根据隧道内的直线进行角度调整？
//		Mat dst = src.clone();
//		if (bangleAdjust)
//		{		
//			Mat cannyimg, grayImg;
//			//颜色转为黑白
//			cv::cvtColor(src, grayImg, cv::COLOR_RGB2GRAY);
//			//去噪
//			cv::GaussianBlur(grayImg, grayImg, cv::Size(5, 5), 0, 0);
//			cv::medianBlur(grayImg, grayImg, 5);
//			//直线提取
//			cv::Canny(grayImg, cannyimg, 180, 30, 3);
//			// 霍夫直线
//			std::vector<cv::Vec4i> lines(10000);
//			cv::HoughLinesP(cannyimg, lines, 2, CV_PI / 180, 80, 100, 20);
//
//			double anglesum = 0;
//			int nsize = 0;
//			double maxdist = 0;
//			double maxangle = 0;
//			cv::Vec4i maxline;
//
//			int colorx = 10;
//			//最长直线的角度
//			for (auto val : lines)
//			{
//				//line(showimg2, Point(val[0], val[1]), Point(val[2], val[3]), Scalar(colorx++, 0, 0), 3);
//				double anglex = atan2f(val[3] - val[1], val[2] - val[0]);
//				//double anglex = val[1];
//				//anglesum += anglex * 180.0 / kkk::pai;
//				nsize++;
//				double dist = (val[0] - val[2])*(val[0] - val[2]) + (val[3] - val[1])*(val[3] - val[1]);
//				if (maxdist < dist)
//				{
//					maxdist = dist;
//					maxangle = anglex * 180.0 / PI;
//					maxline = val;
//				}
//
//			}
//
//			//根据角度旋转图像
//			
//			if (nsize)
//			{
//				std::cout << "angle=" << maxangle << std::endl;
//				//std::cout << "dist=" << maxdist << std::endl;
//				line(src, Point(maxline[0], maxline[1]), Point(maxline[2], maxline[3]), Scalar(255, 0, 255), 3);
//
//				if (maxangle > 90.0) maxangle -= 180.0f;
//				else if (maxangle < -90.0) maxangle += 180.0f;
//
//				bool bhori = true;
//				if (maxangle > 45.0f)
//				{
//					maxangle = 90.0f - maxangle;
//					maxangle *= -1.0f;
//					bhori = false;
//				}
//				else if (maxangle < -45.0f)
//				{
//					maxangle = -90.0f - maxangle;
//					maxangle *= -1.0f;
//					bhori = false;
//				}
//				//std::cout << "angle=" << maxangle << std::endl;
//				Point2f center(src.cols / 2.0f, src.rows / 2.0f);//旋转中心  
//				Mat M = getRotationMatrix2D(center, maxangle, 1.0);//计算旋转的仿射变换矩阵 
//				cv::warpAffine(src, dst, M, Size(src.cols, src.rows));//仿射变换
//			}
//
//		}
//		// 相机的垂直张角
//		Vangle = src.rows* Hangle / (double)src.cols;
//		// 等效焦距,内投
//		myf = src.cols*0.5 / tan(TO_RAD(Hangle*0.5));
//		//外投
//		//myf = src.cols*0.5 / sin(TO_RAD(Hangle*0.5));
//
//		// 柱面图像宽度
//		//retmatWidth = (int)(TO_RAD(Hangle) * myf + 0.5);
//		//retmatHeight = (int)(tan(TO_RAD(Vangle*0.5)) * myf * 2 + 0.5);//src.rows
//
//
//		retmatWidth = (int)(TO_RAD(Hangle) * myf + 0.5);// src.cols;//src.cols;
//		retmatHeight = src.rows;//src.rows
//
//		// 新建柱面图片
//		Mat retMat(retmatHeight, retmatWidth, CV_8UC3);
//
//	
//		//std::cout << "Vangle=" << Vangle << std::endl;
//		//std::cout << "myf=" << myf << std::endl;
//		
//#if 0	// 柱面投影 ，有平面投影柱面	
//		for (int i = 0; i < src.rows; i++)
//		{
//			for (int j = 0; j < src.cols; j++)
//			{
//
//				///////////////正运算 ，可能会有空洞
//				float theta = atan((j - src.cols*0.5f) / myf);
//
//				//注意这里用width / 2，用f * alpha两边缝隙不均匀（只有右边有黑缝隙）
//				int xpos = retMat.cols * 0.5f + myf * theta; 
//				int ypos = myf * (i - src.rows*0.5f) / sqrt((j - src.cols*0.5f) * (j - src.cols*0.5f) + 
//								myf * myf) + retMat.rows*0.5f;
//				if (xpos >= 0 && xpos < retMat.cols && ypos >= 0 && ypos < retMat.rows)
//				{
//					retMat.at<Vec3b>(ypos, xpos) = dst.at<Vec3b>(i, j);
//
//				}
//
//			}
//		}
//#else ////逆运算 ,不会有空洞
//		double xdlt = 1.0;
//		
//	
//		for (int i = 0; i < retMat.rows; i++)
//		{
//			for (int j = 0; j < retMat.cols; j++)
//			{
//
//				int retxj = j - retMat.cols * 0.5;
//				int retyi = i - retMat.rows * 0.5;
//
//				//内投
//				//double xj = tan(retxj / myf) * myf + src.cols * 0.5;
//				//  yi = (retyi)* sqrt((retxj * xdlt) * (retxj * xdlt) + myf * myf) / myf + src.rows * 0.5;
//				//double yi = tan(retyi/myf) * myf / cos(retxj / myf) +  src.rows * 0.5;
//
//				//外投
//				//double xg = myf * cos(TO_RAD(Hangle*0.5));// 为图片在圆柱面弦高
//				//double xj = tan(retxj / myf) *  xg + src.cols * 0.5;
//				//double yi = tan(retyi / myf) *  xg / cos( retxj / myf ) + src.rows * 0.5;
//
//				//简单算法
//				float k = myf / sqrt(myf * myf + retxj * retxj);
//				float xj = retxj / k + src.cols * 0.5;;
//				float yi = retyi / k + src.rows * 0.5;;
//
//
//				if (xj >= 0 && xj < src.cols && yi >= 0 && yi < src.rows)
//				{
//					retMat.at<Vec3b>(i, j) = dst.at<Vec3b>(yi, xj);
//
//				}
//
//
//			}
//		}
//#endif
//		return retMat;
//	}
//
//	void getFeaturePoints( int ntype,Mat imagesrc) 
//	{
//		int  nfeatures = 10000;	//10000
//		
//		Ptr<Feature2D> orb;		
//		if (ntype == 0) 
//		{
//			orb = AKAZE::create();
//		}
//		else if (ntype == 1) 
//		{
//			orb = SURF::create();
//		}
//		else
//		{
//			orb = ORB::create(nfeatures);
//		}
//
//		//orb = ORB::create(nfeatures);
//		//Ptr<Feature2D> orb = AKAZE::create();
//		//Ptr<Feature2D> orb = SURF::create();
//		//Ptr<Feature2D> orb = SIFT::create();	
//
//		orb->detectAndCompute(imagesrc, noArray(), kpRef, descRef);
//	}
//
//private:
//	Mat extcameraMatrix;//外参
//	Mat cameraMatrix;//内参
//	Mat distCoeffs; //畸变校正
//	Mat map1;
//	Mat map2;
//	bool bangleAdjust = false;  //是否角度调整，通过已知目标方向，调整图像方向
//	bool bfirst = true;			//是否第一次运行，加载变形数据
//	bool bcfgValid = false;		//配置数据是否有效
//	std::string configpath;
//	// 图像尺寸
//	int image_width;
//	int image_height;
//	//水平张角
//	double Hangle = 45.0;
//	double Vangle = 30.0;
//	//等效焦距
//	double myf ;
//	// 柱面图像宽度高度
//	int retmatWidth;
//	int retmatHeight;
//public:
//	//特征点
//	std::vector<KeyPoint> kpRef;
//	Mat descRef;
//	Mat image;
//
//};



//bool TunnelMosaic::loadconfig()
//{
//
//	for (auto x : imagemaps)
//	{
//		delete x;
//		x = NULL;
//	}
//	imagemaps.clear();
//	
//	for (int i = 0; i < sumRadial; i++)
//	{
//		std::ostringstream ostrm;
//		ostrm << i;
//		std::string filename = outputpath + "\\" + ostrm.str();
//		//输入目录
//		//CreateDir(filename.c_str());
//		std::string cmrcfg = configpath + "outcamera"+ ostrm.str()+".xml";		
//		// 加载相机的参数
//		ImageMap *imgmap = new ImageMap(cmrcfg);
//		// 相机水平张角设置
//		imgmap->setHoriangle(45.0);
//		imagemaps.push_back(imgmap);
//	}
//
//	char namebuf[256];
//	sprintf(namebuf, "%sout\\", outputpath.c_str());
//	//CreateDir(namebuf);
//	//可能意外中断，
//	//加载上次处理的位置
//	return false;
//}

//void TunnelMosaic::runAxial()
//{
//	//中文路径	
//	//setlocale(LC_ALL, "chs");
//
//}



TunnelMosaic::~TunnelMosaic()
{
	//if (Radialout != NULL) delete Radialout;
	//for (std::vector<ImageMap*>::const_iterator iter = imagemaps.begin(); 
	//				iter != imagemaps.end(); ++iter)
	//{
	//	delete (*iter);
	//}
	//imagemaps.clear();
}


//void TunnelMosaic::runRadial(std::vector<std::string>  imgpathlist,int ipos)
//{
//
//	std::vector<Mat> imgsets;
//	char namebuf[256];
//	bool withRotation = false;
//	bool withScale = false;
//
//	Ptr<DescriptorMatcher> matcher = DescriptorMatcher::create("BruteForce-Hamming");
//	int medpos = sumRadial / 2;
//	std::vector<Mat>  homolist;
//
//	for (int i = 0; i < sumRadial; i++)
//	{
//		Mat frame = cv::imread(imgpathlist[i]);
//		if (frame.empty()) continue;
//		
//		//1 标定矫正，是否需要?
//		//Mat dedistmat = imagemaps[i]->deDistort(frame);
//		//if (dedistmat.empty()) continue;
//		
//		//2 柱面投影? 内投和外投
//		//Mat outmat = imagemaps[i]->tocylinder(dedistmat);
//		//if (outmat.empty()) continue;	
//		//imgsets.push_back(outmat);
//		//imgsets.push_back(frame);
//		
//		//3 调整整体亮度？
//
//		//4 提取特征点
//		imagemaps[i]->getFeaturePoints(-1, frame);
//
//		//kpRef, descRef
//		//保存特征点到目标文件，，用于轴向的拼接以及周向的拼接
//		sprintf(namebuf, "%s%d\\outimgfeat%.2d.yaml", outputpath.c_str(),i,ipos);
//		FileStorage fs(namebuf, FileStorage::WRITE);		
//		fs << "feature_points" << imagemaps[i]->kpRef;
//		fs << "points_description" << imagemaps[i]->descRef;
//		fs.release();
//		//特征点输出结束
//		//变形后的图像加入列表，用于轴向的拼接
//		imagemaps[i]->image = frame.clone();
//
//		Mat outmat;
//		//outmat = frame.clone();
//		//drawKeypoints(frame, imagemaps[i]->kpRef, outmat);		
//		//投影后图的输入目录
//		sprintf(namebuf, "%s%d\\outimg%.2d.jpg", outputpath.c_str(),i,ipos);
//		//cv::imwrite(namebuf, frame);		
//		if (i > 0)
//		{
//			//5 环向的拼接进行匹配
//			std::vector<DMatch> matchesAll, matchesGMS;
//			matcher->match(imagemaps[i]->descRef, imagemaps[i - 1]->descRef, matchesAll);			
//			//精确匹配
//			matchGMS(imagemaps[i]->image.size(), imagemaps[i-1]->image.size(), 
//				     imagemaps[i]->kpRef, imagemaps[i-1]->kpRef, 
//				     matchesAll, matchesGMS, withRotation, withScale);
//			//绘制匹配  用于调试
//			//Mat dst; 
//			//drawMatches(imagemaps[i  ]->image, imagemaps[i    ]->kpRef, 
//			//	        imagemaps[i-1]->image, imagemaps[i - 1]->kpRef, matchesGMS, dst);
//			//
//			//6 计算单应矩阵
//			// 将精确匹配点取出，imagePoints1为第一幅图片，imagePoints2为第二幅图像
//			std::vector<Point2f> imagePoints1, imagePoints2;
//			for (int j = 0; j < matchesGMS.size(); j++)
//			{
//				imagePoints2.push_back(imagemaps[i  ]->kpRef[matchesGMS[j].queryIdx].pt);//i
//				imagePoints1.push_back(imagemaps[i-1]->kpRef[matchesGMS[j].trainIdx].pt);//i-1
//			}			
//			Mat homo;
//			//four_corners_t corners;
//			//if (i <= medpos ) 
//			//{
//			//	homo = findHomography(imagePoints1, imagePoints2, cv::RANSAC);
//
//			//	//invert(h12, h21, DECOMP_LU);
//			//	
//			//	CalcCorners(homo, imagemaps[i-1]->image.size(), corners);
//			//}
//			//else
//			//{
//			//	homo = findHomography(imagePoints2, imagePoints1, cv::RANSAC);
//			//}
//			
//#if 0
//			// 计算 1图到2图的单应矩阵
//			homo = findHomography(imagePoints1, imagePoints2, cv::RANSAC);
//			//图像矫正 输出为调试
//			//warpPerspective(imagemaps[i-1]->image, result, homo, Size(2 * imagemaps[i]->image.cols, imagemaps[i]->image.rows));
//#else
//			// 计算 后图到前图的单应矩阵 ， 
//			homo = findHomography(imagePoints2, imagePoints1, cv::RANSAC);			
//			//图像矫正  输出为调试
//			//warpPerspective(imagemaps[i]->image, result, homo, Size(2 * imagemaps[i]->image.cols, imagemaps[i]->image.rows));
//#endif
//
//			homolist.push_back(homo);
//		}
//		
//	}
//
//	Mat result;
//	four_corners_t corners;
//	Mat tempLH = Mat::eye(3, 3, CV_64FC1);
//	//全白图
//	Mat mymask(imagemaps[0]->image.size(), CV_8U, Scalar(255));
//	Mat xtmp;
//	Mat xtmpmask;
//	xtmp = imagemaps[2]->image.clone();
//	for (int i = 0; i < homolist.size(); i++)
//	{
//		tempLH = tempLH * homolist[i];
//		//tempLH = tempLH / tempLH.at<double>(2, 2);
//		//std::cout << "H:" << i << ":" << std::endl;
//		//std::cout << tempLH << std::endl;
//
//		four_corners_t corners;
//		CalcCorners(tempLH, imagemaps[i+1]->image.size(), corners);
//		//图像矫正
//		Mat result;
//		Mat resultmask;
//		warpPerspective(imagemaps[i+1]->image, result, tempLH, Size(4 * imagemaps[i]->image.cols, imagemaps[i]->image.rows));
//		warpPerspective(mymask, resultmask, tempLH, Size(4 * mymask.cols, mymask.rows));
//
//		if (i == 0)
//		{
//			Mat half(result, cv::Rect(0, 0, imagemaps[i]->image.cols, imagemaps[i]->image.rows));
//			Mat halfm(resultmask, cv::Rect(0, 0, imagemaps[i]->image.cols, imagemaps[i]->image.rows));
//			imagemaps[i]->image.copyTo(half);
//			mymask.copyTo(halfm);
//			xtmpmask = resultmask.clone();
//		}
//		else
//		{	
//
//			//char namebuf[256];
//			//sprintf(namebuf, "%soutimg%.2dori_new.jpg", outputpath.c_str(), i);
//			//cv::imwrite(namebuf, resultmask);
//			//sprintf(namebuf, "%soutimg%.2dori.jpg", outputpath.c_str(), i);			
//			//cv::imwrite(namebuf, xtmpmask);
//			Mat element2 = getStructuringElement(MORPH_RECT, Size(5, 5));
//			cv::erode(xtmpmask, xtmpmask, element2);
//			xtmp.copyTo(result, xtmpmask);
//			/*Mat tempresult;
//			Point center(xtmpmask.cols / 2, xtmpmask.rows / 2);
//			seamlessClone(xtmp, result, xtmpmask, center, tempresult, MONOCHROME_TRANSFER);
//			result = tempresult.clone();*/
//			resultmask.copyTo(xtmpmask, resultmask);		
//		}
//	 
//		//变行后图的输入目录
//		//char namebuf[256];
//		sprintf(namebuf, "%soutimg%.2d.jpg", outputpath.c_str(),i);
//		cv::imwrite(namebuf, result);
//		//sprintf(namebuf, "%soutimg%.2dmask.jpg", outputpath.c_str(), i);
//		//cv::imwrite(namebuf, resultmask);
//		xtmp = result.clone();
//	}
//}
//
//   //606-996
//void TunnelMosaic::runRadial2(std::vector<std::string>  imgpathlist, int ipos)
//{
//
//	std::vector<Mat> imgsets;
//	char namebuf[256];
//	bool withRotation = false;  // true
//	bool withScale = false;   // true
//
//	Ptr<DescriptorMatcher> matcher = DescriptorMatcher::create("FlannBased");  //BruteForce-Hamming   FlannBased
//	// orb描述子深度是CV_8U  可以用汉明距离的
//	// sift，surf，为CV32F  float类型的匹配方式有：FlannBasedMatcher，BruteForce<L2>，BruteForce<SL2>，BruteForce<L1>。uchar类型的匹配方式有：BruteForce，BruteForce
//	// ORB和BRIEF特征描述子只能使用BruteForce匹配法。
//	medpos = sumRadial / 2;
//
//	std::vector<Mat>  homolist;
//
//
//
//	for (int i = 0; i < sumRadial; i++)
//	{
//		//Mat frame = cv::imread(imgpathlist[i]);
//		Mat frame = cv::imread(imgpathlist[i], IMREAD_GRAYSCALE);
//		//cv::resize(frame, frame, cv::Size(1000, 1000));
//		//blur(frame, frame, Size(5, 5), Point(-1, -1));
//		//GaussianBlur(frame,frame,Size(5,5),11,11);
//		//cout << frame.type() << endl;
//		////cv::imshow("temp", frame);
//		//Mat tem;
//		//frame.convertTo(tem, CV_8UC1);
//		//cout << tem.type() << endl;
//		//cvtColor(frame, frame, COLOR_GRAY2RGB);
//		if (frame.empty()) continue;
//
//		//1 标定矫正，是否需要?
//		//Mat dedistmat = imagemaps[i]->deDistort(frame);
//		//if (dedistmat.empty()) continue;
//
//		//2 柱面投影? 内投和外投
//		//Mat outmat = imagemaps[i]->tocylinder(dedistmat);
//		//if (outmat.empty()) continue;	
//		//imgsets.push_back(outmat);
//		//imgsets.push_back(frame);
//
//		//3 调整整体亮度？
//
//		//////////////////////
//		//颜色校正？
//
//
//
//
//
//		//4 提取特征点
//		imagemaps[i]->getFeaturePoints(1, frame);
//
//		//kpRef, descRef
//		//保存特征点到目标文件，，用于轴向的拼接以及周向的拼接
//		sprintf(namebuf, "%s%d\\outimgfeat%.2d.yaml", outputpath.c_str(), i, ipos);
//		FileStorage fs(namebuf, FileStorage::WRITE);
//		fs << "feature_points" << imagemaps[i]->kpRef;
//		fs << "points_description" << imagemaps[i]->descRef;
//		fs.release();
//		//特征点输出结束
//		//变形后的图像加入列表，用于轴向的拼接
//		imagemaps[i]->image = frame.clone();
//		Mat outmat;
//		//outmat = frame.clone();
//		//drawKeypoints(frame, imagemaps[i]->kpRef, outmat);		
//		//投影后图的输入目录
//		//sprintf(namebuf, "%s%d\\outimg%.2d.jpg", outputpath.c_str(), i, ipos);
//		//cv::imwrite(namebuf, frame);
//
//		if (i > 0)
//		{
//			//5 环向的拼接进行匹配
//			std::vector<DMatch> matchesAll, matchesGMS;
//			matcher->match(imagemaps[i]->descRef, imagemaps[i - 1]->descRef, matchesAll);
//
//			//test
//			//cv::imwrite("E:\\Temp\\out\\2.png", imagemaps[i]->image);
//			//cv::imwrite("E:\\Temp\\out\\1.png", imagemaps[i - 1]->image);
//			Mat tt;
//			Mat tt1;
//			cv::drawKeypoints(imagemaps[i - 1]->image, imagemaps[i - 1]->kpRef, tt);
//			cv::drawKeypoints(imagemaps[i]->image, imagemaps[i]->kpRef, tt1);
//			cv::imwrite("E:\\Temp\\out\\test.png", tt);
//			cv::imwrite("E:\\Temp\\out\\test1.png", tt1);
//
//			//精确匹配
//			matchGMS(imagemaps[i]->image.size(), imagemaps[i - 1]->image.size(),
//				imagemaps[i]->kpRef, imagemaps[i - 1]->kpRef,
//				matchesAll, matchesGMS, withRotation, withScale);
//			Mat finalMatches;
//			cv::drawMatches(imagemaps[i]->image, imagemaps[i]->kpRef, imagemaps[i - 1]->image, imagemaps[i - 1]->kpRef, 
//				matchesGMS, finalMatches, Scalar::all(-1), Scalar::all(-1), std::vector<char>(), DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS);
//			cv::imwrite("E:\\Temp\\out\\MatchesGMS.png", finalMatches);
//
//			//绘制匹配  用于调试
//			//Mat dst; 
//			//drawMatches(imagemaps[i  ]->image, imagemaps[i    ]->kpRef, 
//			//	        imagemaps[i-1]->image, imagemaps[i - 1]->kpRef, matchesGMS, dst);
//			//
//
//			//6 计算单应矩阵
//			// 将精确匹配点取出，imagePoints1为第一幅图片，imagePoints2为第二幅图像
//			std::vector<Point2f> imagePoints1, imagePoints2;
//			for (int j = 0; j < matchesGMS.size(); j++)
//			{
//				imagePoints2.push_back(imagemaps[i]->kpRef[matchesGMS[j].queryIdx].pt);//i
//				imagePoints1.push_back(imagemaps[i - 1]->kpRef[matchesGMS[j].trainIdx].pt);//i-1
//			}
//
//
//			Mat homo;
//	
//			if (i <= medpos ) 
//			{
//				homo = findHomography(imagePoints1, imagePoints2, cv::RANSAC);
//		
//			}
//			else
//			{
//				homo = findHomography(imagePoints2, imagePoints1, cv::RANSAC);
//			}			
//			//invert(h12, h21, DECOMP_LU);
//			homolist.push_back(homo);
//		}
//	}
//
//	Mat result;
//	four_corners_t corners;
//
//	//
//	vector<int> maxsizelist(sumRadial,0);
//
//	//计算拼接后图像的尺寸，以medpos图像为基准	
//	Mat tempLH = Mat::eye(3, 3, CV_64FC1);
//	for (int i = medpos; i < homolist.size(); i++)
//	{
//		tempLH = tempLH * homolist[i];
//		CalcCorners(tempLH, imagemaps[medpos]->image.size(), corners);//
//		maxsizelist[i] = corners.right_bottom.x;//
//	}
//
//	CalcCorners(tempLH, imagemaps[medpos]->image.size(), corners);
//	int maxRight = corners.right_bottom.x;
//
//	tempLH = Mat::eye(3, 3, CV_64FC1);
//	for (int i = medpos-1; i >=0; i--)
//	{
//		tempLH = tempLH * homolist[i];
//		CalcCorners(tempLH, imagemaps[medpos]->image.size(), corners);//
//		maxsizelist[i] = corners.left_bottom.x;/* imagemaps[medpos]->image.cols + abs(corners.left_bottom.x)*/
//	}
//
//	CalcCorners(tempLH, imagemaps[medpos]->image.size(), corners);
//	int minleft = corners.left_bottom.x;
//
//	/*Mat tempresult;
//	Point center(xtmpmask.cols / 2, xtmpmask.rows / 2);
//	seamlessClone(xtmp, result, xtmpmask, center, tempresult, MONOCHROME_TRANSFER);
//	result = tempresult.clone();*/
//	
//	tempLH = Mat::eye(3, 3, CV_64FC1);
//	Mat mymask(imagemaps[medpos]->image.size(), CV_8U, Scalar(255));//全白图
//	Mat xtmp;
//	Mat xtmpmask;
//	//xtmp = imagemaps[medpos]->image;//.clone()
//	//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\xtmp0.jpg", xtmp);
//
//	//for (int i = sumRadial-1; i >= 0; i--)
//	//{
//	//	imagemaps[i]->image.convertTo(imagemaps[i]->image, CV_32FC3);
//	//}
//
//	 //3->5图
//	for (int i = medpos; i < homolist.size(); i++)
//	{
//		////
//		//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\imagemaps0.jpg", imagemaps[0]->image);
//
//		tempLH = tempLH * homolist[i];
//		//tempLH = tempLH / tempLH.at<double>(2, 2);
//		//std::cout << "H:" << i << ":" << std::endl;
//		//std::cout << tempLH << std::endl;
//		//four_corners_t corners;
//		//CalcCorners(tempLH, imagemaps[i + 1]->image.size(), corners);
//		//图像矫正
//		Mat result;
//		Mat resultmask;
//		result = imagemaps[i + 1]->image.clone();
//		//////
//		//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\xtmp0.jpg", xtmp);
//		//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\result0.jpg", result);
//
//		//修正
//		//if (i == medpos)
//		//{
//		//	Color_correction(&imagemaps[i]->image/*xtmp*/, &imagemaps[i + 1]->image/*result*/);
//		//}
//		//else
//		//{
//		//	Color_correction(&xtmp, &imagemaps[i + 1]->image/*result*/);
//		//}
//
//		//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\xtmp1.jpg", xtmp);
//		//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\result1.jpg", result);
//		//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\mymask1.jpg", mymask);
//		//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\resutmask1.jpg", resultmask);
//		warpPerspective(imagemaps[i + 1]->image, result, tempLH, Size(maxsizelist[i]/*maxRight*/, imagemaps[i]->image.rows));
//		warpPerspective(mymask, resultmask, tempLH, Size(maxsizelist[i]/*maxRight*/, mymask.rows));
//		//if(i==4)cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\mymask.jpg", mymask);
//		//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\resultmask.jpg", resultmask);
//		//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\result.jpg", result);
//
//		//Color_correction(&xtmp, &result);
//
//		if (i == medpos)
//		{
//			Mat half(result, cv::Rect(0, 0, imagemaps[i]->image.cols, imagemaps[i]->image.rows));
//			Mat halfm(resultmask, cv::Rect(0, 0, imagemaps[i]->image.cols, imagemaps[i]->image.rows));
//
//
//			imagemaps[i]->image.copyTo(half);
//			mymask.copyTo(halfm);
//			xtmpmask = resultmask.clone();
//			//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\half.jpg", half);
//			//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\halfm.jpg", halfm);
//			//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\xtmpmask.jpg", xtmpmask);//
//			//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\result2.jpg", result);
//		}
//		else
//		{
//			//char namebuf[256];
//			//sprintf(namebuf, "%soutimg%.2dori_new.jpg", outputpath.c_str(), i);
//			//cv::imwrite(namebuf, resultmask);
//			//sprintf(namebuf, "%soutimg%.2dori.jpg", outputpath.c_str(), i);			
//			//cv::imwrite(namebuf, xtmpmask);
//			Mat element2 = getStructuringElement(MORPH_RECT, Size(5, 5));
//			cv::erode(xtmpmask, xtmpmask, element2);
//			//
//			Mat half(result, cv::Rect(0, 0, maxsizelist[i-1], imagemaps[i]->image.rows));//
//			Mat halfm(resultmask, cv::Rect(0, 0, maxsizelist[i-1], imagemaps[i]->image.rows));//
//			xtmp.copyTo(half);
//			xtmpmask.copyTo(halfm);
//			xtmpmask = resultmask.clone();
//
//			//这里亮度平衡调整？
//			//图像数据位置？？
//			//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\xtmp0.jpg", xtmp);
//			//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\result0.jpg", result);
//			//Color_correction(&xtmp, &result);
//			//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\xtmp1.jpg", xtmp);
//			//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\result1.jpg", result);
//
//			//原
//			//xtmp.copyTo(result, xtmpmask);
//			//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\result.jpg", result);
//			///*Mat blendout;
//			//Point center(0, 0);
//			//seamlessClone(xtmp, result, xtmpmask, center, blendout, NORMAL_CLONE);
//			//result = blendout.clone();*/			
//			//resultmask.copyTo(xtmpmask, resultmask);
//			//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\xtmpmask.jpg", xtmpmask);
//		}
//		//变行后图的输入目录
//		//char namebuf[256];
//		//sprintf(namebuf, "%soutimg%.2d.jpg", outputpath.c_str(), i);
//		//cv::imwrite(namebuf, result);
//		//sprintf(namebuf, "%soutimg%.2dmask.jpg", outputpath.c_str(), i);
//		//cv::imwrite(namebuf, resultmask);
//		xtmp = result.clone();
//		//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\xtmp345.jpg", xtmp);
//		
//	}
//
//	Mat rightimg = xtmp.clone();
//
//	tempLH = Mat::eye(3, 3, CV_64FC1);
//	xtmp = imagemaps[medpos]->image.clone();
//	int imgmaxwidth = imagemaps[medpos]->image.cols + abs(minleft);
//	tempLH.at<double>(0, 2)= abs(minleft);//
//  //0->3图
//	for (int i = medpos - 1; i >= 0; i--)
//	{
//		//if (i == 2)
//		//{
//		//	tempLH.at<double>(0, 2) = abs(maxsizelist[2]);
//		//}
//		//else
//		//{
//		//	tempLH.at<double>(0, 2) = imagemaps[3]->image.cols + abs(maxsizelist[i]) /*- imagemaps[i]->image.cols*/;
//		//}
//		tempLH = tempLH * homolist[i];
//
//		//tempLH = tempLH / tempLH.at<double>(2, 2);
//		//std::cout << "H:" << i << ":" << std::endl;
//		//std::cout << tempLH << std::endl;	
//		//图像矫正
//		Mat_<Vec3b> result;//数据类型注意调整Vec3b//Vec3f
//		Mat resultmask;
//
//		if (i == medpos - 1)
//		{
//			//Color_correction(&imagemaps[i]->image/*xtmp*/, &imagemaps[i + 1]->image/*result*/);
//
//		}
//		else
//		{
//			if (i == 0)
//			{
//				Mat temp_correction(xtmp, cv::Rect(xtmp.cols - abs(maxsizelist[i + 1]) - imagemaps[medpos]->image.cols, 0, abs(maxsizelist[i + 1]) + imagemaps[medpos]->image.cols, imagemaps[medpos]->image.rows));
//				Mat tempz = temp_correction.clone();
//				//Color_correction(&imagemaps[i]->image, &tempz);
//				tempz.copyTo(temp_correction);
//			}
//			else
//			{
//				Mat temp_correction(xtmp, cv::Rect(xtmp.cols - abs(maxsizelist[i + 1]) - imagemaps[medpos]->image.cols, 0, abs(maxsizelist[i + 1]) + imagemaps[medpos]->image.cols, imagemaps[medpos]->image.rows));
//				Mat tempz = temp_correction.clone();
//				//Color_correction(&imagemaps[i]->image, &tempz);
//				tempz.copyTo(temp_correction);
//			}
//			//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\image.jpg", imagemaps[i]->image);
//			//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\temp_correection.jpg", temp_correction);
//
//		}
//
//	    warpPerspective(imagemaps[i]->image, result, tempLH, Size(/*imagemaps[3]->image.cols + abs(maxsizelist[i])*/imgmaxwidth, imagemaps[i]->image.rows));
//		warpPerspective(mymask, resultmask, tempLH, Size(/*imagemaps[3]->image.cols + abs(maxsizelist[i])*/imgmaxwidth, mymask.rows));
//
//		//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\xtmp.jpg", xtmp);
//		//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\result.jpg", result);
//		//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\mymask.jpg", mymask);
//		//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\resutmask.jpg", resultmask);
//
//		if (i == medpos - 1)
//		{
//			//Mat tmp = result.clone();
//			//result = tmp.clone();
//			Mat half(result, cv::Rect(/*abs(maxsizelist[i])*/abs(minleft), 0, imagemaps[i]->image.cols, imagemaps[i]->image.rows));
//			Mat halfm(resultmask, cv::Rect(/*abs(maxsizelist[i])*/abs(minleft), 0, imagemaps[i]->image.cols, imagemaps[i]->image.rows));
//			//Mat_<Vec3b> tmp1;
//			//imagemaps[i + 1]->image.convertTo(tmp1, result.type());
//			//tmp1.copyTo(half);
//			imagemaps[i + 1]->image.copyTo(half);
//			mymask.copyTo(halfm);
//			xtmpmask = resultmask.clone();
//			//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\result.jpg", result);
//			//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\half.jpg", half);
//			//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\imagemapi+1.jpg", imagemaps[i + 1]->image);
//			//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\xtmpmask.jpg", xtmpmask);//
//		}
//		else
//		{
//			//char namebuf[256];
//			//sprintf(namebuf, "%soutimg%.2dori_new.jpg", outputpath.c_str(), i);
//			//cv::imwrite(namebuf, resultmask);
//			//sprintf(namebuf, "%soutimg%.2dori.jpg", outputpath.c_str(), i);
//			//cv::imwrite(namebuf, xtmpmask);
//			Mat element2 = getStructuringElement(MORPH_RECT, Size(5, 5));
//			cv::erode(xtmpmask, xtmpmask, element2);
//			//
//			Mat half(result, cv::Rect(xtmp.cols - abs(maxsizelist[i+1]) - imagemaps[medpos]->image.cols, 0, abs(maxsizelist[i+1]) + imagemaps[medpos]->image.cols, imagemaps[medpos]->image.rows));
//			Mat halfm(resultmask, cv::Rect(xtmp.cols - abs(maxsizelist[i+1]) - imagemaps[medpos]->image.cols, 0, abs(maxsizelist[i+1]) + imagemaps[medpos]->image.cols, imagemaps[medpos]->image.rows));
//			Mat temp_correction(xtmp, cv::Rect(xtmp.cols - abs(maxsizelist[i+1]) - imagemaps[medpos]->image.cols, 0, abs(maxsizelist[i+1]) + imagemaps[medpos]->image.cols, imagemaps[medpos]->image.rows));
//			temp_correction.copyTo(half);
//			xtmpmask.copyTo(halfm);
//			xtmpmask = resultmask.clone();
//
//
//			//xtmp.copyTo(result, xtmpmask);
//			//resultmask.copyTo(xtmpmask, resultmask);
//
//			//
//			//Mat half(result, cv::Rect(0/*abs(maxsizelist[i] - maxsizelist[i + 1])*/, 0, imagemaps[i]->image.cols, imagemaps[i]->image.rows));
//			//Mat halfm(resultmask, cv::Rect(0/*abs(maxsizelist[i] - maxsizelist[i + 1])*/, 0, imagemaps[i]->image.cols, imagemaps[i]->image.rows));
//			//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\xtmpmask0.jpg", xtmpmask);//
//			//xtmp.copyTo(half);
//			//xtmpmask.copyTo(halfm);
//			//xtmpmask = resultmask.clone();
//			//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\half1.jpg", half);
//			//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\halfm1.jpg", halfm);
//			//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\xtmpmask1.jpg", xtmpmask);//
//		}
//
//		//变行后图的输入目录
//		//char namebuf[256];
//		//sprintf(namebuf, "%soutimg%.2d.jpg", outputpath.c_str(), i);
//		//cv::imwrite(namebuf, result);
//		//sprintf(namebuf, "%soutimg%.2dmask.jpg", outputpath.c_str(), i);
//		//cv::imwrite(namebuf, resultmask);
//		xtmp = result.clone();
//		//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\xtmp0123.jpg", xtmp);
//	}
//	
//	Radialout = new ImageOut();
//	imgmaxwidth = rightimg.cols + xtmp.cols - imagemaps[medpos]->image.cols;
//	Radialout->outimage = Mat( rightimg.rows, imgmaxwidth, rightimg.type());
//
//	Mat temp_correction(xtmp, cv::Rect(0, 0, abs(maxsizelist[0])+100, rightimg.rows));
//	Mat tempz = temp_correction.clone();
//	//Color_correction(&tempz, &rightimg);
//	tempz.copyTo(temp_correction);
//
//	lmedleft = imgmaxwidth - rightimg.cols;
//	lwidth = imgmaxwidth;
//	lheight = rightimg.rows;
//
//	Mat halfleft(Radialout->outimage, cv::Rect(0, 0, xtmp.cols, xtmp.rows));
//	Mat halfright(Radialout->outimage, cv::Rect(imgmaxwidth - rightimg.cols, 0, rightimg.cols, rightimg.rows));
//
//	xtmp.copyTo(halfleft);
//	rightimg.copyTo(halfright);
//
//	sprintf(namebuf, "%sout\\%.2d.jpg", outputpath.c_str(),ipos);
//	cv::imwrite(namebuf, Radialout->outimage);
//
//}
//
//void TunnelMosaic::runRadial23(std::vector<std::string>  imgpathlist, int ipos)
//{
//
//	std::vector<Mat> imgsets;
//	char namebuf[256];
//	bool withRotation = false;
//	bool withScale = false;
//
//	Ptr<DescriptorMatcher> matcher = DescriptorMatcher::create("BruteForce-Hamming");
//	medpos = sumRadial / 2;
//
//	std::vector<Mat>  homolist;
//
//
//
//	for (int i = 0; i < sumRadial; i++)
//	{
//		Mat frame = cv::imread(imgpathlist[i]);
//		if (frame.empty()) continue;
//
//		//1 标定矫正，是否需要?
//		//Mat dedistmat = imagemaps[i]->deDistort(frame);
//		//if (dedistmat.empty()) continue;
//
//		//2 柱面投影? 内投和外投
//		//Mat outmat = imagemaps[i]->tocylinder(dedistmat);
//		//if (outmat.empty()) continue;	
//		//imgsets.push_back(outmat);
//		//imgsets.push_back(frame);
//
//		//3 调整整体亮度？
//
//		//////////////////////
//		//颜色校正？
//
//
//
//
//
//		//4 提取特征点
//		imagemaps[i]->getFeaturePoints(-1, frame);
//
//		//kpRef, descRef
//		//保存特征点到目标文件，，用于轴向的拼接以及周向的拼接
//		sprintf(namebuf, "%s%d\\outimgfeat%.2d.yaml", outputpath.c_str(), i, ipos);
//		FileStorage fs(namebuf, FileStorage::WRITE);
//		fs << "feature_points" << imagemaps[i]->kpRef;
//		fs << "points_description" << imagemaps[i]->descRef;
//		fs.release();
//		//特征点输出结束
//		//变形后的图像加入列表，用于轴向的拼接
//		imagemaps[i]->image = frame.clone();
//		Mat outmat;
//		//outmat = frame.clone();
//		//drawKeypoints(frame, imagemaps[i]->kpRef, outmat);		
//		//投影后图的输入目录
//		//sprintf(namebuf, "%s%d\\outimg%.2d.jpg", outputpath.c_str(), i, ipos);
//		//cv::imwrite(namebuf, frame);
//
//		if (i > 0)
//		{
//			//5 环向的拼接进行匹配
//			std::vector<DMatch> matchesAll, matchesGMS;
//			matcher->match(imagemaps[i]->descRef, imagemaps[i - 1]->descRef, matchesAll);
//
//			//精确匹配
//			matchGMS(imagemaps[i]->image.size(), imagemaps[i - 1]->image.size(),
//				imagemaps[i]->kpRef, imagemaps[i - 1]->kpRef,
//				matchesAll, matchesGMS, withRotation, withScale);
//
//			//绘制匹配  用于调试
//			//Mat dst; 
//			//drawMatches(imagemaps[i  ]->image, imagemaps[i    ]->kpRef, 
//			//	        imagemaps[i-1]->image, imagemaps[i - 1]->kpRef, matchesGMS, dst);
//			//
//
//			//6 计算单应矩阵
//			// 将精确匹配点取出，imagePoints1为第一幅图片，imagePoints2为第二幅图像
//			std::vector<Point2f> imagePoints1, imagePoints2;
//			for (int j = 0; j < matchesGMS.size(); j++)
//			{
//				imagePoints2.push_back(imagemaps[i]->kpRef[matchesGMS[j].queryIdx].pt);//i
//				imagePoints1.push_back(imagemaps[i - 1]->kpRef[matchesGMS[j].trainIdx].pt);//i-1
//			}
//
//
//			Mat homo;
//
//			if (i <= medpos)
//			{
//				homo = findHomography(imagePoints1, imagePoints2, cv::RANSAC);
//
//			}
//			else
//			{
//				homo = findHomography(imagePoints2, imagePoints1, cv::RANSAC);
//			}
//			//invert(h12, h21, DECOMP_LU);
//			homolist.push_back(homo);
//		}
//	}
//
//	Mat result;
//	four_corners_t corners;
//
//	//
//	vector<int> maxsizelist(sumRadial, 0);//-1
//
//	//计算拼接后图像的尺寸，以medpos图像为基准	
//	Mat tempLH = Mat::eye(3, 3, CV_64FC1);
//	for (int i = medpos; i < homolist.size(); i++)
//	{
//		tempLH = tempLH * homolist[i];
//		CalcCorners(tempLH, imagemaps[medpos]->image.size(), corners);//
//		maxsizelist[i + 1] = corners.right_bottom.x;//
//	}
//
//	CalcCorners(tempLH, imagemaps[medpos]->image.size(), corners);
//	int maxRight = corners.right_bottom.x;
//
//	tempLH = Mat::eye(3, 3, CV_64FC1);
//	for (int i = medpos - 1; i >= 0; i--)
//	{
//		tempLH = tempLH * homolist[i];
//		CalcCorners(tempLH, imagemaps[medpos]->image.size(), corners);//
//		maxsizelist[i] = corners.left_bottom.x;/* imagemaps[medpos]->image.cols + abs(corners.left_bottom.x)*/
//	}
//
//	//换初始点
//	/*maxsizelist[medpos] = maxsizelist[medpos] - maxsizelist[0];*/
//	//for (int i = sumRadial - 1; i > medpos; i--)
//	//{
//	//	maxsizelist[i] = maxsizelist[i] - maxsizelist[0] - imagemaps[i]->image.cols;
//	//}
//	for (int i = medpos;i >= 0; i--)
//	{
//		maxsizelist[i]=maxsizelist[i] - maxsizelist[0];
//	}
//
//
//	CalcCorners(tempLH, imagemaps[medpos]->image.size(), corners);
//	int minleft = corners.left_bottom.x;
//
//	/*Mat tempresult;
//	Point center(xtmpmask.cols / 2, xtmpmask.rows / 2);
//	seamlessClone(xtmp, result, xtmpmask, center, tempresult, MONOCHROME_TRANSFER);
//	result = tempresult.clone();*/
//
//	tempLH = Mat::eye(3, 3, CV_64FC1);
//	Mat mymask(imagemaps[medpos]->image.size(), CV_8U, Scalar(255));//全白图
//	Mat xtmp;
//	Mat xtmpmask;
//	int imgmaxwidth = imagemaps[medpos]->image.cols + abs(minleft);
//	tempLH.at<double>(0, 2) = abs(minleft);
//
//	//计算序列
//	vector<Mat> tempLHlist(5);
//	for (int i = medpos - 1; i >= 0; i--)
//	{
//		tempLHlist[i]= tempLH * homolist[i];
//		tempLH = tempLH * homolist[i];
//	}
//	tempLH = Mat::eye(3, 3, CV_64FC1);
//	//tempLH.at<double>(0, 2) = maxRight;
//	for (int i = medpos; i < sumRadial - 1; i++)
//	{
//		tempLHlist[i] = tempLH * homolist[i];
//		tempLH = tempLH * homolist[i];
//	}
//
//	 ////0->3
//	for (int i = 0; i <= medpos - 1; i++)
//	{
//		//tempLH = tempLH * homolist[i];
//		Mat_<Vec3b> result;//数据类型注意调整Vec3b//Vec3f
//		Mat resultmask;
//
//		if (i == 0)
//		{
//			Color_correction(&imagemaps[i]->image/*xtmp*/, &imagemaps[i + 1]->image/*result*/);
//
//			warpPerspective(imagemaps[i]->image, result, tempLHlist[i], Size(maxsizelist[medpos] + imagemaps[medpos]->image.cols, imagemaps[i]->image.rows));
//			xtmp = result.clone();
//
//			warpPerspective(imagemaps[i + 1]->image, result, tempLHlist[i + 1], Size(maxsizelist[medpos] + imagemaps[medpos]->image.cols, imagemaps[i]->image.rows));
//			warpPerspective(mymask, resultmask, tempLHlist[i + 1], Size(maxsizelist[medpos] + imagemaps[medpos]->image.cols, imagemaps[i]->image.rows));
//			//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\result.jpg", result);
//			//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\resultmask.jpg", resultmask);
//			Mat half(xtmp, cv::Rect(maxsizelist[i + 1], 0, imagemaps[i + 1]->image.cols, imagemaps[i]->image.rows));
//			Mat half1(result, cv::Rect(maxsizelist[i + 1], 0, imagemaps[i + 1]->image.cols, imagemaps[i]->image.rows));
//			Mat halfm(resultmask, cv::Rect(maxsizelist[i + 1], 0, imagemaps[i + 1]->image.cols, imagemaps[i]->image.rows));
//			half1.copyTo(half, halfm);
//			//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\xtmp.jpg", xtmp);
//		}
//		else
//		{
//			Mat temp_correction(xtmp, cv::Rect(0, 0, maxsizelist[i]+imagemaps[i]->image.cols - 20, imagemaps[medpos]->image.rows));
//			Mat tempz = temp_correction.clone();
//			Color_correction(&tempz, &imagemaps[i + 1]->image);
//			tempz.copyTo(temp_correction);
//			//////Color_correction(&xtmp, &imagemaps[i + 1]->image/*result*/);;////待续
//
//			if (i == medpos - 1)
//			{
//				Mat half(xtmp, cv::Rect(maxsizelist[i + 1], 0, imagemaps[i + 1]->image.cols, imagemaps[i]->image.rows));
//				imagemaps[i + 1]->image.copyTo(half);
//			}
//			else
//			{
//				warpPerspective(imagemaps[i + 1]->image, result, tempLHlist[i + 1], Size(maxsizelist[sumRadial - 1] + imagemaps[medpos]->image.cols, imagemaps[i]->image.rows));
//				warpPerspective(mymask, resultmask, tempLHlist[i + 1], Size(maxsizelist[sumRadial - 1] + imagemaps[medpos]->image.cols, imagemaps[i]->image.rows));
//				//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\result.jpg", result);
//				//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\resultmask.jpg", resultmask);
//				Mat half(xtmp, cv::Rect(maxsizelist[i + 1], 0, imagemaps[i + 1]->image.cols, imagemaps[i]->image.rows));
//				Mat half1(result, cv::Rect(maxsizelist[i + 1], 0, imagemaps[i + 1]->image.cols, imagemaps[i]->image.rows));
//				Mat halfm(resultmask, cv::Rect(maxsizelist[i + 1], 0, imagemaps[i + 1]->image.cols, imagemaps[i]->image.rows));
//				half1.copyTo(half, halfm);
//			}
//		}
//
//		//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\xtmp.jpg", xtmp);
//	}
//
//
//	//// 4->5
//	for (int i = medpos; i < sumRadial - 1; i++)
//	{
//		Mat_<Vec3b> result;//数据类型注意调整Vec3b//Vec3f
//		Mat resultmask;
//
//		if (i == medpos)
//		{
//			Mat temp_correction(xtmp, cv::Rect(0, 0, xtmp.cols - 20, imagemaps[medpos]->image.rows));
//			Mat tempz = temp_correction.clone();
//			Color_correction(&tempz, &imagemaps[i + 1]->image);
//			tempz.copyTo(temp_correction);
//		}
//		else {
//			Mat temp_correction(xtmp, cv::Rect(0, 0, xtmp.cols - 20, imagemaps[medpos]->image.rows));
//			Mat tempz = temp_correction.clone();
//			Color_correction(&tempz, &imagemaps[i + 1]->image);
//			tempz.copyTo(temp_correction);
//		}
//
//
//		//Color_correction(&xtmp, &imagemaps[i + 1]->image);
//
//		warpPerspective(imagemaps[i + 1]->image, result, tempLHlist[i], Size(maxsizelist[i + 1], imagemaps[i]->image.rows));
//		warpPerspective(mymask, resultmask, tempLHlist[i], Size(maxsizelist[i + 1], imagemaps[i]->image.rows));
//		//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\result.jpg", result);
//		//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\resultmask.jpg", resultmask);
//		if (i == medpos)
//		{
//			cv::copyMakeBorder(xtmp, xtmp, 0, 0, 0, maxsizelist[i + 1] - imagemaps[medpos]->image.cols, CV_8UC3, cv::Scalar(0));
//		}
//		else{
//			cv::copyMakeBorder(xtmp, xtmp, 0, 0, 0, maxsizelist[i + 1] - maxsizelist[i], CV_8UC3, cv::Scalar(0));
//		}
//
//		//Mat element2 = getStructuringElement(MORPH_RECT, Size(5, 5)); 
//		//cv::erode(resultmask, resultmask, element2);
//
//		Mat half(xtmp, cv::Rect(xtmp.cols - imagemaps[i]->image.cols+20, 0, imagemaps[i]->image.cols-20, imagemaps[i]->image.rows));
//		Mat half1(result, cv::Rect(result.cols - imagemaps[i]->image.cols+20, 0, imagemaps[i]->image.cols-20, imagemaps[i]->image.rows));
//		Mat halfm(resultmask, cv::Rect(resultmask.cols - imagemaps[i]->image.cols+20, 0, imagemaps[i]->image.cols-20, imagemaps[i]->image.rows));
//		//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\half.jpg", half);
//		//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\half1.jpg", half1);
//		//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\halfm.jpg", halfm);
//		half1.copyTo(half,halfm);
//		//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\xtmp.jpg", xtmp);
//	}
//
//	Radialout = new ImageOut();
//	
//	Radialout->outimage = xtmp.clone();
//
//	sprintf(namebuf, "%sout\\%.2d.jpg", outputpath.c_str(), ipos);
//	cv::imwrite(namebuf, Radialout->outimage);
//
//}
//
//
//
//void TunnelMosaic::runAxial(TunnelMosaic& preRadia)
//{
//
//	int medpos = preRadia.medpos;
//
//	Ptr<DescriptorMatcher> matcher = DescriptorMatcher::create("BruteForce-Hamming");
//	std::vector<DMatch> matchesAll, matchesGMS;
//	matcher->match(imagemaps[medpos]->descRef, preRadia.imagemaps[medpos]->descRef, matchesAll);
//
//	//精确匹配
//	matchGMS(imagemaps[medpos]->image.size(), preRadia.imagemaps[medpos]->image.size(),
//		imagemaps[medpos]->kpRef, preRadia.imagemaps[medpos]->kpRef,
//		matchesAll, matchesGMS, false, false);
//
//	std::vector<Point2f> imagePoints1, imagePoints2;
//	for (int j = 0; j < matchesGMS.size(); j++)
//	{
//		imagePoints2.push_back(imagemaps[medpos]->kpRef[matchesGMS[j].queryIdx].pt);//i
//		imagePoints1.push_back(preRadia.imagemaps[medpos]->kpRef[matchesGMS[j].trainIdx].pt);//i-1
//	}
//	Mat homo;
//	homo = findHomography(imagePoints2, imagePoints1, cv::RANSAC);
//
//	four_corners_t corners;
//	CalcCorners(homo, imagemaps[medpos]->image.size(), corners);
//
//	Mat result;
//	warpPerspective(Radialout->outimage, result, homo, Size(preRadia.Radialout->outimage.cols, Radialout->outimage.rows + corners.left_top.y));
//
//	Mat halfleft(result, cv::Rect(0, 0, preRadia.Radialout->outimage.cols, preRadia.Radialout->outimage.rows));
//	preRadia.Radialout->outimage.copyTo(halfleft);
//
//	char namebuf[256];
//	sprintf(namebuf, "%stwoutwarp.jpg", outputpath.c_str());
//	cv::imwrite(namebuf, result);
//
//}

void TunnelMosaic::OptimizeSeam()
{
	//POSSION泊松融合
	//MASK白色为src，det为黑色
	//seamlessClone(SRC, DET, mask, center, blendImg, NORMAL_CLONE);

}

class LaplacianBlending {
private:
    Mat left;
    Mat right;
    Mat blendMask;

    //Laplacian Pyramids
    std::vector<Mat> leftLapPyr, rightLapPyr, resultLapPyr;
    Mat leftHighestLevel, rightHighestLevel, resultHighestLevel;
    //mask为三通道方便矩阵相乘
	std::vector<Mat> maskGaussianPyramid;

    int levels;

    void buildPyramids()
    {
        buildLaplacianPyramid(left, leftLapPyr, leftHighestLevel);
        buildLaplacianPyramid(right, rightLapPyr, rightHighestLevel);
        buildGaussianPyramid();
    }

    void buildGaussianPyramid()
    {
        //金字塔内容为每一层的掩模
        assert(leftLapPyr.size()>0);

        maskGaussianPyramid.clear();
        Mat currentImg;
        cvtColor(blendMask, currentImg, cv::COLOR_GRAY2BGR);
        //保存mask金字塔的每一层图像
        maskGaussianPyramid.push_back(currentImg); //0-level

        currentImg = blendMask;
        for (int l = 1; l<levels + 1; l++) {
            Mat _down;
            if (leftLapPyr.size() > l)
                pyrDown(currentImg, _down, leftLapPyr[l].size());
            else
                pyrDown(currentImg, _down, leftHighestLevel.size()); //lowest level

            Mat down;
            cvtColor(_down, down, COLOR_GRAY2BGR);
            //add color blend mask into mask Pyramid
            maskGaussianPyramid.push_back(down);
            currentImg = _down;
        }
    }

    void buildLaplacianPyramid(const Mat& img, std::vector<Mat>& lapPyr, Mat& HighestLevel)
    {
        lapPyr.clear();
        Mat currentImg = img;
        for (int l = 0; l<levels; l++) {
            Mat down, up;
            pyrDown(currentImg, down);
            pyrUp(down, up, currentImg.size());
            Mat lap = currentImg - up;
            lapPyr.push_back(lap);
            currentImg = down;
        }
        currentImg.copyTo(HighestLevel);
    }

    Mat reconstructImgFromLapPyramid()
    {
        //将左右laplacian图像拼成的resultLapPyr金字塔中每一层
        //从上到下插值放大并与残差相加，即得blend图像结果
        Mat currentImg = resultHighestLevel;
        for (int l = levels - 1; l >= 0; l--)
        {
            Mat up;
            pyrUp(currentImg, up, resultLapPyr[l].size());
            currentImg = up + resultLapPyr[l];
        }
        return currentImg;
    }

    void blendLapPyrs()
    {
        //获得每层金字塔中直接用左右两图Laplacian变换拼成的图像resultLapPyr
        resultHighestLevel = leftHighestLevel.mul(maskGaussianPyramid.back()) +
            rightHighestLevel.mul(Scalar(1.0, 1.0, 1.0) - maskGaussianPyramid.back());
        for (int l = 0; l<levels; l++)
        {
            Mat A = leftLapPyr[l].mul(maskGaussianPyramid[l]);
            Mat antiMask = Scalar(1.0, 1.0, 1.0) - maskGaussianPyramid[l];
            Mat B = rightLapPyr[l].mul(antiMask);
            Mat blendedLevel = A + B;

            resultLapPyr.push_back(blendedLevel);
        }
    }

public:
	//construct function, used in LaplacianBlending lb(l,r,m,4);
    LaplacianBlending(const Mat& _left, const Mat& _right, const Mat& _blendMask, int _levels) :
        left(_left), right(_right), blendMask(_blendMask), levels(_levels)
    {
        assert(_left.size() == _right.size());
        assert(_left.size() == _blendMask.size());
        //创建拉普拉斯金字塔和高斯金字塔
        buildPyramids();
        //每层金字塔图像合并为一个
        blendLapPyrs();
    };

    Mat blend()
    {
		//重建拉普拉斯金字塔
        return reconstructImgFromLapPyramid();
    }
};

Mat LaplacianBlend(const Mat &left, const Mat &right, const Mat &mask)
{
    LaplacianBlending laplaceBlend(left, right, mask, 10);
    return laplaceBlend.blend();
}

//void TunnelMosaic::

void TunnelMosaic::my_stitching(std::vector<std::string>  imgpathlist, Mat *dst)
{
	std::vector<Mat> imglist;
	std::vector<Mat> outimglist;
	
	Mat tempimg;

	Mat frame;
	//柱面校正 -> 拼接 -> 切割 -> 亮度调整(另)
	//校正
	for (size_t i = 0; i < imgpathlist.size(); i++)
	{
		frame = cv::imread(imgpathlist[i]);
		//cout << "imgpathlist:" << imgpathlist[i] << endl;

		if (frame.empty()) continue;

		if (i == 0) {
			Mat srcimg, dstimg;
			srcimg = frame.clone();

			//
			//char input_win[] = "input image"; //定义一个变量input_win来表示窗口名
			//cvtColor(src, src, CV_BGR2GRAY);
			//namedWindow(input_win, CV_WINDOW_AUTOSIZE);
			//imshow(input_win, src);

			// contrast and brigthtness changes
			int height = srcimg.rows;
			int width = srcimg.cols;

			//cout << "height" << src.rows << endl;

			dstimg = Mat::zeros(srcimg.size(), srcimg.type());
			float alpha = 1.8;
			float beta = 5;

			Mat m1;
			srcimg.convertTo(m1, CV_32F); //强制转化浮点型，提高精度
			for (int row = 0; row < height; row++) {
				for (int col = 0; col < width; col++) {
					if (srcimg.channels() == 3) {
						float b = m1.at<Vec3f>(row, col)[0];// blue
						float g = m1.at<Vec3f>(row, col)[1]; // green
						float r = m1.at<Vec3f>(row, col)[2]; // red

						// output
						dstimg.at<Vec3b>(row, col)[0] = saturate_cast<uchar>(b * alpha + beta);
						dstimg.at<Vec3b>(row, col)[1] = saturate_cast<uchar>(g * alpha + beta);
						dstimg.at<Vec3b>(row, col)[2] = saturate_cast<uchar>(r * alpha + beta);
					}
					else if (srcimg.channels() == 1) {
						float v = srcimg.at<uchar>(row, col);
						dstimg.at<uchar>(row, col) = saturate_cast<uchar>(v * alpha + beta);
					}
				}
			}

			frame = dstimg.clone();


		}
		else {
			if (i == imgpathlist.size() - 1) {
				Mat srcimg, dstimg;
				srcimg = frame.clone();

				//
				//char input_win[] = "input image"; //定义一个变量input_win来表示窗口名
				//cvtColor(src, src, CV_BGR2GRAY);
				//namedWindow(input_win, CV_WINDOW_AUTOSIZE);
				//imshow(input_win, src);

				// contrast and brigthtness changes
				int height = srcimg.rows;
				int width = srcimg.cols;

				//cout << "height" << src.rows << endl;

				dstimg = Mat::zeros(srcimg.size(), srcimg.type());
				float alpha = 1.8;
				float beta = 5;

				Mat m1;
				srcimg.convertTo(m1, CV_32F); //强制转化浮点型，提高精度
				for (int row = 0; row < height; row++) {
					for (int col = 0; col < width; col++) {
						if (srcimg.channels() == 3) {
							float b = m1.at<Vec3f>(row, col)[0];// blue
							float g = m1.at<Vec3f>(row, col)[1]; // green
							float r = m1.at<Vec3f>(row, col)[2]; // red

							// output
							dstimg.at<Vec3b>(row, col)[0] = saturate_cast<uchar>(b * alpha + beta);
							dstimg.at<Vec3b>(row, col)[1] = saturate_cast<uchar>(g * alpha + beta);
							dstimg.at<Vec3b>(row, col)[2] = saturate_cast<uchar>(r * alpha + beta);
						}
						else if (srcimg.channels() == 1) {
							float v = srcimg.at<uchar>(row, col);
							dstimg.at<uchar>(row, col) = saturate_cast<uchar>(v * alpha + beta);
						}
					}
				}
				Mat TEMP;
				dstimg.convertTo(TEMP, frame.type());
				frame = dstimg.clone();
			}

		}

		//cv::imshow("111", frame);
		//cv::waitKey(1);
		cylindrical_projection(&frame);



		imglist.push_back(frame);

		//cout << "imglist:" << imglist.size() << endl;
		frame.release();

	}

	


	//拼接
	if (!imglist.empty()) 
	{
		//cout << "imgpathlist:" << imgpathlist.size() << endl;

		if (imglist.size() >= 2) {
			my_stitching_method(&imglist[0], &imglist[1], &tempimg, 0, 2, 1);
		}
		if (imglist.size() >= 3) {
			my_stitching_method(&tempimg, &imglist[2], &tempimg, -61, 11, 1);
		}
		if (imglist.size() >= 4) {
			my_stitching_method(&tempimg, &imglist[3], &tempimg, 0, 0, 1);
		}
		if (imglist.size() >= 5) {
			my_stitching_method(&tempimg, &imglist[4], &tempimg, -20, -30, 1);
		}
		if (imglist.size() >= 6) {
			my_stitching_method(&tempimg, &imglist[5], &tempimg, -230, -10, 1);
		}

		//cv::imwrite("out.jpg", tempimg);

		Mat cut(tempimg, Rect(0/*tempimg.cols * 0.12*/, tempimg.rows * 0.12, tempimg.cols /** 0.76*/, tempimg.rows * 0.76));

		*dst = cut.clone();

	}





}

void TunnelMosaic::start_stitching(string rawPath, string resultPath)
{


	string bootpath = rawPath;
	vector<string> imgpath = { "\\Camera1\\","\\Camera2\\","\\Camera3\\","\\Camera4\\","\\Camera5\\","\\Camera6\\" };
	string outimgpath = resultPath;

	vector<string> fileslist0;
	vector<string> fileslist1;
	vector<string> fileslist2;
	vector<string> fileslist3;
	vector<string> fileslist4;
	vector<string> fileslist5;

	//获取该路径下的图片文件  
	getFiles(bootpath + imgpath[0], ".jpg", fileslist0);
	getFiles(bootpath + imgpath[1], ".jpg", fileslist1);
	getFiles(bootpath + imgpath[2], ".jpg", fileslist2);
	getFiles(bootpath + imgpath[3], ".jpg", fileslist3);
	getFiles(bootpath + imgpath[4], ".jpg", fileslist4);
	getFiles(bootpath + imgpath[5], ".jpg", fileslist5);

	Mat outimg;


	for (size_t i = 0/*1235*/; i < fileslist0.size(); i++)
	{
		//cout << "正在处理第" << i + 1 << "张···" << endl;
		
		vector<string> names;
		names.push_back(bootpath + imgpath[0] + fileslist0[i]);
		names.push_back(bootpath + imgpath[1] + fileslist1[i]);
		names.push_back(bootpath + imgpath[2] + fileslist2[i]);
		names.push_back(bootpath + imgpath[3] + fileslist3[i]);
		names.push_back(bootpath + imgpath[4] + fileslist4[i]);
		names.push_back(bootpath + imgpath[5] + fileslist5[i]);

		/*for (int j = 0; j < names.size(); j++)
		{
			cout << names[j] << endl;
		}*/

		my_stitching(names, &outimg);

		//cv::imwrite(outimgpath + "\\" + fileslist[i], outimg);
		//string file = outimgpath + "\\" + fileslist[i];
		//sprintf(file, "%s\\%s", outimgpath, fileslist[i]);
		//if (!bmp.ReadBmp(file.c_str())) {
		//	cout << "文件打开失败！" << file << endl;
		//	continue;
		//}
		//cout << "正在进行亮度校正" << file << endl;
		//if (bmp.bi.biBitCount < 24) {
		//	bmp.GrayDrawHistogram();
		//	bmp.GrayHistogramEqualization();
		//	bmp.GrayDrawHistogram();
		//}
		//else {
		//	bmp.ColorDrawHistogram();
		//	bmp.ColorHistogramEqualization();
		//	bmp.ColorDrawHistogram();
		//}
		////sprintf_s(file, "%d-out.jpg", i);
		//bmp.SaveBmp(file.c_str());
		/*Mat dst(outimg.rows, outimg.cols, outimg.type());*/
		//cv::imwrite(outimgpath + "\\1.jpg", dst);
		//cout << "正在处理第" << i+1 << "张···" << endl;

		/*Mat src, dst;
		src = outimg.clone();*/
		
		////
		////char input_win[] = "input image"; //定义一个变量input_win来表示窗口名
		////cvtColor(src, src, CV_BGR2GRAY);
		////namedWindow(input_win, CV_WINDOW_AUTOSIZE);
		////imshow(input_win, src);

		//// contrast and brigthtness changes
		//int height = src.rows;
		//int width = src.cols;

		////cout << "height" << src.rows << endl;

		//dst = Mat::zeros(src.size(), src.type());
		//float alpha = 1.8;
		//float beta = 5;

		//Mat m1;
		//src.convertTo(m1, CV_32F); //强制转化浮点型，提高精度
		//for (int row = 0; row < height; row++) {
		//	for (int col = 0; col < width; col++) {
		//		if (src.channels() == 3) {
		//			float b = m1.at<Vec3f>(row, col)[0];// blue
		//			float g = m1.at<Vec3f>(row, col)[1]; // green
		//			float r = m1.at<Vec3f>(row, col)[2]; // red

		//			// output
		//			dst.at<Vec3b>(row, col)[0] = saturate_cast<uchar>(b * alpha + beta);
		//			dst.at<Vec3b>(row, col)[1] = saturate_cast<uchar>(g * alpha + beta);
		//			dst.at<Vec3b>(row, col)[2] = saturate_cast<uchar>(r * alpha + beta);
		//		}
		//		else if (src.channels() == 1) {
		//			float v = src.at<uchar>(row, col);
		//			dst.at<uchar>(row, col) = saturate_cast<uchar>(v * alpha + beta);
		//		}
		//	}
		//}
		////

		//char output_title[] = "contrast and brightness change demo";
		//namedWindow(output_title, CV_WINDOW_AUTOSIZE);

		//Mat out;
		//dst.convertTo(out, CV_8U);
		//char outpath[256];
		//sprintf(outpath, "%s\\%s", outimgpath, fileslist[i]);
		string outpath = outimgpath + "\\" + fileslist0[i];
		//cout << "r" << dst.rows << endl;

		//cout<<"r"<<out.rows << endl;
		//cv::imwrite(outpath, out);
		cv::imwrite(outpath, outimg);

		num = i + 1;
		//waitKey(0);

	}

	
}

void TunnelMosaic::startMosaicAlgorithm(string rawPath, string resultPath)
{
	setconfigpath(rawPath);
	settemppath(rawPath);
	setoutputpath(rawPath);
	//loadconfig();

	srcPath = rawPath;
	dstPath = resultPath;

	num = 0;

	start_stitching(rawPath, resultPath);

}

int TunnelMosaic::getCurrImageNum()
{
	if (dstPath.size() == 0) {
		cout << "输入路径缺失！" << endl;
		return 0;
	}

	return num;
}
