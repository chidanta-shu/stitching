// mosaicTest.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include <sstream>
#include <string>
#include <CameraCalib.h>
#include <TunnelMosaic.h>
//#include <dwrite.h>

using namespace std;

void calibTest()
{
	std::string cfgpath = "C:\\Users\\10048\\Desktop\\mosaic(1)\\bin\\x64\\config\\";
	int nid = 0;//相机通道编号
	std::ostringstream ostrm;
	ostrm << nid;
	std::string cfgfile = cfgpath + "calibcfg" + ostrm.str() + ".xml";

	////先生成标定配置文件
	////文件目录为E:\\rail\\mosaic\\bin\\x64\\config\\
	////可以针对某一相机拍摄的标定照片标定，照片路径为"
	CameraCalib::initcalibfile("C:\\Users\\10048\\Desktop\\mosaic(1)\\bin\\x64\\config\\", "C:\\Users\\10048\\Desktop\\mosaic(1)\\bin\\x64\\calibimages\\0\\", nid);
	////输出文件名为cfgfile
	////在这个基础上进行运行标定程序
	CameraCalib::calibmain(cfgfile);

}


int main(int argc, char **argv)
{

	
	string rawPath;
	rawPath = "E:\\隧道巡检-公司内部数据\\20210419_014016_41千米400米_RUD_15型_DFRD_2020_01_15";
	string resultPath;
	resultPath = "E:\\隧道巡检-公司内部数据\\20210419_014016_41千米400米_RUD_15型_DFRD_2020_01_15\\out";
	
	TunnelMosaic test;
	test.startMosaicAlgorithm(rawPath, resultPath);



	//std::string cfgpath = "C:\\Users\\10048\Desktop\\x64\\config\\";
	////原始图像总路径，通道分为1，2，3，4，5....
	//std::string	imagelistpath;

	//TunnelMosaic   myTunnelMosaic1, myTunnelMosaic2;

	//myTunnelMosaic1.setconfigpath(cfgpath);
	//myTunnelMosaic1.settemppath("E:\\Temp\\");
	//myTunnelMosaic1.setoutputpath("E:\\Temp\\");
	//myTunnelMosaic1.loadconfig();

	//myTunnelMosaic2.setconfigpath(cfgpath);
	//myTunnelMosaic2.settemppath("E:\\Temp\\");
	//myTunnelMosaic2.setoutputpath("E:\\Temp\\");
	//myTunnelMosaic2.loadconfig();

	//加载上次运行的位置，待完成

	//char bufname[256];
	//std::vector<std::string> names;
	//
	//names.push_back("E:\\Temp\\files\\1\\1.jpg");
	//names.push_back("E:\\Temp\\files\\1\\2.jpg");
	//names.push_back("E:\\Temp\\files\\1\\3.jpg");
	//names.push_back("E:\\Temp\\files\\1\\4.jpg");
	//names.push_back("E:\\Temp\\files\\1\\5.jpg");
	//names.push_back("E:\\Temp\\files\\1\\6.jpg");
	//myTunnelMosaic1.my_stitching(names, 1);


	//names.clear();
	//names.push_back("E:\\Temp\\files\\2\\1.jpg");
	//names.push_back("E:\\Temp\\files\\2\\2.jpg");
	//names.push_back("E:\\Temp\\files\\2\\3.jpg");
	//names.push_back("E:\\Temp\\files\\2\\4.jpg");
	//names.push_back("E:\\Temp\\files\\2\\5.jpg");
	//names.push_back("E:\\Temp\\files\\2\\6.jpg");
	//myTunnelMosaic2.runRadial2(names, 2);


	////保存当前的位置，待完成	
	//myTunnelMosaic2.runAxial(myTunnelMosaic1);







}

