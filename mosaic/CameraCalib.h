#pragma once
#include <iostream>
#include"mosaic.h"

class MOSAIC_API CameraCalib
{
public:
	static int initcalibfile(std::string calibconfigpath, std::string imagelistpath,int outid=0);
	static int calibmain(std::string name, int winsize = 11, float dsize = -1);
};



