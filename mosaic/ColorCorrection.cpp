#include <stdafx.h>
#include <ColorCorrection.h>

// get over mean yuv (GOMY)
//GOMY get_overlap_mean_yuv(Mat img1, Mat img2, int shift_x, int shift_y)
//{
//	GOMY GOMY;
//	GOMY.w = img1.cols - shift_x;
//	GOMY.h = img1.rows - shift_y;
//
//	Mat	overlap_img1;
//	Mat overlap_img2;
//	Mat	overlap_img1_yuv;
//	Mat overlap_img2_yuv;
//
//	//vector<Mat> channelrgb1;
//	//vector<Mat> channelrgb2;
//
//	float	Ysum1 = 0;
//	float	Ysum2 = 0;
//	float	Cbsum1 = 0;
//	float	Cbsum2 = 0;
//	float	Crsum1 = 0;
//	float	Crsum2 = 0;
//	int	count = 0;
//
//	//
//	if (shift_y >= 0)
//	{
//		overlap_img1 = img1(Rect(shift_x, 0, GOMY.w, GOMY.h));
//		//overlap_im1 = img1(1:overlap_h, (shift_x + 1) : size(im1, 2), : );
//
//		overlap_img2 = img2(Rect(0, 0, GOMY.w, GOMY.h));
//		//overlap_im2 = im2((size(im2, 1) - overlap_h + 1) : size(im2, 1), 1 : overlap_w, : );
//		////////////////////////////////////////////////
//		//cout << (int)img2.at<Vec3b>(0, 0)[2] << endl;
//		////////////////////////////////////////////////
//
//
//	}
//	else
//	{
//		overlap_img1 = img1(Rect(GOMY.w - 1, 0, shift_x, GOMY.h));
//		//overlap_im1 = im1(1:overlap_h, (shift_x + 1) : size(im1, 2), : );
//
//		overlap_img2 = img2(Rect(0, 0, GOMY.w, GOMY.h));
//		//overlap_im2 = im2((size(im2, 1) - overlap_h + 1) : size(im2, 1), 1 : overlap_w, : );
//	}
//
//
//	//split(overlap_img1, channelrgb1);
//	//split(overlap_img1, channelrgb2);
//	//YUV1 = RGB2YUVZ(overlap_img1);
//	//YUV2 = RGB2YUVZ(overlap_img2);
//	//cvtColor(overlap_img1, overlap_img1_yuv, COLOR_RGB2YUV);
//	//cvtColor(overlap_img2, overlap_img2_yuv, COLOR_RGB2YUV);
//
//	overlap_img1_yuv = RGB2YUV0(overlap_img1);
//	overlap_img2_yuv = RGB2YUV0(overlap_img2);
//
//	//////////////////////////////////////////////////
//	//cout << overlap_img2_yuv.at<Vec3f>(0,0)[0] << endl;
//	//////////////////////////////////////////////////
//
//	//overlap_img1_yuv.convertTo(overlap_img1_yuv, CV_32FC3);
//	//overlap_img2_yuv.convertTo(overlap_img2_yuv, CV_32FC3);
//
//	// R-> Y ,B-> Cb ,G-> Cr 
//
//	for (int i = 0; i < GOMY.h; i++)
//	{
//		for (int j = 0; j < GOMY.w; j++)
//		{
//			count = count + 1;
//			Ysum1 = Ysum1 + overlap_img1_yuv.at<Vec3f>(i, j)[0] * overlap_img1_yuv.at<Vec3f>(i, j)[0];
//			Ysum2 = Ysum2 + overlap_img2_yuv.at<Vec3f>(i, j)[0] * overlap_img2_yuv.at<Vec3f>(i, j)[0];
//
//
//
//			//Ysum1 = Ysum1 +  YCBCR1 (i,j,1)^2;%^2.2;
//			//Ysum2 = Ysum2 + YCBCR2(i, j, 1) ^ 2; %^ 2.2;
//
//			Cbsum1 = Cbsum1 + overlap_img1_yuv.at<Vec3f>(i, j)[1];
//			Cbsum2 = Cbsum2 + overlap_img2_yuv.at<Vec3f>(i, j)[1];
//			//Cbsum1 = Cbsum1 + YCBCR1(i, j, 2);
//			//Cbsum2 = Cbsum2 + YCBCR2(i, j, 2);
//
//			Crsum1 = Crsum1 + overlap_img1_yuv.at<Vec3f>(i, j)[2];
//			Crsum2 = Crsum2 + overlap_img2_yuv.at<Vec3f>(i, j)[2];
//			//Crsum1 = Crsum1 + YCBCR1(i, j, 3);
//			//Crsum2 = Crsum2 + YCBCR2(i, j, 3);
//
//		}
//	}
//
//	GOMY.Bi.push_back(log(Ysum1 / count));
//	GOMY.Bi.push_back(log(Ysum2 / count));
//
//	GOMY.Sicb.push_back(Cbsum1 / count);
//	GOMY.Sicb.push_back(Cbsum2 / count);
//
//	GOMY.Sicr.push_back(Crsum1 / count);
//	GOMY.Sicr.push_back(Crsum2 / count);
//
//	return GOMY;
//}

//已完成
//Mat RGB2YUV0(Mat src)
//{
//	src.convertTo(src, CV_32FC3);
//
//	Mat dst = src.clone();
//
//	int c = src.cols;
//	int r = src.rows;
//
//	for (int i = 0; i < r; i++)
//	{
//		for (int j = 0; j < c; j++)
//		{
//			dst.at<Vec3f>(i, j)[0] = floor(0.299 * src.at<Vec3f>(i, j)[2] + 0.587 * src.at<Vec3f>(i, j)[1] + 0.114 * src.at<Vec3f>(i, j)[0]);
//			dst.at<Vec3f>(i, j)[1] = floor(-0.14713 * src.at<Vec3f>(i, j)[2] - 0.28886 * src.at<Vec3f>(i, j)[1] + 0.436 * src.at<Vec3f>(i, j)[0]);
//			dst.at<Vec3f>(i, j)[2] = floor(0.615 * src.at<Vec3f>(i, j)[2] - 0.51499 * src.at<Vec3f>(i, j)[1] - 0.10001 * src.at<Vec3f>(i, j)[0]);
//
//		}
//	}
//	//y = floor(0.299 * r + 0.587 * g + 0.114 * b);
//	//u = floor(-0.14713 * r - 0.28886 * g + 0.436 * b);
//	//v = floor(0.615 * r - 0.51499 * g - 0.10001 * b);
//
//	return dst;
//}

////已完成
//Mat YUV2RGB0(Mat src)
//{
//	src.convertTo(src, CV_32FC3);
//	//cout << src.at<Vec3f>(0, 0)[2] << endl;
//	Mat dst(src.rows, src.cols, CV_32FC3);
//
//	for (int i = 0; i < src.rows; i++)
//	{
//		for (int j = 0; j < src.cols; j++)
//		{
//			dst.at<Vec3f>(i, j)[2] = src.at<Vec3f>(i, j)[0] + 1.13983 * src.at<Vec3f>(i, j)[2];
//			//if (dst.at<Vec3f>(i, j)[2] > 255)dst.at<Vec3f>(i, j)[2] = 255;
//			//if (dst.at<Vec3f>(i, j)[2] < 0)dst.at<Vec3f>(i, j)[2] = 0;
//			dst.at<Vec3f>(i, j)[1] = src.at<Vec3f>(i, j)[0] - 0.39465 * src.at<Vec3f>(i, j)[1] - 0.5806 * src.at<Vec3f>(i, j)[2];
//			//if (dst.at<Vec3f>(i, j)[1] > 255)dst.at<Vec3f>(i, j)[1] = 255;
//			//if (dst.at<Vec3f>(i, j)[1] < 0)dst.at<Vec3f>(i, j)[1] = 0;
//			dst.at<Vec3f>(i, j)[0] = src.at<Vec3f>(i, j)[0] + 2.03211 * src.at<Vec3f>(i, j)[1];
//			//if (dst.at<Vec3f>(i, j)[0] > 255)dst.at<Vec3f>(i, j)[0] = 255;
//			//if (dst.at<Vec3f>(i, j)[0] < 0)dst.at<Vec3f>(i, j)[0] = 0;
//		}
//	}
//
//	//cout << dst.at<Vec3f>(0, 0)[2] << endl;
//
//	//dst.convertTo(dst, CV_8UC3);
//
//	//r = y + 1.13983 * v;
//	//g = y - 0.39465 * u - 0.5806 * v;
//	//b = y + 2.03211 * u;
//
//
//	return dst;
//}

////  gamma ֵ 
//vector<float> getCorrParameters(vector<vector<float> > B)
//{
//	float sigma_n = 2; //  /255;
//	float sigma_g = 0.5; //  /255;
//
//	vector<float> gamma0;
//
//	vector<float> a = { 0 };
//	vector<float> b = { (B[0][0] * B[0][0]) / (sigma_n * sigma_n) + 1 / (sigma_g * sigma_g) };
//	vector<float> c = { (-B[0][0] * B[0][1]) / (sigma_n * sigma_n) };
//	vector<float> d = { 1 / (sigma_g * sigma_g) };
//
//	//a.push_back(0);
//
//	vector<float> bn;
//	vector<float> an;
//	vector<float> dn;
//
//	bn.push_back(b[0]);
//	an.push_back(a[0]);
//	dn.push_back(d[0]);
//
//
//	// 初值 i = 2 ?
//	for (int i = 1; i < num_img; i++)
//	{
//		a.push_back(-B[i - 1][0] * B[i - 1][1] / (sigma_n * sigma_n));
//		if (i == num_img - 1)
//		{
//			b.push_back((B[i - 1][1] * B[i - 1][1]) / (sigma_n * sigma_n) + 1 / (sigma_g * sigma_g));
//			c.push_back(0);
//		}
//		else
//		{
//			b.push_back((B[i - 1][1] * B[i - 1][1] + B[i][0] * B[i][0]) / (sigma_n * sigma_n) + 1 / (sigma_g * sigma_g));
//			c.push_back(-B[i][0] * B[i][1] / (sigma_n * sigma_n));
//		}
//
//		d.push_back(1 / (sigma_g * sigma_g));
//
//		an.push_back(a[i] / bn[i - 1]);
//		bn.push_back(b[i] - a[i] / bn[i - 1] * c[i - 1]);
//		dn.push_back(d[i] - an[i] * dn[i - 1]);
//	}
//	gamma0.resize(num_img);
//	gamma0[num_img - 1] = dn[num_img - 1] / bn[num_img - 1];
//
//	for (int i = (num_img - 2); i >= 0; i--)
//	{
//		gamma0[i] = (dn[i] - c[i] * gamma0[i + 1]) / bn[i];
//	}
//
//	return gamma0;
//}

////行值
//minhang minhh(vector<vector<float> > a, int m)
//{
//	minhang b;
//	int w = a[0].size();
//	int h = a.size();
//	vector<float> temp;
//
//	if (m < h)
//	{
//		for (int i = 0; i < w; i++)
//		{
//			temp.push_back(a[m][i]);
//		}
//
//		float min = temp[1];
//
//		for (int i = 1; i < w; i++)
//		{
//			if (min > temp[i])
//			{
//				min = temp[i];
//				b.mindiff = temp[i];
//				b.seamindex = i;
//			}
//		}
//
//
//	}
//
//	//b.mindiff = *min_element(temp.begin(), temp.end());
//
//	return b;
//}

////color correction 
////559行bug 边缘像素0导致跳不出判断 极限位置  注意此函数输入图像的边界 
////修起来太麻烦 不整了
//void Color_correction(Mat *im1p,Mat *im2p)
//{
//	//std::vector<std::string> imgpathlist;
//	////std::vector<std::string> num;
//
//	//imgpathlist.push_back("E:\\Temp\\files\\1\\1.jpg");
//	//imgpathlist.push_back("E:\\Temp\\files\\1\\2.jpg");
//	//imgpathlist.push_back("E:\\Temp\\files\\1\\3.jpg");
//	//imgpathlist.push_back("E:\\Temp\\files\\1\\4.jpg");
//	//imgpathlist.push_back("E:\\Temp\\files\\1\\5.jpg");
//	//imgpathlist.push_back("E:\\Temp\\files\\1\\6.jpg");
//
//	Mat im1 = *im1p;
//	Mat im2 = *im2p;
//
//	int W = im1.cols;
//	int H = im1.rows;
//
//	//im1.convertTo(im1, CV_32FC3);
//	int	shift_x = W - /*更改值*/200; // 暂定重合宽度200  //偏移调整0测试
//	int	shift_y = 0;
//
//	GOMY GOMY;
//	vector<vector<float> >	B;
//	vector<vector<float> > Scb;
//	vector<vector<float> > Scr;
//	int	overlap_w;
//	int overlap_h;
//	//overlap_w.resize(im1.cols);
//	//overlap_h.resize(im1.rows);
//
//	vector<float> gamma;
//	vector<float> acb;
//	vector<float> acr;
//
//	string js[num_img];
//
//
//
//	for (int i = 0; i < num_img; i++)
//	{
//		js[i] = i + 1;
//	}
//
//	B.resize(num_img - 1);
//	Scb.resize(num_img - 1);
//	Scr.resize(num_img - 1);
//
//
//	for (int ii = 1; ii < num_img; ii++)
//	{
//		//im2 = cv::imread(imgpathlist[ii]);
//		//im2.convertTo(im2, CV_32FC3);
//		//cout << (int)im2.at<Vec3b>(0, 0)[2] << endl;
//		GOMY = get_overlap_mean_yuv(im1, im2, shift_x, shift_y);
//
//		if (ii == 1)
//		{
//			int rrr = GOMY.Bi.size();
//			for (int k = 0; k < num_img - 1; k++)
//			{
//				B[k].resize(rrr);
//				Scb[k].resize(rrr);
//				Scr[k].resize(rrr);
//			}
//		}
//
//		for (int j = 0; j < 2; j++)
//		{
//			B[ii - 1][j] = GOMY.Bi[j];
//			Scb[ii - 1][j] = GOMY.Sicb[j];
//			Scr[ii - 1][j] = GOMY.Sicr[j];
//
//		}
//
//		overlap_w = GOMY.w;
//		overlap_h = GOMY.h;
//
//		//im1 = im2.clone();
//
//	}
//
//
//
//	gamma = getCorrParameters(B);
//	acb = getCorrParameters(Scb);
//	acr = getCorrParameters(Scr);
//	//////////////////////////////////////////////////
//	//cout << gamma[0] << endl;
//	//cout << gamma[1] << endl;
//	//cout << gamma[2] << endl;
//	//cout << gamma[3] << endl;
//	//cout << gamma[4] << endl;
//	//cout << gamma[5] << endl;
//	//////////////////////////////////////////////////
//
//
//	//色彩和亮度调整
//	//im1 = cv::imread(imgpathlist[0]);
//
//	Mat im1YUV;
//	Mat im2YUV;
//	//cv::cvtColor(im1, im1YUV, COLOR_BGR2YUV);
//	im1YUV = RGB2YUV0(im1);
//
//
//	im1YUV.convertTo(im1YUV, CV_32FC3);
//
//	for (int i = 0; i < im1.rows; i++)
//	{
//		for (int j = 0; j < im1.cols; j++)
//		{
//			im1YUV.at<Vec3f>(i, j)[0] = pow(im1YUV.at<Vec3f>(i, j)[0], gamma[0]);
//			im1YUV.at<Vec3f>(i, j)[1] = im1YUV.at<Vec3f>(i, j)[1] * acb[0];
//			im1YUV.at<Vec3f>(i, j)[2] = im1YUV.at<Vec3f>(i, j)[2] * acr[0];
//		}
//	}
//
//	im1 = YUV2RGB0(im1YUV);
//	//////////////////////////////////////////////////
//	//cout << im1.at<Vec3b>(0, 0)[2] << endl;
//	//////////////////////////////////////////////////
//	//im1YUV.convertTo(im1YUV,CV_8UC3); //???
//	//cv::cvtColor(im1YUV, im1, COLOR_YUV2BGR);
//
//	//
//	//im1 = YUV2YUVMat(YCBCR);
//	//im1.convertTo(im1, CV_32FC3);
//	//RGB rgb_im1 = YUV2RGB(im1);
//	//im1 = RGB2RGBMat(rgb_im1);
//
//	//cout << imgpathlist.size() << endl;
//
//	for (int ii = 1; ii < num_img; ii++)
//	{
//		//im2 = imread(imgpathlist[ii]);
//		//im2.convertTo(im2, CV_32FC3);
//		//cout << (int)im2.at<Vec3b>(0, 0)[2] << endl;
//		//cout << (int)im2.at<Vec3b>(0, 0)[1] << endl;
//		//cout << (int)im2.at<Vec3b>(0, 0)[0] << endl;
//
//
//
//		im2YUV = RGB2YUV0(im2);
//		//cout << im2YUV.at<Vec3f>(0, 0)[0] << endl;
//		//cout << im2YUV.at<Vec3f>(0, 0)[1] << endl;
//		//cout << im2YUV.at<Vec3f>(0, 0)[2] << endl;
//
//
//		//cv::cvtColor(im2, im2YUV, COLOR_BGR2YUV);
//		//im2YUV.convertTo(im2YUV, CV_32FC3);
//
//		//cout << (int)im2YUV.at<Vec3b>(0, 0)[0] << endl;
//
//		for (int i = 0; i < im2.rows; i++)
//		{
//			for (int j = 0; j < im2.cols; j++)
//			{
//				im2YUV.at<Vec3f>(i, j)[0] = pow(im2YUV.at<Vec3f>(i, j)[0], gamma[ii]);
//				//cout << (int)im2YUV.at<Vec3b>(0, 0)[0] << endl;
//				im2YUV.at<Vec3f>(i, j)[1] = im2YUV.at<Vec3f>(i, j)[1] * acb[ii];
//				im2YUV.at<Vec3f>(i, j)[2] = im2YUV.at<Vec3f>(i, j)[2] * acr[ii];
//			}
//		}
//		//cout << im2YUV.at<Vec3f>(0, 0)[0] << endl;
//		//cout << im2YUV.at<Vec3f>(0, 0)[1] << endl;
//		//cout << im2YUV.at<Vec3f>(0, 0)[2] << endl;
//
//
//
//		//im2YUV.convertTo(im2YUV, CV_8UC3);
//		//cv::cvtColor(im2YUV, im2, COLOR_YUV2BGR);
//
//		//cout << (int)im2YUV.at<Vec3b>(0, 0)[0] << endl;
//
//		im2 = YUV2RGB0(im2YUV);
//
//		//im1.convertTo(im1, CV_32FC3);
//		//im2.convertTo(im2, CV_32FC3);
//
//		////////////////////RGB YUV 数值不匹配？？   im2    
//		//cout << (int)im2.at<Vec3b>(0, 0)[2] << endl;
//		//cout << (int)im2.at<Vec3b>(0, 0)[1] << endl;
//		//cout << (int)im2.at<Vec3b>(0, 0)[0] << endl;
//
//		//cout << im2.empty() << endl;
//
//		//seam出溢出 -200
//		//shift_x = W - 200;
//		//shift_y = 0;
//
//		overlap_w = im1.cols - shift_x;
//		overlap_h = im1.rows - shift_y;
//
//		Mat	overlap_img1;
//		Mat overlap_img2;
//
//		//cout << overlap_img1.rows << endl;  //行数√
//		//cout << overlap_img1.at<Vec3b>(202, 0)[0] << endl;
//
//
//		if (shift_y >= 0)
//		{
//			overlap_img1 = im1(Rect(shift_x, 0, overlap_w, overlap_h))/*.clone()*/;
//			//overlap_im1 = img1(1:overlap_h, (shift_x + 1) : size(im1, 2), : );
//
//
//			overlap_img2 = im2(Rect(0, 0, overlap_w, overlap_h))/*.clone()*/;
//			//overlap_im2 = im2((size(im2, 1) - overlap_h + 1) : size(im2, 1), 1 : overlap_w, : );
//
//		}
//		else
//		{
//			overlap_img1 = im1(Rect(shift_x, 0, GOMY.w, GOMY.h)).clone();
//			//overlap_im1 = im1(1:overlap_h, (shift_x + 1) : size(im1, 2), : );
//
//			overlap_img2 = im2(Rect(0, 0, GOMY.w, GOMY.h)).clone();
//			//overlap_im2 = im2((size(im2, 1) - overlap_h + 1) : size(im2, 1), 1 : overlap_w, : );
//		}
//
//		//cout << (int)overlap_img2.at<Vec3b>(0, 0)[2] << endl;
//		//cout << (int)overlap_img2.at<Vec3b>(0, 0)[1] << endl;
//		//cout << (int)overlap_img2.at<Vec3b>(0, 0)[0] << endl;
//		//cout << overlap_img1.at<Vec3b>(0, 0)[0] <<
//		//	overlap_img1.size() << endl;
//
//		//overlap_img1.convertTo(overlap_img1, CV_32FC3);
//		//overlap_img2.convertTo(overlap_img2, CV_32FC3);
//
//
//
//		vector<vector<float> > squared_diff(overlap_h+1, vector<float>(overlap_w+1, 0));
//		vector<vector<float> > accumulate_D(overlap_h+1, vector<float>(overlap_w+1, 0));
//		vector<vector<float> > color_diff(overlap_h+1, vector<float>(3, 0));
//		vector<float> neighbor;
//
//		for (int i = 1; i <= overlap_h; i++)
//		{
//			for (int j = 1; j <= overlap_w; j++)
//			{
//				squared_diff[i][j] = pow(/*(int)*/overlap_img1.at<Vec3f>(i - 1, j - 1)[2] - /*(int)*/overlap_img2.at<Vec3f>(i - 1, j - 1)[2], 2);
//
//				//cout << (int)overlap_img1.at<Vec3b>(0, 0)[2] << endl;
//				//cout << (int)overlap_img2.at<Vec3b>(0, 0)[2] << endl;
//
//				if (i == 1)
//				{
//					accumulate_D[i][j] = squared_diff[i][j];
//
//				}
//				else
//				{
//					if (j == 1)
//					{
//						neighbor.push_back(accumulate_D[i - 1][j]);
//						neighbor.push_back(accumulate_D[i - 1][j + 1]);
//					}
//					else
//					{
//						if (j == overlap_w)
//						{
//							neighbor.push_back(accumulate_D[i - 1][j - 1]);
//							neighbor.push_back(accumulate_D[i - 1][j]);
//						}
//						else
//						{
//							neighbor.push_back(accumulate_D[i - 1][j - 1]);
//							neighbor.push_back(accumulate_D[i - 1][j]);
//							neighbor.push_back(accumulate_D[i - 1][j + 1]);
//						}
//					}
//					accumulate_D[i][j] = squared_diff[i][j] + *min_element(neighbor.begin(), neighbor.end());
//					neighbor.clear();
//				}
//			}
//		}
//
//
//		//for (int j = 0; j < overlap_w; j++)
//		//{
//		//	/*cout << squared_diff[i][j] << endl;*/
//		//	cout << accumulate_D[overlap_h - 1][j] << endl;
//		//}
//
//		minhang mindst;
//		float mindiff;
//		vector<int> seamindex;
//		seamindex.resize(overlap_h + 1, 0);
//
//		mindst = minhh(accumulate_D, overlap_h);
//		mindiff = mindst.mindiff;
//		seamindex[overlap_h] = mindst.seamindex;
//
//		//cout << mindst.mindiff << endl;
//		//cout << mindst.seamindex << endl;
//
//		////有bug   //边缘取值达到0 
//		for (int i = overlap_h - 1; i >= 1; i--)
//		{
//			if (seamindex[i + 1] == 1)
//			{
//				if ((accumulate_D[i][seamindex[i + 1]] + squared_diff[i + 1][seamindex[i + 1]]) == mindiff)
//				{
//					seamindex[i] = seamindex[i + 1];
//					mindiff = accumulate_D[i][seamindex[i]];
//				}
//				else
//				{
//					if ((accumulate_D[i][seamindex[i + 1] + 1] + squared_diff[i + 1][seamindex[i + 1]]) == mindiff)
//					{
//						seamindex[i] = seamindex[i + 1] + 1;
//						mindiff = accumulate_D[i][seamindex[i]];
//					}
//				}
//			}
//			else
//			{
//				if (seamindex[i + 1] == overlap_w)
//				{
//					if ((accumulate_D[i][seamindex[i + 1] - 1] + squared_diff[i + 1][seamindex[i + 1]]) == mindiff)
//					{
//						seamindex[i] = seamindex[i + 1] - 1;
//						mindiff = accumulate_D[i][seamindex[i]];
//					}
//					else
//					{
//						if ((accumulate_D[i][seamindex[i + 1]] + squared_diff[i + 1][seamindex[i + 1]]) == mindiff)
//						{
//							seamindex[i] = seamindex[i + 1];
//							mindiff = accumulate_D[i][seamindex[i]];
//						}
//					}
//				}
//				else
//				{
//					if (accumulate_D[i][seamindex[i + 1] - 1] + squared_diff[i + 1][seamindex[i + 1]] == mindiff)
//					{
//						seamindex[i] = seamindex[i + 1] - 1;
//						mindiff = accumulate_D[i][seamindex[i]];
//					}
//					else
//					{
//						if ((accumulate_D[i][seamindex[i + 1]] + squared_diff[i + 1][seamindex[i + 1]]) == mindiff)
//						{
//							seamindex[i] = seamindex[i + 1];
//							mindiff = accumulate_D[i][seamindex[i]];
//						}
//						else
//						{
//							if ((accumulate_D[i][seamindex[i + 1] + 1] + squared_diff[i + 1][seamindex[i + 1]]) == mindiff)
//							{
//								seamindex[i] = seamindex[i + 1] + 1;
//								mindiff = accumulate_D[i][seamindex[i]];
//							}
//						}
//					}
//				}
//				
//			}
//		}
//
//
//		Mat temp(overlap_h, overlap_w, CV_32FC3);
//
//		for (int i = 0; i < overlap_h; i++)
//		{
//			for (int j = 0; j < 3; j++)
//			{
//				for (int k = 0; k < seamindex[i+1] - 1; k++)
//				{
//					temp.at<Vec3f>(i, k)[j] = overlap_img1.at<Vec3f>(i, k)[j];
//				}
//
//				//cout << "seamindex:" << seamindex[i] << endl;
//
//				temp.at<Vec3f>(i, seamindex[i+1])[j] = 0;
//
//				for (int k = seamindex[i+1] + 1; k < overlap_w; k++)
//				{
//					temp.at<Vec3f>(i, k)[j] = overlap_img2.at<Vec3f>(i, k)[j];
//				}
//			}
//
//			for (int j = 0; j < 3; j++)
//			{
//				color_diff[i][j] = /*(int)*/overlap_img1.at<Vec3f>(i, seamindex[i+1])[j] - /*(int)*/overlap_img2.at<Vec3f>(i, seamindex[i])[j];
//				//cout << "color_diff:" << color_diff[i][j] << endl;
//
//
//			}
//		}
//
//		//double total = 0;
//		//double blendvalue[3][1] = { 0 };
//
//		vector<float> dist;
//		dist.resize(overlap_h);
//
//		for (int i = 0; i < overlap_h; i++)
//		{
//			for (int j = seamindex[i+1] + 1; j < im2.cols; j++)
//			{
//				float total = 0;
//				float blendvalue[3] = { 0 };
//
//				for (int k = 0; k < overlap_h; k++)
//				{
//					dist[k] = abs(i - k) + abs(j - (seamindex[k+1] /*- 1*/));
//					//cout << "dist[k]:" << dist[k] << endl;
//					//cout << "color_diff:" << color_diff[k][0] << endl;
//
//					blendvalue[0] = blendvalue[0] + 1 / dist[k] * color_diff[k][0];
//					blendvalue[1] = blendvalue[1] + 1 / dist[k] * color_diff[k][1];
//					blendvalue[2] = blendvalue[2] + 1 / dist[k] * color_diff[k][2];
//
//					total = total + 1 / dist[k];
//
//
//				}
//
//				//cout << "total:" << total << endl;
//				//cout << (int)im2.at<Vec3b>(i, j)[0] << endl;
//				//cout <<"blendvalue:"<< blendvalue[2] << endl;
//				im2.at<Vec3f>(i, j)[0] = /*(int)*/im2.at<Vec3f>(i, j)[0] + /*(int)*/(blendvalue[2] / total);
//				im2.at<Vec3f>(i, j)[1] = /*(int)*/im2.at<Vec3f>(i, j)[1] + /*(int)*/(blendvalue[1] / total);
//				im2.at<Vec3f>(i, j)[2] = /*(int)*/im2.at<Vec3f>(i, j)[2] + /*(int)*/(blendvalue[0] / total);
//
//				//RGB 2 1 0
//
//			}
//		}
//
//		im1.convertTo(im1, CV_8UC3);
//		im2.convertTo(im2, CV_8UC3);
//
//		*im1p = im1.clone();
//		*im2p = im2.clone();
//
//		//if (ii == 1)cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\1.jpg", im1);
//		//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\2.jpg", im2);
//
//	}
//}



GOMY get_overlap_mean_yuv(Mat img1, Mat img2, int shift_x, int shift_y)
{
	GOMY G;
	G.w = img1.cols - shift_x;
	G.h = img1.rows - shift_y;

	Mat	overlap_img1;
	Mat overlap_img2;
	Mat	overlap_img1_yuv;
	Mat overlap_img2_yuv;

	//vector<Mat> channelrgb1;
	//vector<Mat> channelrgb2;

	float	Ysum1 = 0;
	float	Ysum2 = 0;
	float	Cbsum1 = 0;
	float	Cbsum2 = 0;
	float	Crsum1 = 0;
	float	Crsum2 = 0;
	int	count = 0;

	//
	if (shift_y >= 0)
	{
		overlap_img1 = img1(Rect(shift_x, 0, G.w, G.h));
		//overlap_im1 = img1(1:overlap_h, (shift_x + 1) : size(im1, 2), : );

		overlap_img2 = img2(Rect(0, 0, G.w, G.h));
		//overlap_im2 = im2((size(im2, 1) - overlap_h + 1) : size(im2, 1), 1 : overlap_w, : );
		////////////////////////////////////////////////
		//cout << (int)img2.at<Vec3b>(0, 0)[2] << endl;
		////////////////////////////////////////////////


	}
	else
	{
		overlap_img1 = img1(Rect(G.w - 1, 0, shift_x, G.h));
		//overlap_im1 = im1(1:overlap_h, (shift_x + 1) : size(im1, 2), : );

		overlap_img2 = img2(Rect(0, 0, G.w, G.h));
		//overlap_im2 = im2((size(im2, 1) - overlap_h + 1) : size(im2, 1), 1 : overlap_w, : );
	}


	//split(overlap_img1, channelrgb1);
	//split(overlap_img1, channelrgb2);
	//YUV1 = RGB2YUVZ(overlap_img1);
	//YUV2 = RGB2YUVZ(overlap_img2);
	//cvtColor(overlap_img1, overlap_img1_yuv, COLOR_RGB2YUV);
	//cvtColor(overlap_img2, overlap_img2_yuv, COLOR_RGB2YUV);

	overlap_img1_yuv = RGB2YUV0(overlap_img1);
	overlap_img2_yuv = RGB2YUV0(overlap_img2);

	//////////////////////////////////////////////////
	//cout << overlap_img2_yuv.at<Vec3f>(0,0)[0] << endl;
	//////////////////////////////////////////////////

	//overlap_img1_yuv.convertTo(overlap_img1_yuv, CV_32FC3);
	//overlap_img2_yuv.convertTo(overlap_img2_yuv, CV_32FC3);

	// R-> Y ,B-> Cb ,G-> Cr 

	for (int i = 0; i < G.h; i++)
	{
		for (int j = 0; j < G.w; j++)
		{
			count = count + 1;
			Ysum1 = Ysum1 + overlap_img1_yuv.at<Vec3f>(i, j)[0] * overlap_img1_yuv.at<Vec3f>(i, j)[0];
			Ysum2 = Ysum2 + overlap_img2_yuv.at<Vec3f>(i, j)[0] * overlap_img2_yuv.at<Vec3f>(i, j)[0];



			//Ysum1 = Ysum1 +  YCBCR1 (i,j,1)^2;%^2.2;
			//Ysum2 = Ysum2 + YCBCR2(i, j, 1) ^ 2; %^ 2.2;

			Cbsum1 = Cbsum1 + overlap_img1_yuv.at<Vec3f>(i, j)[1];
			Cbsum2 = Cbsum2 + overlap_img2_yuv.at<Vec3f>(i, j)[1];
			//Cbsum1 = Cbsum1 + YCBCR1(i, j, 2);
			//Cbsum2 = Cbsum2 + YCBCR2(i, j, 2);

			Crsum1 = Crsum1 + overlap_img1_yuv.at<Vec3f>(i, j)[2];
			Crsum2 = Crsum2 + overlap_img2_yuv.at<Vec3f>(i, j)[2];
			//Crsum1 = Crsum1 + YCBCR1(i, j, 3);
			//Crsum2 = Crsum2 + YCBCR2(i, j, 3);

		}
	}

	G.Bi.push_back(log(Ysum1 / count));
	G.Bi.push_back(log(Ysum2 / count));

	G.Sicb.push_back(Cbsum1 / count);
	G.Sicb.push_back(Cbsum2 / count);

	G.Sicr.push_back(Crsum1 / count);
	G.Sicr.push_back(Crsum2 / count);

	return G;
}

Mat RGB2YUV0(Mat src)
{
	src.convertTo(src, CV_32FC3);

	Mat dst = src.clone();

	int c = src.cols;
	int r = src.rows;

	for (int i = 0; i < r; i++)
	{
		for (int j = 0; j < c; j++)
		{
			dst.at<Vec3f>(i, j)[0] = floor(0.299 * src.at<Vec3f>(i, j)[2] + 0.587 * src.at<Vec3f>(i, j)[1] + 0.114 * src.at<Vec3f>(i, j)[0]);
			dst.at<Vec3f>(i, j)[1] = floor(-0.14713 * src.at<Vec3f>(i, j)[2] - 0.28886 * src.at<Vec3f>(i, j)[1] + 0.436 * src.at<Vec3f>(i, j)[0]);
			dst.at<Vec3f>(i, j)[2] = floor(0.615 * src.at<Vec3f>(i, j)[2] - 0.51499 * src.at<Vec3f>(i, j)[1] - 0.10001 * src.at<Vec3f>(i, j)[0]);

		}
	}
	//y = floor(0.299 * r + 0.587 * g + 0.114 * b);
	//u = floor(-0.14713 * r - 0.28886 * g + 0.436 * b);
	//v = floor(0.615 * r - 0.51499 * g - 0.10001 * b);

	return dst;
}

Mat YUV2RGB0(Mat src)
{
	src.convertTo(src, CV_32FC3);
	//cout << src.at<Vec3f>(0, 0)[2] << endl;
	Mat dst(src.rows, src.cols, CV_32FC3);

	for (int i = 0; i < src.rows; i++)
	{
		for (int j = 0; j < src.cols; j++)
		{
			dst.at<Vec3f>(i, j)[2] = src.at<Vec3f>(i, j)[0] + 1.13983 * src.at<Vec3f>(i, j)[2];
			//if (dst.at<Vec3f>(i, j)[2] > 255)dst.at<Vec3f>(i, j)[2] = 255;
			//if (dst.at<Vec3f>(i, j)[2] < 0)dst.at<Vec3f>(i, j)[2] = 0;
			dst.at<Vec3f>(i, j)[1] = src.at<Vec3f>(i, j)[0] - 0.39465 * src.at<Vec3f>(i, j)[1] - 0.5806 * src.at<Vec3f>(i, j)[2];
			//if (dst.at<Vec3f>(i, j)[1] > 255)dst.at<Vec3f>(i, j)[1] = 255;
			//if (dst.at<Vec3f>(i, j)[1] < 0)dst.at<Vec3f>(i, j)[1] = 0;
			dst.at<Vec3f>(i, j)[0] = src.at<Vec3f>(i, j)[0] + 2.03211 * src.at<Vec3f>(i, j)[1];
			//if (dst.at<Vec3f>(i, j)[0] > 255)dst.at<Vec3f>(i, j)[0] = 255;
			//if (dst.at<Vec3f>(i, j)[0] < 0)dst.at<Vec3f>(i, j)[0] = 0;
		}
	}

	//cout << dst.at<Vec3f>(0, 0)[2] << endl;

	//dst.convertTo(dst, CV_8UC3);

	//r = y + 1.13983 * v;
	//g = y - 0.39465 * u - 0.5806 * v;
	//b = y + 2.03211 * u;


	return dst;
}

vector<float> getCorrParameters(vector<vector<float>> B)
{
	float sigma_n = 2; //  /255;
	float sigma_g = 0.5; //  /255;

	vector<float> gamma0;

	vector<float> a = { 0 };
	vector<float> b = { (B[0][0] * B[0][0]) / (sigma_n * sigma_n) + 1 / (sigma_g * sigma_g) };
	vector<float> c = { (-B[0][0] * B[0][1]) / (sigma_n * sigma_n) };
	vector<float> d = { 1 / (sigma_g * sigma_g) };

	//a.push_back(0);

	vector<float> bn;
	vector<float> an;
	vector<float> dn;

	bn.push_back(b[0]);
	an.push_back(a[0]);
	dn.push_back(d[0]);


	// 初值 i = 2 ?
	for (int i = 1; i < num_img; i++)
	{
		a.push_back(-B[i - 1][0] * B[i - 1][1] / (sigma_n * sigma_n));
		if (i == num_img - 1)
		{
			b.push_back((B[i - 1][1] * B[i - 1][1]) / (sigma_n * sigma_n) + 1 / (sigma_g * sigma_g));
			c.push_back(0);
		}
		else
		{
			b.push_back((B[i - 1][1] * B[i - 1][1] + B[i][0] * B[i][0]) / (sigma_n * sigma_n) + 1 / (sigma_g * sigma_g));
			c.push_back(-B[i][0] * B[i][1] / (sigma_n * sigma_n));
		}

		d.push_back(1 / (sigma_g * sigma_g));

		an.push_back(a[i] / bn[i - 1]);
		bn.push_back(b[i] - a[i] / bn[i - 1] * c[i - 1]);
		dn.push_back(d[i] - an[i] * dn[i - 1]);
	}
	gamma0.resize(num_img);
	gamma0[num_img - 1] = dn[num_img - 1] / bn[num_img - 1];

	for (int i = (num_img - 2); i >= 0; i--)
	{
		gamma0[i] = (dn[i] - c[i] * gamma0[i + 1]) / bn[i];
	}

	return gamma0;
}

minhang minhh(vector<vector<float>> a, int m)
{
	minhang b;
	int w = a[0].size();
	int h = a.size();
	vector<float> temp;

	if (m < h)
	{
		for (int i = 0; i < w; i++)
		{
			temp.push_back(a[m][i]);
		}

		float min = temp[1];

		for (int i = 1; i < w; i++)
		{
			if (min > temp[i])
			{
				min = temp[i];
				b.mindiff = temp[i];
				b.seamindex = i;
			}
		}
	}

	//b.mindiff = *min_element(temp.begin(), temp.end());

	return b;
}

void Color_correction(Mat* im1p, Mat* im2p)
{
	//std::vector<std::string> imgpathlist;
	////std::vector<std::string> num;

	//imgpathlist.push_back("E:\\Temp\\files\\1\\1.jpg");
	//imgpathlist.push_back("E:\\Temp\\files\\1\\2.jpg");
	//imgpathlist.push_back("E:\\Temp\\files\\1\\3.jpg");
	//imgpathlist.push_back("E:\\Temp\\files\\1\\4.jpg");
	//imgpathlist.push_back("E:\\Temp\\files\\1\\5.jpg");
	//imgpathlist.push_back("E:\\Temp\\files\\1\\6.jpg");

	Mat im1 = *im1p;
	Mat im2 = *im2p;

	int W = im1.cols;
	int H = im1.rows;

	//im1.convertTo(im1, CV_32FC3);
	int	shift_x = W - /*更改值*/200; // 暂定重合宽度200  //偏移调整0测试
	int	shift_y = 0;

	GOMY gomy;
	vector<vector<float> >	B;
	vector<vector<float> > Scb;
	vector<vector<float> > Scr;
	int	overlap_w;
	int overlap_h;
	//overlap_w.resize(im1.cols);
	//overlap_h.resize(im1.rows);

	vector<float> gamma;
	vector<float> acb;
	vector<float> acr;

	string js[num_img];



	for (int i = 0; i < num_img; i++)
	{
		js[i] = i + 1;
	}

	B.resize(num_img - 1);
	Scb.resize(num_img - 1);
	Scr.resize(num_img - 1);


	for (int ii = 1; ii < num_img; ii++)
	{
		//im2 = cv::imread(imgpathlist[ii]);
		//im2.convertTo(im2, CV_32FC3);
		//cout << (int)im2.at<Vec3b>(0, 0)[2] << endl;
		gomy = get_overlap_mean_yuv(im1, im2, shift_x, shift_y);

		if (ii == 1)
		{
			int rrr = gomy.Bi.size();
			for (int k = 0; k < num_img - 1; k++)
			{
				B[k].resize(rrr);
				Scb[k].resize(rrr);
				Scr[k].resize(rrr);
			}
		}

		for (int j = 0; j < 2; j++)
		{
			B[ii - 1][j] = gomy.Bi[j];
			Scb[ii - 1][j] = gomy.Sicb[j];
			Scr[ii - 1][j] = gomy.Sicr[j];

		}

		overlap_w = gomy.w;
		overlap_h = gomy.h;

		//im1 = im2.clone();

	}



	gamma = getCorrParameters(B);
	acb = getCorrParameters(Scb);
	acr = getCorrParameters(Scr);
	//////////////////////////////////////////////////
	//cout << gamma[0] << endl;
	//cout << gamma[1] << endl;
	//cout << gamma[2] << endl;
	//cout << gamma[3] << endl;
	//cout << gamma[4] << endl;
	//cout << gamma[5] << endl;
	//////////////////////////////////////////////////


	//色彩和亮度调整
	//im1 = cv::imread(imgpathlist[0]);

	Mat im1YUV;
	Mat im2YUV;
	//cv::cvtColor(im1, im1YUV, COLOR_BGR2YUV);
	im1YUV = RGB2YUV0(im1);


	im1YUV.convertTo(im1YUV, CV_32FC3);

	for (int i = 0; i < im1.rows; i++)
	{
		for (int j = 0; j < im1.cols; j++)
		{
			im1YUV.at<Vec3f>(i, j)[0] = pow(im1YUV.at<Vec3f>(i, j)[0], gamma[0]);
			im1YUV.at<Vec3f>(i, j)[1] = im1YUV.at<Vec3f>(i, j)[1] * acb[0];
			im1YUV.at<Vec3f>(i, j)[2] = im1YUV.at<Vec3f>(i, j)[2] * acr[0];
		}
	}

	im1 = YUV2RGB0(im1YUV);
	//////////////////////////////////////////////////
	//cout << im1.at<Vec3b>(0, 0)[2] << endl;
	//////////////////////////////////////////////////
	//im1YUV.convertTo(im1YUV,CV_8UC3); //???
	//cv::cvtColor(im1YUV, im1, COLOR_YUV2BGR);

	//
	//im1 = YUV2YUVMat(YCBCR);
	//im1.convertTo(im1, CV_32FC3);
	//RGB rgb_im1 = YUV2RGB(im1);
	//im1 = RGB2RGBMat(rgb_im1);

	//cout << imgpathlist.size() << endl;

	for (int ii = 1; ii < num_img; ii++)
	{
		//im2 = imread(imgpathlist[ii]);
		//im2.convertTo(im2, CV_32FC3);
		//cout << (int)im2.at<Vec3b>(0, 0)[2] << endl;
		//cout << (int)im2.at<Vec3b>(0, 0)[1] << endl;
		//cout << (int)im2.at<Vec3b>(0, 0)[0] << endl;



		im2YUV = RGB2YUV0(im2);
		//cout << im2YUV.at<Vec3f>(0, 0)[0] << endl;
		//cout << im2YUV.at<Vec3f>(0, 0)[1] << endl;
		//cout << im2YUV.at<Vec3f>(0, 0)[2] << endl;


		//cv::cvtColor(im2, im2YUV, COLOR_BGR2YUV);
		//im2YUV.convertTo(im2YUV, CV_32FC3);

		//cout << (int)im2YUV.at<Vec3b>(0, 0)[0] << endl;

		for (int i = 0; i < im2.rows; i++)
		{
			for (int j = 0; j < im2.cols; j++)
			{
				im2YUV.at<Vec3f>(i, j)[0] = pow(im2YUV.at<Vec3f>(i, j)[0], gamma[ii]);
				//cout << (int)im2YUV.at<Vec3b>(0, 0)[0] << endl;
				im2YUV.at<Vec3f>(i, j)[1] = im2YUV.at<Vec3f>(i, j)[1] * acb[ii];
				im2YUV.at<Vec3f>(i, j)[2] = im2YUV.at<Vec3f>(i, j)[2] * acr[ii];
			}
		}
		//cout << im2YUV.at<Vec3f>(0, 0)[0] << endl;
		//cout << im2YUV.at<Vec3f>(0, 0)[1] << endl;
		//cout << im2YUV.at<Vec3f>(0, 0)[2] << endl;



		//im2YUV.convertTo(im2YUV, CV_8UC3);
		//cv::cvtColor(im2YUV, im2, COLOR_YUV2BGR);

		//cout << (int)im2YUV.at<Vec3b>(0, 0)[0] << endl;

		im2 = YUV2RGB0(im2YUV);

		//im1.convertTo(im1, CV_32FC3);
		//im2.convertTo(im2, CV_32FC3);

		////////////////////RGB YUV 数值不匹配？？   im2    
		//cout << (int)im2.at<Vec3b>(0, 0)[2] << endl;
		//cout << (int)im2.at<Vec3b>(0, 0)[1] << endl;
		//cout << (int)im2.at<Vec3b>(0, 0)[0] << endl;

		//cout << im2.empty() << endl;

		//seam出溢出 -200
		//shift_x = W - 200;
		//shift_y = 0;

		overlap_w = im1.cols - shift_x;
		overlap_h = im1.rows - shift_y;

		Mat	overlap_img1;
		Mat overlap_img2;

		//cout << overlap_img1.rows << endl;  //行数√
		//cout << overlap_img1.at<Vec3b>(202, 0)[0] << endl;


		if (shift_y >= 0)
		{
			overlap_img1 = im1(Rect(shift_x, 0, overlap_w, overlap_h))/*.clone()*/;
			//overlap_im1 = img1(1:overlap_h, (shift_x + 1) : size(im1, 2), : );


			overlap_img2 = im2(Rect(0, 0, overlap_w, overlap_h))/*.clone()*/;
			//overlap_im2 = im2((size(im2, 1) - overlap_h + 1) : size(im2, 1), 1 : overlap_w, : );

		}
		else
		{
			overlap_img1 = im1(Rect(shift_x, 0, gomy.w, gomy.h)).clone();
			//overlap_im1 = im1(1:overlap_h, (shift_x + 1) : size(im1, 2), : );

			overlap_img2 = im2(Rect(0, 0, gomy.w, gomy.h)).clone();
			//overlap_im2 = im2((size(im2, 1) - overlap_h + 1) : size(im2, 1), 1 : overlap_w, : );
		}

		//cout << (int)overlap_img2.at<Vec3b>(0, 0)[2] << endl;
		//cout << (int)overlap_img2.at<Vec3b>(0, 0)[1] << endl;
		//cout << (int)overlap_img2.at<Vec3b>(0, 0)[0] << endl;
		//cout << overlap_img1.at<Vec3b>(0, 0)[0] <<
		//	overlap_img1.size() << endl;

		//overlap_img1.convertTo(overlap_img1, CV_32FC3);
		//overlap_img2.convertTo(overlap_img2, CV_32FC3);



		vector<vector<float> > squared_diff(overlap_h + 1, vector<float>(overlap_w + 1, 0));
		vector<vector<float> > accumulate_D(overlap_h + 1, vector<float>(overlap_w + 1, 0));
		vector<vector<float> > color_diff(overlap_h + 1, vector<float>(3, 0));
		vector<float> neighbor;

		for (int i = 1; i <= overlap_h; i++)
		{
			for (int j = 1; j <= overlap_w; j++)
			{
				squared_diff[i][j] = pow(/*(int)*/overlap_img1.at<Vec3f>(i - 1, j - 1)[2] - /*(int)*/overlap_img2.at<Vec3f>(i - 1, j - 1)[2], 2);

				//cout << (int)overlap_img1.at<Vec3b>(0, 0)[2] << endl;
				//cout << (int)overlap_img2.at<Vec3b>(0, 0)[2] << endl;

				if (i == 1)
				{
					accumulate_D[i][j] = squared_diff[i][j];

				}
				else
				{
					if (j == 1)
					{
						neighbor.push_back(accumulate_D[i - 1][j]);
						neighbor.push_back(accumulate_D[i - 1][j + 1]);
					}
					else
					{
						if (j == overlap_w)
						{
							neighbor.push_back(accumulate_D[i - 1][j - 1]);
							neighbor.push_back(accumulate_D[i - 1][j]);
						}
						else
						{
							neighbor.push_back(accumulate_D[i - 1][j - 1]);
							neighbor.push_back(accumulate_D[i - 1][j]);
							neighbor.push_back(accumulate_D[i - 1][j + 1]);
						}
					}
					accumulate_D[i][j] = squared_diff[i][j] + *min_element(neighbor.begin(), neighbor.end());
					neighbor.clear();
				}
			}
		}


		//for (int j = 0; j < overlap_w; j++)
		//{
		//	/*cout << squared_diff[i][j] << endl;*/
		//	cout << accumulate_D[overlap_h - 1][j] << endl;
		//}

		minhang mindst;
		float mindiff;
		vector<int> seamindex;
		seamindex.resize(overlap_h + 1, 0);

		mindst = minhh(accumulate_D, overlap_h);
		mindiff = mindst.mindiff;
		seamindex[overlap_h] = mindst.seamindex;

		//cout << mindst.mindiff << endl;
		//cout << mindst.seamindex << endl;

		////有bug   //边缘取值达到0 
		for (int i = overlap_h - 1; i >= 1; i--)
		{
			if (seamindex[i + 1] == 1)
			{
				if ((accumulate_D[i][seamindex[i + 1]] + squared_diff[i + 1][seamindex[i + 1]]) == mindiff)
				{
					seamindex[i] = seamindex[i + 1];
					mindiff = accumulate_D[i][seamindex[i]];
				}
				else
				{
					if ((accumulate_D[i][seamindex[i + 1] + 1] + squared_diff[i + 1][seamindex[i + 1]]) == mindiff)
					{
						seamindex[i] = seamindex[i + 1] + 1;
						mindiff = accumulate_D[i][seamindex[i]];
					}
				}
			}
			else
			{
				if (seamindex[i + 1] == overlap_w)
				{
					if ((accumulate_D[i][seamindex[i + 1] - 1] + squared_diff[i + 1][seamindex[i + 1]]) == mindiff)
					{
						seamindex[i] = seamindex[i + 1] - 1;
						mindiff = accumulate_D[i][seamindex[i]];
					}
					else
					{
						if ((accumulate_D[i][seamindex[i + 1]] + squared_diff[i + 1][seamindex[i + 1]]) == mindiff)
						{
							seamindex[i] = seamindex[i + 1];
							mindiff = accumulate_D[i][seamindex[i]];
						}
					}
				}
				else
				{
					if (accumulate_D[i][seamindex[i + 1] - 1] + squared_diff[i + 1][seamindex[i + 1]] == mindiff)
					{
						seamindex[i] = seamindex[i + 1] - 1;
						mindiff = accumulate_D[i][seamindex[i]];
					}
					else
					{
						if ((accumulate_D[i][seamindex[i + 1]] + squared_diff[i + 1][seamindex[i + 1]]) == mindiff)
						{
							seamindex[i] = seamindex[i + 1];
							mindiff = accumulate_D[i][seamindex[i]];
						}
						else
						{
							if ((accumulate_D[i][seamindex[i + 1] + 1] + squared_diff[i + 1][seamindex[i + 1]]) == mindiff)
							{
								seamindex[i] = seamindex[i + 1] + 1;
								mindiff = accumulate_D[i][seamindex[i]];
							}
						}
					}
				}

			}
		}


		Mat temp(overlap_h, overlap_w, CV_32FC3);

		for (int i = 0; i < overlap_h; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				for (int k = 0; k < seamindex[i + 1] - 1; k++)
				{
					temp.at<Vec3f>(i, k)[j] = overlap_img1.at<Vec3f>(i, k)[j];
				}

				//cout << "seamindex:" << seamindex[i] << endl;

				temp.at<Vec3f>(i, seamindex[i + 1])[j] = 0;

				for (int k = seamindex[i + 1] + 1; k < overlap_w; k++)
				{
					temp.at<Vec3f>(i, k)[j] = overlap_img2.at<Vec3f>(i, k)[j];
				}
			}

			for (int j = 0; j < 3; j++)
			{
				color_diff[i][j] = /*(int)*/overlap_img1.at<Vec3f>(i, seamindex[i + 1])[j] - /*(int)*/overlap_img2.at<Vec3f>(i, seamindex[i])[j];
				//cout << "color_diff:" << color_diff[i][j] << endl;


			}
		}

		//double total = 0;
		//double blendvalue[3][1] = { 0 };

		vector<float> dist;
		dist.resize(overlap_h);

		for (int i = 0; i < overlap_h; i++)
		{
			for (int j = seamindex[i + 1] + 1; j < im2.cols; j++)
			{
				float total = 0;
				float blendvalue[3] = { 0 };

				for (int k = 0; k < overlap_h; k++)
				{
					dist[k] = abs(i - k) + abs(j - (seamindex[k + 1] /*- 1*/));
					//cout << "dist[k]:" << dist[k] << endl;
					//cout << "color_diff:" << color_diff[k][0] << endl;

					blendvalue[0] = blendvalue[0] + 1 / dist[k] * color_diff[k][0];
					blendvalue[1] = blendvalue[1] + 1 / dist[k] * color_diff[k][1];
					blendvalue[2] = blendvalue[2] + 1 / dist[k] * color_diff[k][2];

					total = total + 1 / dist[k];


				}

				//cout << "total:" << total << endl;
				//cout << (int)im2.at<Vec3b>(i, j)[0] << endl;
				//cout <<"blendvalue:"<< blendvalue[2] << endl;
				im2.at<Vec3f>(i, j)[0] = /*(int)*/im2.at<Vec3f>(i, j)[0] + /*(int)*/(blendvalue[2] / total);
				im2.at<Vec3f>(i, j)[1] = /*(int)*/im2.at<Vec3f>(i, j)[1] + /*(int)*/(blendvalue[1] / total);
				im2.at<Vec3f>(i, j)[2] = /*(int)*/im2.at<Vec3f>(i, j)[2] + /*(int)*/(blendvalue[0] / total);

				//RGB 2 1 0

			}
		}

		im1.convertTo(im1, CV_8UC3);
		im2.convertTo(im2, CV_8UC3);

		*im1p = im1.clone();
		*im2p = im2.clone();

		//if (ii == 1)cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\1.jpg", im1);
		//cv::imwrite("C:\\Users\\10048\\Desktop\\pic\\2.jpg", im2);

	}

}
