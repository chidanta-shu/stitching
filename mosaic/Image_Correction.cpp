#include <Image_Correction.h>
#include <stdafx.h>
// get over mean yuv (GOMY)
inline GOMY get_overlap_mean_yuv(Mat img1, Mat img2, int shift_x, int shift_y)
{
	GOMY GOMY;
	GOMY.w = img1.cols - shift_x;
	GOMY.h = img1.rows - shift_y;

	Mat	overlap_img1;
	Mat overlap_img2;

	YUV YUV1;
	YUV YUV2;

	//vector<Mat> channelrgb1;
	//vector<Mat> channelrgb2;

	double	Ysum1 = 0;
	double	Ysum2 = 0;
	double	Cbsum1 = 0;
	double	Cbsum2 = 0;
	double	Crsum1 = 0;
	double	Crsum2 = 0;
	int	count = 0;

	//
	if (shift_y >= 0)
	{
		overlap_img1 = img1(Rect( GOMY.w - 1, 0, shift_x, GOMY.h));
		//overlap_im1 = img1(1:overlap_h, (shift_x + 1) : size(im1, 2), : );

		overlap_img2 = img2(Rect( 0, 0, GOMY.w, GOMY.h));
		//overlap_im2 = im2((size(im2, 1) - overlap_h + 1) : size(im2, 1), 1 : overlap_w, : );
	}
	else
	{
		overlap_img1 = img1(Rect(GOMY.w - 1, 0, shift_x, GOMY.h));
		//overlap_im1 = im1(1:overlap_h, (shift_x + 1) : size(im1, 2), : );

		overlap_img2 = img2(Rect(0, 0, GOMY.w, GOMY.h));
	    //overlap_im2 = im2((size(im2, 1) - overlap_h + 1) : size(im2, 1), 1 : overlap_w, : );
	}

	//split(overlap_img1, channelrgb1);
	//split(overlap_img1, channelrgb2);

	
	YUV1 = RGB2YUV(overlap_img1);
	YUV2 = RGB2YUV(overlap_img2);


	//YCBCR1 = RGB2YUV(overlap_img1);
	//YCBCR2 = RGB2YUV(overlap_img2);

	//YCBCR1 = double(YCBCR1);
	//YCBCR2 = double(YCBCR2);


	// R-> Y ,B-> Cb ,G-> Cr 

	for (int i = 0; i < GOMY.h; i++)
	{
		for (int j = 0; j < GOMY.w; j++)
		{
			count = count + 1;
			Ysum1 = Ysum1 + YUV1.Y[i][j] * YUV1.Y[i][j];
			Ysum2 = Ysum2 + YUV2.Y[i][j] * YUV2.Y[i][j];
			//Ysum1 = Ysum1 +  YCBCR1 (i,j,1)^2;%^2.2;
			//Ysum2 = Ysum2 + YCBCR2(i, j, 1) ^ 2; %^ 2.2;

			Cbsum1 = Cbsum1 + YUV1.U[i][j];
			Cbsum2 = Cbsum2 + YUV2.U[i][j];
			//Cbsum1 = Cbsum1 + YCBCR1(i, j, 2);
			//Cbsum2 = Cbsum2 + YCBCR2(i, j, 2);

			Crsum1 = Crsum1 + YUV1.V[i][j];
			Crsum2 = Crsum2 + YUV2.V[i][j];
			//Crsum1 = Crsum1 + YCBCR1(i, j, 3);
			//Crsum2 = Crsum2 + YCBCR2(i, j, 3);

		}
	}

	GOMY.Bi[0] = log(Ysum1 / count);
	GOMY.Bi.push_back(log(Ysum2 / count));

	GOMY.Sicb[0] = Cbsum1 / count;
	GOMY.Sicb.push_back(Cbsum2 / count);

	GOMY.Sicr[0] = Crsum1 / count;
	GOMY.Sicr.push_back(Crsum1 / count);

		

	return GOMY;
}

//  gamma ֵ 
vector<double> getCorrParameters(vector<vector<double> > B)
{
	//
	double sigma_n = 2; //  /255;
	double sigma_g = 0.5; //  /255;

	vector<double> gamma0;

	vector<double> a = { 0 };
	vector<double> b = { ( B[1][1] * B[1][1] ) / ( sigma_n * sigma_n ) + 1 / ( sigma_g * sigma_g ) };
	vector<double> c = { ( -B[1][1] * B[1][2] ) / ( sigma_n * sigma_n ) };
	vector<double> d = { 1 / ( sigma_g * sigma_g ) };
	
	//a.push_back(0);

	vector<double> bn;
	vector<double> an;
	vector<double> dn;

	bn[1] = b[1];
	an[1] = a[1];
	dn[1] = d[1];


	// 初值 i = 2 ?
	for (int i = 1; i < num_img; i++)
	{
		a[i] = -B[i - 1][1] * B[i - 1][2] / ( sigma_n * sigma_n );
		if ( i==num_img )
		{
			b[i] = ( B[i - 1][2] * B[i - 1][2] ) / ( sigma_n * sigma_n ) + 1 / ( sigma_g * sigma_g );
			c[i] = 0;
		}
		else
		{
			b[i] = ( B[i - 1][2] * B[ i- 1][2] + B[i][1] * B[i][1] ) / ( sigma_n * sigma_n ) + 1 / ( sigma_g * sigma_g );
			c[i] = -B[i][1] * B[i][2] / ( sigma_n * sigma_n );
		}

		d[i] = 1 / ( sigma_g * sigma_g );

		an[i] = a[i] / bn[i - 1];
		bn[i] = b[i] - a[i] / bn[i - 1] * c[i - 1];
		dn[i] = d[i] - an[i] * dn[i - 1];
	}
	
	gamma0[num_img] = dn[num_img] / bn[num_img];

	for (int i = (num_img - 1); i >= 0; i--)
	{
		gamma0[i] = ( dn[i] - c[i] * gamma0[i + 1] ) / bn[i];
	}

	return gamma0;
}

Mat RGB2RGBMat(RGB src)
{
	cv::Mat	dst3(sizeof(src.R[0]) / sizeof(src.R[0][0]), sizeof(src.R) / sizeof(src.R[0]), CV_32FC3);

	Mat RR(sizeof(src.R[0]) / sizeof(src.R[0][0]), sizeof(src.R) / sizeof(src.R[0]), CV_32FC1);
	Mat GG(sizeof(src.R[0]) / sizeof(src.R[0][0]), sizeof(src.R) / sizeof(src.R[0]), CV_32FC1);
	Mat BB(sizeof(src.R[0]) / sizeof(src.R[0][0]), sizeof(src.R) / sizeof(src.R[0]), CV_32FC1);

	for (int i = 0; i < sizeof(src.R) / sizeof(src.R[0]); i++)
	{
		for (int j = 0; j < sizeof(src.R[0]) / sizeof(src.R[0][0]); j++)
		{
			RR.at<double>(i, j) = src.R[i][j];
			GG.at<double>(i, j) = src.G[i][j];
			BB.at<double>(i, j) = src.B[i][j];

		}
	}

	Mat dst[3] = { BB, GG, RR };
	merge(dst, 3, dst3);

	return dst3;
}

Mat YUV2YUVMat(YUV YCBCR)
{
	cv::Mat	YUV_H(sizeof(YCBCR.Y[0]) / sizeof(YCBCR.Y[0][0]), sizeof(YCBCR.Y) / sizeof(YCBCR.Y[0]), CV_32FC3);
    Mat YY(sizeof(YCBCR.Y[0]) / sizeof(YCBCR.Y[0][0]), sizeof(YCBCR.Y) / sizeof(YCBCR.Y[0]), CV_32FC1);
    Mat UU(sizeof(YCBCR.Y[0]) / sizeof(YCBCR.Y[0][0]), sizeof(YCBCR.Y) / sizeof(YCBCR.Y[0]), CV_32FC1);
    Mat VV(sizeof(YCBCR.Y[0]) / sizeof(YCBCR.Y[0][0]), sizeof(YCBCR.Y) / sizeof(YCBCR.Y[0]), CV_32FC1);

    for (int i = 0; i < sizeof(YCBCR.Y) / sizeof(YCBCR.Y[0]); i++)
    {
	    for (int j = 0; j < sizeof(YCBCR.Y[0]) / sizeof(YCBCR.Y[0][0]); j++)
	    {
		    YY.at<double>(i, j) = YCBCR.Y[i][j];
	    	UU.at<double>(i, j) = YCBCR.Y[i][j];
	    	VV.at<double>(i, j) = YCBCR.Y[i][j];

	    }
    }
    Mat YUVH[3] = { YY, UU, VV };


    merge(YUVH, 3, YUV_H);


	return YUV_H;
}

minhang minhh(vector<vector<double> > a,int m)
{
	minhang b;
	int w = sizeof(a[0]) / sizeof(a[0][0]);
	int h = sizeof(a) / sizeof(a[0]);
	vector<double> temp;
	if (m < h)
	{
		for (int i = 0; i < w; i++)
		{
			temp[i] = a[m - 1][i];
		}
	}
	b.mindiff = *min_element(temp.begin(), temp.end());
	double min = temp[0];

	for (int i = 0; i < w; i++)
	{
		if (min < temp[i])
		{
			min = temp[i];
			b.seamindex = i;
		}
	}
	return b;
}

// YUV
inline YUV RGB2YUV(Mat src)
{
	YUV YUV;
	vector<Mat> dst;
	vector<vector<double> > R;
	vector<vector<double> > G;
    vector<vector<double> > B;

	split(src, dst);

	R = dst.at(2);
	G = dst.at(1);
    B = dst.at(0);
	int c = src.cols;
	int r = src.rows;


	for (int i = 0; i < r; i++)
	{
		for (int j = 0; j < c; j++)
		{
            YUV.Y[i][j] = 0.299 * R[i][j] + 0.587 * G[i][j] + 0.114 * B[i][j];
			YUV.U[i][j] = -0.14713 * R[i][j] - 0.28886 * G[i][j] + 0.436 * B[i][j];
			YUV.V[i][j] = 0.615 * R[i][j] - 0.51499 * G[i][j] - 0.10001 * B[i][j];

		}
	}

	//y = floor(0.299 * r + 0.587 * g + 0.114 * b);
	//u = floor(-0.14713 * r - 0.28886 * g + 0.436 * b);
	//v = floor(0.615 * r - 0.51499 * g - 0.10001 * b);

	///////////////////////////////////////////
	//cv::Mat imageY(src.rows, src.cols, 1);
	//cv::Mat imageU(src.rows, src.cols, 1);
	//cv::Mat imageV(src.rows, src.cols, 1);

	//cv::Mat imageYUV;
	//cv::cvtColor(src, imageYUV, COLOR_BGR2YUV);
	//std::vector<Mat> mv;
	//split(src, (vector<Mat>&)mv);

	//imageY = mv[0].clone();
	//imageU = mv[1].clone();
	//imageV = mv[2].clone();
	////////////////////////////////////////////
	return YUV;
}

// RGB 未完成
inline RGB YUV2RGB(Mat src)
{
	RGB RGB;
	vector<Mat> dst;
	vector<vector<double> > Y;
	vector<vector<double> > U;
	vector<vector<double> > V;

	split(src, dst);

	Y = dst.at(0);
	U = dst.at(1);
	V = dst.at(2);

	int c = src.cols;
	int r = src.rows;

	for (int i = 0; i < c; i++)
	{
		for (int j = 0; j < r; j++)
		{
			RGB.R[i][j] =  Y[i][j] + 1.13983 * V[i][j];
			RGB.G[i][j] = Y[i][j] - 0.39465 * U[i][j] - 0.5806 * V[i][j];
			RGB.B[i][j] = Y[i][j] + 2.03211 * U[i][j];

		}
	}

	//cv::Mat imageR(src.rows, src.cols, 1);
	//cv::Mat imageG(src.rows, src.cols, 1);
	//cv::Mat imageB(src.rows, src.cols, 1);

	//cv::Mat imageRGB;
	//cv::cvtColor(src, imageRGB, COLOR_YUV420sp2RGB);
	//std::vector<Mat> mv;
	//split(src, (vector<Mat>&)mv);

	//imageR = mv[0].clone();
	//imageG = mv[1].clone();
	//imageB = mv[2].clone();



	return RGB;
}
