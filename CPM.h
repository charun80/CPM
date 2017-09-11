/*

Code of the Coarse-to-Fine PatchMatch, published at CVPR 2016 in
"Efficient Coarse-to-Fine PatchMatch for Large Displacement Optical Flow"
by Yinlin.Hu, Rui Song and Yunsong Li.

Email: huyinlin@gmail.com

Version 1.2

Copyright (C) 2016 Yinlin.Hu

Usages:

The program "cpm.exe" has been built and tested on Windows 7.

USAGE: cpm.exe img1Name img2Name outMatchName

Explanations:

The output of the program is a text file, which is in the format of "x1,y1,x2,y2"
corresponding to one match per line.

*/

#ifndef _CPM_H_
#define _CPM_H_


#ifdef PYLIB
    #undef _MATLAB
    #define _NO_IMAGE_IO
#endif

#include "ImagePyramid.h"


namespace cpm
{


struct sCPMParameters
{
	int m_Step_i;

	bool m_IsStereo_b;

	int   m_maxIters_i;
	float m_StopIterRatio_f;

	float m_pydRatio_f;

	int m_maxDisplacement_i;
	int m_checkThreshold_i;

	int m_borderWidth_i;

    // default parameters
	sCPMParameters()
        : m_Step_i( 3 )
        , m_IsStereo_b(false)
        , m_maxIters_i( 8 )
        , m_StopIterRatio_f( 0.05 )
        , m_pydRatio_f( 0.5 )
        , m_maxDisplacement_i( 400 )
        , m_checkThreshold_i( 3 )
        , m_borderWidth_i( 5 )
    {}
};



class CPM
{
public:

	CPM();
	explicit CPM( const sCPMParameters &f_Param );

	~CPM();

	int Matching(FImage& img1, FImage& img2, FImage& outMatches);


private:
	void CrossCheck(IntImage& seeds, FImage& seedsFlow, FImage& seedsFlow2, IntImage& kLabel2, int* valid, float th);
	float MatchCost(FImage& img1, FImage& img2, UCImage* im1f, UCImage* im2f, int x1, int y1, int x2, int y2);

	// a good initialization is already stored in bestU & bestV
	int Propogate(FImagePyramid& pyd1, FImagePyramid& pyd2, UCImage* pyd1f, UCImage* pyd2f, int level, float* radius, int iterCnt, IntImage* pydSeeds, IntImage& neighbors, FImage* pydSeedsFlow, float* bestCosts);
	void PyramidRandomSearch(FImagePyramid& pyd1, FImagePyramid& pyd2, UCImage* im1f, UCImage* im2f, IntImage* pydSeeds, IntImage& neighbors, FImage* pydSeedsFlow);
	void OnePass(FImagePyramid& pyd1, FImagePyramid& pyd2, UCImage* im1f, UCImage* im2f, IntImage& seeds, IntImage& neighbors, FImage* pydSeedsFlow);
	void UpdateSearchRadius(IntImage& neighbors, FImage* pydSeedsFlow, int level, float* outRadius);

	// minimum circle
	struct Point{
		double x, y;
	};
	double dist(Point a, Point b);
	Point intersection(Point u1, Point u2, Point v1, Point v2);
	Point circumcenter(Point a, Point b, Point c);
	// return the radius of the minimal circle
	float MinimalCircle(float* x, float*y, int n, float* centerX = NULL, float* centerY = NULL);

	//
	const sCPMParameters m_Param;


	IntImage _kLabels, _kLabels2;

	FImagePyramid _pyd1;
	FImagePyramid _pyd2;

	UCImage* _im1f;
	UCImage* _im2f;

	FImage* _pydSeedsFlow;
	FImage* _pydSeedsFlow2;

	IntImage _seeds;
	IntImage _seeds2;
	IntImage _neighbors;
	IntImage _neighbors2;

};

} //namespace cpm

#endif // _CPM_H_
