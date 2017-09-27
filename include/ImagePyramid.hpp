#ifndef _GaussianPyramid_hpp
#define _GaussianPyramid_hpp

#include "ImagePyramid.hpp"

namespace cpm
{


//---------------------------------------------------------------------------------------
// function to construct the pyramid
// this is the slow way
//---------------------------------------------------------------------------------------
/*void GaussianPyramid::ConstructPyramid(const DImage &image, float ratio, int minWidth)
{
// the ratio cannot be arbitrary numbers
if(ratio>0.98 || ratio<0.4)
ratio=0.75;
// first decide how many levels
nLevels=std::log((float)minWidth/image.width())/std::log(ratio);
if(ImPyramid!=NULL)
delete []ImPyramid;
ImPyramid=new DImage[nLevels];
ImPyramid[0].copyData(image);
float baseSigma=(1/ratio-1);
for(int i=1;i<nLevels;i++)
{
DImage foo;
float sigma=baseSigma*i;
image.GaussianSmoothing(foo,sigma,sigma*2.5);
foo.imresize(ImPyramid[i],std::pow(ratio,i));
}
}//*/

//---------------------------------------------------------------------------------------
// function to construct the pyramid
// this is the fast way
//---------------------------------------------------------------------------------------
template <class T>
void ImagePyramid<T>::ConstructPyramid(const FImage& image, float ratio /*= 0.8*/, int minWidth /*= 30*/)
{
	// the ratio cannot be arbitrary numbers
	if (ratio>0.98 || ratio<0.4)
		ratio = 0.75;
	// first decide how many levels
	nLevels = std::log((float)minWidth / image.width()) / std::log(ratio);
	fRatio = ratio;
	if (ImPyramid != NULL)
		delete[]ImPyramid;
	ImPyramid = new FImage[nLevels];
	ImPyramid[0].copyData(image);
	float baseSigma = (1 / ratio - 1);
	int n = std::log(0.25) / std::log(ratio);
	float nSigma = baseSigma*n;
	for (int i = 1; i<nLevels; i++)
	{
		FImage foo;
		if (i <= n)
		{
			float sigma = baseSigma*i;
			image.GaussianSmoothing(foo, sigma, sigma * 3);
			foo.imresize(ImPyramid[i], std::pow(ratio, i));
		}
		else
		{
			ImPyramid[i - n].GaussianSmoothing(foo, nSigma, nSigma * 3);
			float rate = (float)std::pow(ratio, i)*image.width() / foo.width();
			foo.imresize(ImPyramid[i], rate);
		}
	}
}

template <class T>
void ImagePyramid<T>::ConstructPyramidLevels(const FImage& image, float ratio /*= 0.8*/, int _nLevels /*= 2*/)
{
	// the ratio cannot be arbitrary numbers
	if (ratio>0.98 || ratio<0.4)
		ratio = 0.75;
	nLevels = _nLevels;
	fRatio = ratio;
	if (ImPyramid != NULL)
		delete[]ImPyramid;
	ImPyramid = new FImage[nLevels];
	ImPyramid[0].copyData(image);
	float baseSigma = (1 / ratio - 1);
	int n = std::log(0.25) / std::log(ratio);
	float nSigma = baseSigma*n;
	for (int i = 1; i<nLevels; i++)
	{
		FImage foo;
		if (i <= n)
		{
			float sigma = baseSigma*i;
			image.GaussianSmoothing(foo, sigma, sigma * 3);
			foo.imresize(ImPyramid[i], std::pow(ratio, i));
		}
		else
		{
			ImPyramid[i - n].GaussianSmoothing(foo, nSigma, nSigma * 3);
			float rate = (float)std::pow(ratio, i)*image.width() / foo.width();
			foo.imresize(ImPyramid[i], rate);
		}
	}
}

}  // namespace cpm

#endif // _GaussianPyramid_hpp
