#ifndef _IMAGEFEATURES_H
#define _IMAGEFEATURES_H

#include "Image.h"
#include <cmath>
#include <memory.h>
#include <vector>


namespace cpm
{


class ImageFeature
{
public:
	ImageFeature(void);
	~ImageFeature(void);

	template <class T>
	static void imSIFT(const Image<T>& imsrc, UCImage& imsift, int cellSize = 2, int stepSize = 1, bool IsBoundaryIncluded = false, int nBins = 8);

	template <class T>
	static void imSIFT(const Image<T>& imsrc, UCImage& imsift, const std::vector<int> cellSizeVect, int stepSize = 1, bool IsBoundaryIncluded = false, int nBins = 8);

};


}  // namespace cpm


#include "ImageFeature.hpp"

#endif // _IMAGEFEATURES_H
