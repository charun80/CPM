#ifndef _GaussianPyramid_h
#define _GaussianPyramid_h

#include "Image.h"

namespace cpm
{


template <class T>
class ImagePyramid
{
private:
	Image<T>* ImPyramid;
	int nLevels;
	float fRatio;
public:
	ImagePyramid(void){ ImPyramid = NULL; };
	~ImagePyramid(void){if(ImPyramid != NULL) delete[]ImPyramid;};
	inline Image<T>& operator[](int level) { return ImPyramid[level]; };
	void ConstructPyramid(const FImage& image, float ratio = 0.8, int minWidth = 30);
	void ConstructPyramidLevels(const FImage& image, float ratio = 0.8, int _nLevels = 2);
	void displayTop(const char* filename){ ImPyramid[nLevels - 1].imwrite(filename); };
	inline int nlevels() const {return nLevels;};
	inline float ratio() const { return fRatio; };
};

typedef ImagePyramid<float> FImagePyramid;


}  // namespace cpm


#include "ImagePyramid.hpp"

#endif
