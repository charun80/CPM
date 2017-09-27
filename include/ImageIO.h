#ifndef _ImageIO_h
#define _ImageIO_h

#include "opencv2/opencv.hpp"


namespace cpm
{


class ImageIO
{
public:
	enum ImageType{standard, derivative, normalized};
	ImageIO(void);
	~ImageIO(void);
public:
	template <class T>
	static bool loadImage(const char* filename,T*& pImagePlane,int& width,int& height, int& nchannels);
	template <class T>
	static bool saveImage(const char* filename,const T* pImagePlane,int width,int height, int nchannels,ImageType imtype = standard);
	template <class T>
	static void showImage(const char* winname, const T* pImagePlane, int width, int height, int nchannels, ImageType imtype = standard, int waittime = 1);
	template <class T>
	static void showGrayImageAsColor(const char* winname, const unsigned char* pImagePlane, int width, int height, T minV, T maxV, int waittime = 1);
	template <class T>
	static cv::Mat CvmatFromPixels(const T* pImagePlane, int width, int height, int nchannels, ImageType imtype = standard);
	template <class T>
	static void CvmatToPixels(const cv::Mat& cvInImg, T*& pOutImagePlane, int& width, int& height, int& nchannels);
private:
};


}  // namespace cpm


#include "ImageIO.hpp"

#endif  // _ImageIO_h
