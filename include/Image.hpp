#ifndef _IMAGE_HPP
#define _IMAGE_HPP

#include "Image.h"
#include <typeinfo>
#include <cstring>
#include "Memory.h"



namespace cpm
{


//------------------------------------------------------------------------------------------
// constructor
//------------------------------------------------------------------------------------------
template <class T>
Image<T>::Image()
{
	pData=NULL;
	imWidth=imHeight=nChannels=nPixels=nElements=0;
	IsDerivativeImage=false;
	colorType = DATA;
}

//------------------------------------------------------------------------------------------
// constructor with specified dimensions
//------------------------------------------------------------------------------------------
template <class T>
Image<T>::Image(int width,int height,int nchannels)
{
	imWidth=width;
	imHeight=height;
	nChannels=nchannels;
	computeDimension();
	pData=NULL;
	pData = (T*)xmalloc(nElements * sizeof(T));
	if(nElements>0)
		std::memset(pData,0,sizeof(T)*nElements);
	IsDerivativeImage=false;
	colorType = DATA;
}

template <class T>
Image<T>::Image(const T& value,int _width,int _height,int _nchannels)
{
	pData=NULL;
	allocate(_width,_height,_nchannels);
	setValue(value);
}

#ifndef _MATLAB
//template <class T>
//Image<T>::Image(const QImage& image)
//{
//	pData=NULL;
//	imread(image);
//}
#endif

template <class T>
void Image<T>::allocate(int width,int height,int nchannels)
{
	clear();
	imWidth=width;
	imHeight=height;
	nChannels=nchannels;
	computeDimension();
	pData=NULL;
	colorType = DATA;

	if(nElements>0)
	{
		pData = (T*)xmalloc(nElements * sizeof(T));
		std::memset(pData,0,sizeof(T)*nElements);
	}
}

template <class T>
template <class T1>
void Image<T>::allocate(const Image<T1> &other)
{
	allocate(other.width(),other.height(),other.nchannels());
	IsDerivativeImage = other.isDerivativeImage();
	colorType = other.colortype();
}

//------------------------------------------------------------------------------------------
// copy constructor
//------------------------------------------------------------------------------------------
template <class T>
Image<T>::Image(const Image<T>& other)
{
	imWidth=imHeight=nChannels=nElements=0;
	pData=NULL;
	copyData(other);
}

//------------------------------------------------------------------------------------------
// destructor
//------------------------------------------------------------------------------------------
template <class T>
Image<T>::~Image()
{
	if (pData != NULL){
		xfree(pData);
	}
}

//------------------------------------------------------------------------------------------
// clear the image
//------------------------------------------------------------------------------------------
template <class T>
void Image<T>::clear()
{
	if (pData != NULL){
		xfree(pData);
	}
	pData=NULL;
	imWidth=imHeight=nChannels=nPixels=nElements=0;
}

//------------------------------------------------------------------------------------------
// reset the image (reset the buffer to zero)
//------------------------------------------------------------------------------------------
template <class T>
void Image<T>::reset()
{
	if(pData!=NULL)
		std::memset(pData,0,sizeof(T)*nElements);
}

template <class T>
void Image<T>::setValue(const T &value)
{
	for(int i=0;i<nElements;i++)
		pData[i]=value;
}

template <class T>
void Image<T>::setValue(const T& value,int _width,int _height,int _nchannels)
{
	if(imWidth!=_width || imHeight!=_height || nChannels!=_nchannels)
		allocate(_width,_height,_nchannels);
	setValue(value);
}

template <class T>
void Image<T>::setPixel(int row, int col, T* valPtr)
{
	T* ptr = pixPtr(row, col);
	std::memcpy(ptr, valPtr, sizeof(T)*nChannels);
}

//------------------------------------------------------------------------------------------
// copy from other image
//------------------------------------------------------------------------------------------
template <class T>
void Image<T>::copyData(const Image<T>& other)
{
	imWidth=other.imWidth;
	imHeight=other.imHeight;
	nChannels=other.nChannels;
	nPixels=other.nPixels;
	IsDerivativeImage=other.IsDerivativeImage;
	colorType = other.colorType;

	if(nElements!=other.nElements)
	{
		nElements=other.nElements;	
		if(pData!=NULL)
			xfree(pData);
		pData=NULL;
		pData = (T*)xmalloc(nElements * sizeof(T));
	}
	if(nElements>0)
		std::memcpy(pData,other.pData,sizeof(T)*nElements);
}

template <class T>
template <class T1>
void Image<T>::copy(const Image<T1>& other)
{
	clear();

	imWidth=other.width();
	imHeight=other.height();
	nChannels=other.nchannels();
	computeDimension();

	IsDerivativeImage=other.isDerivativeImage();
	colorType = other.colortype();

	pData=NULL;
	pData = (T*)xmalloc(nElements * sizeof(T));
	const T1*& srcData=other.data();
	
	// adjust the data according to the data type
	bool isThisFloat = this->IsFloat();
	bool isOtherFloat = other.IsFloat();
	if (isThisFloat != isOtherFloat){
		if (isThisFloat){
			for (int i = 0; i < nElements; i++)
				pData[i] = srcData[i] / 255.0;
		}else{
			for (int i = 0; i < nElements; i++)
				pData[i] = srcData[i] * 255.0;
		}
	}else{
		for (int i = 0; i < nElements; i++)
			pData[i] = srcData[i];
	}
}

template <class T>
void Image<T>::im2float()
{
	if(IsFloat())
		for(int i=0;i<nElements;i++)
			pData[i]/=255;
}

//------------------------------------------------------------------------------------------
// override equal operator
//------------------------------------------------------------------------------------------
template <class T>
Image<T>& Image<T>::operator=(const Image<T>& other)
{
	copyData(other);
	return *this;
}

template <class T>
bool Image<T>::IsFloat() const
{
	if (typeid(T) == typeid(double) || typeid(T) == typeid(float))
		return true;
	else
		return false;
}

template <class T>
template <class T1>
bool Image<T>::matchDimension(const Image<T1>& image) const
{
	if(imWidth==image.width() && imHeight==image.height() && nChannels==image.nchannels())
		return true;
	else
		return false;
}

template <class T>
bool Image<T>::matchDimension(int width, int height, int nchannels) const
{
	if(imWidth==width && imHeight==height && nChannels==nchannels)
		return true;
	else
		return false;
}

//------------------------------------------------------------------------------------------
// function to move this image to a dest image at (x,y) with specified width and height
//------------------------------------------------------------------------------------------
template <class T>
template <class T1>
void Image<T>::moveto(Image<T1>& image,int x0,int y0,int width,int height)
{
	if(width==0)
		width=imWidth;
	if(height==0)
		height=imHeight;
	int NChannels=std::min(nChannels,image.nchannels());

	int x,y;
	for(int i=0;i<height;i++)
	{
		y=y0+i;
		if(y>=image.height())
			break;
		for(int j=0;j<width;j++)
		{
			x=x0+j;
			if(x>=image.width())
				break;
			for(int k=0;k<NChannels;k++)
				image.data()[(y*image.width()+x)*image.nchannels()+k]=pData[(i*imWidth+j)*nChannels+k];
		}
	}
}


//------------------------------------------------------------------------------------------
// resize the image
//------------------------------------------------------------------------------------------
template <class T>
bool Image<T>::imresize(float ratio, InterType type/* = INTER_LINEAR*/)
{
	if(pData==NULL)
		return false;

	T* pDstData;
	int DstWidth,DstHeight;
	DstWidth=(float)imWidth*ratio;
	DstHeight=(float)imHeight*ratio;
	pDstData = (T*)xmalloc(DstWidth*DstHeight*nChannels * sizeof(T));

	ImageProcessing::ResizeImage(pData,pDstData,imWidth,imHeight,nChannels,ratio,type);

	xfree(pData);
	pData=pDstData;
	imWidth=DstWidth;
	imHeight=DstHeight;
	computeDimension();
	return true;
}

template <class T>
template <class T1>
void Image<T>::imresize(Image<T1>& result, float ratio, InterType type/* = INTER_LINEAR*/) const
{
	int DstWidth,DstHeight;
	DstWidth=(float)imWidth*ratio;
	DstHeight=(float)imHeight*ratio;
	if(result.width()!=DstWidth || result.height()!=DstHeight || result.nchannels()!=nChannels)
		result.allocate(DstWidth,DstHeight,nChannels);
	else
		result.reset();
	ImageProcessing::ResizeImage(pData,result.data(),imWidth,imHeight,nChannels,ratio,type);
}

template <class T>
template <class T1>
void Image<T>::imresize(Image<T1>& result, int DstWidth, int DstHeight, InterType type/* = INTER_LINEAR*/) const
{
	if(result.width()!=DstWidth || result.height()!=DstHeight || result.nchannels()!=nChannels)
		result.allocate(DstWidth,DstHeight,nChannels);
	else
		result.reset();
	ImageProcessing::ResizeImage(pData,result.data(),imWidth,imHeight,nChannels,DstWidth,DstHeight,type);
}


template <class T>
void Image<T>::imresize(int dstWidth, int dstHeight, InterType type/* = INTER_LINEAR*/)
{
	Image foo(dstWidth,dstHeight,nChannels); // kfj: it should be Image instead of FImage
	ImageProcessing::ResizeImage(pData,foo.data(),imWidth,imHeight,nChannels,dstWidth,dstHeight,type);
	copyData(foo);
}

template <class T>
template <class T1>
void Image<T>::upSampleNN(Image<T1>& output,int ratio) const
{
	int width = imWidth*ratio;
	int height = imHeight*ratio;
	if(!output.matchDimension(width,height,nChannels))
		output.allocate(width,height,nChannels);
	for(int i =  0; i <imHeight; i++)
		for(int j = 0; j<imWidth; j++)
		{
			int offset = (i*imWidth+j)*nChannels;
			for(int ii = 0 ;ii<ratio;ii++)
				for(int jj=0;jj<ratio;jj++)
				{
					int offset1 = ((i*ratio+ii)*width+j*ratio+jj)*nChannels;
					for(int k = 0; k<nChannels; k++)
						output.data()[offset1+k] = pData[offset+k];
				}
		}
}

//------------------------------------------------------------------------------------------
// function of reading or writing images (uncompressed)
//------------------------------------------------------------------------------------------
template <class T>
bool Image<T>::saveImage(const char *filename) const
{
	std::ofstream myfile(filename,std::ios::out | std::ios::binary);
	if(myfile.is_open())
	{
		bool foo = saveImage(myfile);
		myfile.close();
		return foo;
	}
	else
		return false;
}

template <class T>
bool Image<T>::saveImage(std::ofstream& myfile) const
{
	char type[16];
	sprintf(type,"%s",typeid(T).name());
	myfile.write(type,16);
	myfile.write((char *)&imWidth,sizeof(int));
	myfile.write((char *)&imHeight,sizeof(int));
	myfile.write((char *)&nChannels,sizeof(int));
	myfile.write((char *)&IsDerivativeImage,sizeof(bool));
	myfile.write((char *)pData,sizeof(T)*nElements);
	return true;
}

template <class T>
bool Image<T>::loadImage(const char *filename)
{
	std::ifstream myfile(filename, std::ios::in | std::ios::binary);
	if(myfile.is_open())
	{
		bool foo = loadImage(myfile);
		myfile.close();
		return foo;
	}
	else
		return false;
}

template <class T>
bool Image<T>::loadImage(std::ifstream& myfile)
{
	char type[16];
	myfile.read(type,16);

	if(strcasecmp(type,"uint16")==0)
		sprintf(type,"unsigned short");
	if(strcasecmp(type,"uint32")==0)
		sprintf(type,"unsigned int");
	if(strcasecmp(type,typeid(T).name())!=0)
	{
		std::cerr<<"The type of the image is different from the type of the object!"<<std::endl;
		return false;
	}

	int width,height,nchannels;
	myfile.read((char *)&width,sizeof(int));
	myfile.read((char *)&height,sizeof(int));
	myfile.read((char *)&nchannels,sizeof(int));
	if(!matchDimension(width,height,nchannels))
		allocate(width,height,nchannels);
	myfile.read((char *)&IsDerivativeImage,sizeof(bool));
	myfile.read((char *)pData,sizeof(T)*nElements);
	
	return true;
}

//------------------------------------------------------------------------------------------
// function to load the image
//------------------------------------------------------------------------------------------
#ifndef _NO_IMAGE_IO
#ifndef _MATLAB

template <class T>
bool Image<T>::imread(const char* filename)
{
	clear();
	if(ImageIO::loadImage(filename,pData,imWidth,imHeight,nChannels))
	{
		computeDimension();
		colorType = BGR; // when we use qt or opencv to load the image, it's often BGR
		return true;
	}
	return false;
}


//template <class T>
//bool Image<T>::imread(const QString &filename)
//{
//	clear();
//	if(ImageIO::loadImage(filename,pData,imWidth,imHeight,nChannels))
//	{
//		computeDimension();
//		return true;
//	}
//	return false;
//}
//
//template <class T>
//void Image<T>::imread(const QImage& image)
//{
//	clear();
//	ImageIO::loadImage(image,pData,imWidth,imHeight,nChannels);
//	computeDimension();
//}
//
 //------------------------------------------------------------------------------------------
 // function to write the image 
 //------------------------------------------------------------------------------------------
template <class T>
bool Image<T>::imwrite(const char* filename) const
{
	ImageIO::ImageType type;
	if(IsDerivativeImage)
		type=ImageIO::derivative;
	else
		type=ImageIO::standard;
	return ImageIO::saveImage(filename,pData,imWidth,imHeight,nChannels,type);
}

template <class T>
bool Image<T>::imwrite(const char* filename,ImageIO::ImageType type) const
{
	return ImageIO::saveImage(filename,pData,imWidth,imHeight,nChannels,type);
}

template <class T>
void Image<T>::imshow(char* winname, int waittime /*= 1*/) const
{
	ImageIO::ImageType type;
	if (IsDerivativeImage){
		type = ImageIO::derivative;
	}else{
		type = ImageIO::standard;
		float minV = 0, maxV = 255;
		if (IsFloat()){
			maxV = 1.0;
		}
		for (int i = 0; i < nElements; i++){
			if (pData[i] < minV || pData[i] > maxV){
				type = ImageIO::normalized;
				break;
			}
		}
	}
	ImageIO::showImage(winname, pData, imWidth, imHeight, nChannels, type, waittime);
}

template <class T>
void Image<T>::imagesc(char* winname, int waittime /*= 1*/) const
{
	if (nChannels == 1){
		UCImage img;
		this->normalize(img);
		ImageIO::showGrayImageAsColor(winname, img.pData, imWidth, imHeight, 0, 255, waittime);
	}else{
		fprintf(stderr, "imagesc() only support one-channel image.");
	}
}

template <class T>
void Image<T>::ToLab() const
{
	if (colorType == BGR){
		ImageProcessing::BGR2Lab(pData, pData, imWidth, imHeight);
		colorType = LAB;
	}
}


template <class T>
template <class T1>
void Image<T>::ToLab(Image<T1>& image) const
{
	if (matchDimension(image) == false)
		image.allocate(imWidth, imHeight, nChannels);
	ImageProcessing::BGR2Lab(pData, image.pData, imWidth, imHeight);
}

//template <class T>
//bool Image<T>::imwrite(const QString &filename, ImageIO::ImageType imagetype, int quality) const
//{
//	return ImageIO::writeImage(filename,(const T*&)pData,imWidth,imHeight,nChannels,imagetype,quality);
//}
//
//template <class T>
//bool Image<T>::imwrite(const QString &filename, T min, T max, int quality) const
//{
//	return ImageIO::writeImage(filename,(const T*&)pData,imWidth,imHeight,nChannels,min,max,quality);
//}

#endif  // _MATLAB
#endif  // _NO_IMAGE_IO

//------------------------------------------------------------------------------------------
// function to get x-derivative of the image
//------------------------------------------------------------------------------------------
template <class T>
template <class T1>
void Image<T>::dx(Image<T1>& result,bool IsAdvancedFilter) const
{
	if(matchDimension(result)==false)
		result.allocate(imWidth,imHeight,nChannels);
	result.reset();
	result.setDerivative();
	T1*& data=result.data();
	int i,j,k,offset;
	if(IsAdvancedFilter==false)
		for(i=0;i<imHeight;i++)
			for(j=0;j<imWidth-1;j++)
			{
				offset=i*imWidth+j;
				for(k=0;k<nChannels;k++)
					data[offset*nChannels+k]=(T1)pData[(offset+1)*nChannels+k]-pData[offset*nChannels+k];
			}
	else
	{
		float xFilter[5]={1,-8,0,8,-1};
		for(i=0;i<5;i++)
			xFilter[i]/=12;
		ImageProcessing::hfiltering(pData,data,imWidth,imHeight,nChannels,xFilter,2);
	}
}

template <class T>
template <class T1>
Image<T1> Image<T>::dx(bool IsAdvancedFilter) const
{
	Image<T1> result;
	dx<T1>(result,IsAdvancedFilter);
	return result;
}

//------------------------------------------------------------------------------------------
// function to get y-derivative of the image
//------------------------------------------------------------------------------------------
template <class T>
template <class T1>
void Image<T>::dy(Image<T1>& result,bool IsAdvancedFilter) const
{
	if(matchDimension(result)==false)
		result.allocate(imWidth,imHeight,nChannels);
	result.setDerivative();
	T1*& data=result.data();
	int i,j,k,offset;
	if(IsAdvancedFilter==false)
		for(i=0;i<imHeight-1;i++)
			for(j=0;j<imWidth;j++)
			{
				offset=i*imWidth+j;
				for(k=0;k<nChannels;k++)
					data[offset*nChannels+k]=(T1)pData[(offset+imWidth)*nChannels+k]-pData[offset*nChannels+k];
			}
	else
	{
		float yFilter[5]={1,-8,0,8,-1};
		for(i=0;i<5;i++)
			yFilter[i]/=12;
		ImageProcessing::vfiltering(pData,data,imWidth,imHeight,nChannels,yFilter,2);
	}
}

template <class T>
template <class T1>
Image<T1> Image<T>::dy(bool IsAdvancedFilter) const
{
	Image<T1> result;
	dy<T1>(result,IsAdvancedFilter);
	return result;
}

//------------------------------------------------------------------------------------------
// function to compute the second order derivative 
//------------------------------------------------------------------------------------------
template <class T>
template <class T1>
void Image<T>::dxx(Image<T1> &image) const
{
	if(!matchDimension(image))
		image.allocate(imWidth,imHeight,nChannels);
	T1* pDstData=image.data();
	if(nChannels==1) // if there is only one image channel
		for(int i=0;i<imHeight;i++)
			for(int j=0;j<imWidth;j++)
			{
				int offset=i*imWidth+j;
				if(j==0)
				{
					pDstData[offset] = -pData[offset] + pData[offset + 1];
					continue;
				}
				if(j==imWidth-1)
				{
					pDstData[offset] = pData[offset - 1] - pData[offset];
					continue;
				}
				pDstData[offset] = pData[offset - 1] - 2 * pData[offset] + pData[offset + 1];
			}
	else
		for(int i=0;i<imHeight;i++)
			for(int j=0;j<imWidth;j++)
			{
				int offset=(i*imWidth+j)*nChannels;
				if(j==0)
				{
					for (int k = 0; k < nChannels; k++)
						pDstData[offset + k] = -pData[offset + k] + pData[offset + nChannels + k];
					continue;
				}
				if(j==imWidth-1)
				{
					for (int k = 0; k < nChannels; k++)
						pDstData[offset + k] = pData[offset - nChannels + k] - pData[offset + k];
					continue;
				}
				for (int k = 0; k < nChannels; k++)
					pDstData[offset + k] = pData[offset - nChannels + k] - 2 * pData[offset + k] + pData[offset + nChannels + k];
			}
}

template <class T>
template <class T1>
void Image<T>::dyy(Image<T1>& image) const
{
	if(!matchDimension(image))
		image.allocate(imWidth,imHeight,nChannels);
	T1* pDstData=image.data();
	if(nChannels==1)
		for(int i=0;i<imHeight;i++)
			for(int j=0;j<imWidth;j++)
			{
				int offset=i*imWidth+j;
				if(i==0)
				{
					pDstData[offset] = -pData[offset] + pData[offset + imWidth];
					continue;
				}
				if(i==imHeight-1)
				{
					pDstData[offset] = pData[offset - imWidth] - pData[offset];
					continue;
				}
				pDstData[offset] = pData[offset - imWidth] - 2 * pData[offset] + pData[offset + imWidth];
			}
	else
		for(int i=0;i<imHeight;i++)
			for(int j=0;j<imWidth;j++)
			{
				int offset=(i*imWidth+j)*nChannels;
				if(i==0)
				{
					for (int k = 0; k < nChannels; k++)
						pDstData[offset + k] = -pData[offset + k] + pData[offset + imWidth*nChannels + k];
					continue;
				}
				if(i==imHeight-1)
				{
					for (int k = 0; k < nChannels; k++)
						pDstData[offset + k] = pData[offset - imWidth*nChannels + k] - pData[offset + k];
					continue;
				}
				for (int k = 0; k < nChannels; k++)
					pDstData[offset + k] = pData[offset - imWidth*nChannels + k] - 2 * pData[offset + k] + pData[offset + imWidth*nChannels + k];
			}
}

template <class T>
template <class T1>
void Image<T>::dxy(Image<T1>& image) const
{
	if(!matchDimension(image))
		image.allocate(imWidth,imHeight,nChannels);
	T1* pDstData=image.data();

	//mask of second derivative
	float pfilter2D[] = {0.25, 0, -0.25, 0, 0, 0, -0.25, 0, 0.25};

	ImageProcessing::filtering(pData, pDstData, imWidth, imHeight, nChannels, pfilter2D, 1);
}

//------------------------------------------------------------------------------------------
// function for fast laplacian computation
//------------------------------------------------------------------------------------------
template <class T>
template <class T1>
void Image<T>::laplacian(Image<T1> &image) const
{
	if(!matchDimension(image))
		image.allocate(*this);
	image.setDerivative(true);
	ImageProcessing::Laplacian(pData,image.data(),imWidth,imHeight,nChannels);
}


//------------------------------------------------------------------------------------------
// function to compute the gradient magnitude of the image
//------------------------------------------------------------------------------------------
template <class T>
template <class T1>
void Image<T>::gradientmag(Image<T1> &image) const
{
	if(image.width()!=imWidth || image.height()!=imHeight)
		image.allocate(imWidth,imHeight);
	FImage Ix,Iy;
	dx(Ix,true);
	dy(Iy,true);
	float temp;
	float* imagedata=image.data();
	const float *Ixdata=Ix.data(),*Iydata=Iy.data();
	for(int i=0;i<nPixels;i++)
	{
		temp=0;
		int offset=i*nChannels;
		for(int k=0;k<nChannels;k++)
		{
			temp+=Ixdata[offset+k]*Ixdata[offset+k];
			temp+=Iydata[offset+k]*Iydata[offset+k];
		}
		imagedata[i]=std::sqrt(temp);
	}
}

//------------------------------------------------------------------------------------------
// function to do Gaussian smoothing
//------------------------------------------------------------------------------------------
template <class T>
void Image<T>::GaussianSmoothing(float sigma,int fsize) 
{
	Image<T> foo;
	GaussianSmoothing(foo,sigma,fsize);
	copy(foo);
}



template <class T>
template <class T1>
void Image<T>::GaussianSmoothing(Image<T1>& image,float sigma,int fsize) const 
{
	// constructing the 1D gaussian filter
	float* gFilter;
	gFilter=new float[fsize*2+1];
	float sum=0;
	sigma=sigma*sigma*2;
	for(int i=-fsize;i<=fsize;i++)
	{
		gFilter[i+fsize]=exp(-(float)(i*i)/sigma);
		sum+=gFilter[i+fsize];
	}
	for(int i=0;i<2*fsize+1;i++)
		gFilter[i]/=sum;

	// apply filtering
	imfilter_hv(image,gFilter,fsize,gFilter,fsize);

	delete []gFilter;
}

//------------------------------------------------------------------------------------------
// function to do Gaussian smoothing
//------------------------------------------------------------------------------------------
template <class T>
template <class T1>
void Image<T>::GaussianSmoothing_transpose(Image<T1>& image,float sigma,int fsize) const 
{
	Image<T1> foo;
	// constructing the 1D gaussian filter
	float* gFilter;
	gFilter=new float[fsize*2+1];
	float sum=0;
	sigma=sigma*sigma*2;
	for(int i=-fsize;i<=fsize;i++)
	{
		gFilter[i+fsize]=exp(-(float)(i*i)/sigma);
		sum+=gFilter[i+fsize];
	}
	for(int i=0;i<2*fsize+1;i++)
		gFilter[i]/=sum;

	// apply filtering
	imfilter_hv_transpose(image,gFilter,fsize,gFilter,fsize);

	delete gFilter;
}


//------------------------------------------------------------------------------------------
// function to smooth the image using a simple 3x3 filter
// the filter is [1 factor 1]/(factor+2), applied horizontally and vertically
//------------------------------------------------------------------------------------------
template <class T>
template <class T1>
void Image<T>::smoothing(Image<T1>& image,float factor)
{
	// build 
	float filter2D[9]={1,0,1,0, 0, 0,1, 0,1};
	filter2D[1]=filter2D[3]=filter2D[5]=filter2D[7]=factor;
	filter2D[4]=factor*factor;
	for(int i=0;i<9;i++)
		filter2D[i]/=(factor+2)*(factor+2);

	if(matchDimension(image)==false)
		image.allocate(imWidth,imHeight,nChannels);
	imfilter<T1>(image,filter2D,1);
}

template <class T>
template <class T1>
Image<T1> Image<T>::smoothing(float factor)
{
	Image<T1> result;
	smoothing(result,factor);
	return result;
}

template <class T>
void Image<T>::smoothing(float factor)
{
	Image<T> result(imWidth,imHeight,nChannels);
	smoothing(result,factor);
	copyData(result);
}

template <class T>
void Image<T>::MedianFiltering(int fsize)
{
	ImageProcessing::Medianfiltering(pData, pData, imWidth, imHeight, nChannels, fsize);
}

//------------------------------------------------------------------------------------------
//	 function of image filtering
//------------------------------------------------------------------------------------------
template <class T>
template <class T1>
void Image<T>::imfilter(Image<T1>& image,const float* filter,int fsize) const
{
	if(matchDimension(image)==false)
		image.allocate(imWidth,imHeight,nChannels);
	ImageProcessing::filtering(pData,image.data(),imWidth,imHeight,nChannels,filter,fsize);
}

template <class T>
template <class T1,class T2>
void Image<T>::imfilter(Image<T1>& image,const Image<T2>& kernel) const
{
	if(kernel.width()!=kernel.height())
	{
		std::cerr<<"Error in Image<T>::imfilter(Image<T1>& image,const Image<T2>& kernel)"<<std::endl;
		exit(-1);
	}
	int winsize = (kernel.width()-1)/2;
	imfilter(image,kernel.data(),winsize);
}

template <class T>
template <class T1>
Image<T1> Image<T>::imfilter(const float *filter, int fsize) const
{
	Image<T1> result;
	imfilter(result,filter,fsize);
	return result;
}

template <class T>
template <class T1>
void Image<T>::imfilter_h(Image<T1>& image,float* filter,int fsize) const
{
	if(matchDimension(image)==false)
		image.allocate(imWidth,imHeight,nChannels);
	ImageProcessing::hfiltering(pData,image.data(),imWidth,imHeight,nChannels,filter,fsize);
}

template <class T>
template <class T1>
void Image<T>::imfilter_v(Image<T1>& image,float* filter,int fsize) const
{
	if(matchDimension(image)==false)
		image.allocate(imWidth,imHeight,nChannels);
	ImageProcessing::vfiltering(pData,image.data(),imWidth,imHeight,nChannels,filter,fsize);
}


template <class T>
template <class T1>
void Image<T>::imfilter_hv(Image<T1> &image, const float *hfilter, int hfsize, const float *vfilter, int vfsize) const
{
	if(matchDimension(image)==false)
		image.allocate(imWidth,imHeight,nChannels);
	T1* pTempBuffer;
	pTempBuffer = (T1*)xmalloc(nElements * sizeof(T1));
	ImageProcessing::hfiltering(pData,pTempBuffer,imWidth,imHeight,nChannels,hfilter,hfsize);
	ImageProcessing::vfiltering(pTempBuffer,image.data(),imWidth,imHeight,nChannels,vfilter,vfsize);
	xfree(pTempBuffer);
}

template <class T>
template <class T1>
void Image<T>::imfilter_hv(Image<T1>& image,const Image<float>& hfilter,const Image<float>& vfilter) const
{
	int hfsize = (std::max(hfilter.width(),hfilter.height())-1)/2;
	int vfsize = (std::max(vfilter.width(),vfilter.height())-1)/2;
	imfilter_hv(image,hfilter.data(),hfsize,vfilter.data(),vfsize);
}

template<class T>
bool Image<T>::BoundaryCheck() const
{
	for(int i = 0;i<nElements;i++)
		if(!(pData[i]<1E10 && pData[i]>-1E10))
		{
			std::cerr<<"Error, bad data!"<<std::endl;
			return false;
		}
	return true;
}


//------------------------------------------------------------------------------------------
//	 function of image filtering transpose
//------------------------------------------------------------------------------------------
template <class T>
template <class T1>
void Image<T>::imfilter_transpose(Image<T1>& image,const float* filter,int fsize) const
{
	if(matchDimension(image)==false)
		image.allocate(imWidth,imHeight,nChannels);
	ImageProcessing::filtering_transpose(pData,image.data(),imWidth,imHeight,nChannels,filter,fsize);
}

template <class T>
template <class T1,class T2>
void Image<T>::imfilter_transpose(Image<T1>& image,const Image<T2>& kernel) const
{
	if(kernel.width()!=kernel.height())
	{
		std::cerr<<"Error in Image<T>::imfilter(Image<T1>& image,const Image<T2>& kernel)"<<std::endl;
		exit(-1);
	}
	int winsize = (kernel.width()-1)/2;
	imfilter_transpose(image,kernel.data(),winsize);
}

template <class T>
template <class T1>
Image<T1> Image<T>::imfilter_transpose(const float *filter, int fsize) const
{
	Image<T1> result;
	imfilter_transpose(result,filter,fsize);
	return result;
}

template <class T>
template <class T1>
void Image<T>::imfilter_h_transpose(Image<T1>& image,float* filter,int fsize) const
{
	if(matchDimension(image)==false)
		image.allocate(imWidth,imHeight,nChannels);
	ImageProcessing::hfiltering_transpose(pData,image.data(),imWidth,imHeight,nChannels,filter,fsize);
}

template <class T>
template <class T1>
void Image<T>::imfilter_v_transpose(Image<T1>& image,float* filter,int fsize) const
{
	if(matchDimension(image)==false)
		image.allocate(imWidth,imHeight,nChannels);
	ImageProcessing::vfiltering_transpose(pData,image.data(),imWidth,imHeight,nChannels,filter,fsize);
}


template <class T>
template <class T1>
void Image<T>::imfilter_hv_transpose(Image<T1> &image, const float *hfilter, int hfsize, const float *vfilter, int vfsize) const
{
	if(matchDimension(image)==false)
		image.allocate(imWidth,imHeight,nChannels);
	Image<T1> temp(imWidth,imHeight,nChannels);
	//imwrite("input.jpg");
	ImageProcessing::vfiltering_transpose(pData,temp.data(),imWidth,imHeight,nChannels,vfilter,vfsize);
	//temp.imwrite("temp.jpg");
	ImageProcessing::hfiltering_transpose(temp.data(),image.data(),imWidth,imHeight,nChannels,hfilter,hfsize);
}

template <class T>
template <class T1>
void Image<T>::imfilter_hv_transpose(Image<T1>& image,const Image<float>& hfilter,const Image<float>& vfilter) const
{
	int hfsize = (std::max(hfilter.width(),hfilter.height())-1)/2;
	int vfsize = (std::max(vfilter.width(),vfilter.height())-1)/2;
	imfilter_hv_transpose(image,hfilter.data(),hfsize,vfilter.data(),vfsize);
}

//------------------------------------------------------------------------------------------
//	 function for desaturation
//------------------------------------------------------------------------------------------
template <class T>
template <class T1>
void Image<T>::desaturate(Image<T1> &image) const
{
	if(nChannels!=3)
	{
		collapse(image);
		return;
	}
	if (!(image.width() == imWidth && image.height() == imHeight && image.nchannels() == 1))
		image.allocate(imWidth,imHeight,1);
	T1* data=image.data();
	int offset;
	for(int i=0;i<nPixels;i++)
	{
		offset=i*3;
		if(colorType == RGB)
			data[i]=(float)pData[offset]*.299+pData[offset+1]*.587+pData[offset+2]*.114;
		else
			data[i]=(float)pData[offset]*.114+pData[offset+1]*.587+pData[offset+2]*.299;
	}
}

template <class T>
void Image<T>::desaturate()
{
	Image<T> temp;
	desaturate(temp);
	copyData(temp);
}

template <class T>
template <class T1>
void Image<T>::collapse(Image<T1> &image,collapse_type type) const
{
	if(!(image.width()==imWidth && image.height()==imHeight && image.nchannels()==1))
		image.allocate(imWidth,imHeight,1);
	image.setDerivative(IsDerivativeImage);
	if(nChannels == 1)
	{
		image.copy(*this);
		return;
	}
	T1* data=image.data();
	int offset;
	float temp;
	for(int i=0;i<nPixels;i++)
	{
		offset=i*nChannels;
		switch(type){
			case collapse_average:
				temp=0;
				for(int j=0;j<nChannels;j++)
					temp+=pData[offset+j];
				data[i]=temp/nChannels;
				break;
			case collapse_max:
				data[i] = pData[offset];
				for(int j=1;j<nChannels;j++)
					data[i] = std::max(data[i],pData[offset+j]);
				break;
			case collapse_min:
				data[i] = pData[offset];
				for(int j = 1;j<nChannels;j++)
					data[i]=std::min(data[i],pData[offset+j]);
				break;
		}
	}
}

template <class T>
void Image<T>::collapse(collapse_type type)
{
	if(nChannels == 1)
		return;
	Image<T> result;
	collapse(result,type);
	copyData(result);
}

//------------------------------------------------------------------------------------------
//  function to concatenate two images
//------------------------------------------------------------------------------------------
template <class T>
template <class T1,class T2>
void Image<T>::concatenate(Image<T1> &destImage, const Image<T2> &addImage) const
{
	if(addImage.width()!=imWidth || addImage.height()!=imHeight)
	{
		destImage.copy(*this);
		return;
	}
	int extNChannels=nChannels+addImage.nchannels();
	if(destImage.width()!=imWidth || destImage.height()!=imHeight || destImage.nchannels()!=extNChannels)
		destImage.allocate(imWidth,imHeight,extNChannels);
	int offset;
	T1*& pDestData=destImage.data();
	const T2*& pAddData=addImage.data();
	for(int i=0;i<imHeight;i++)
		for(int j=0;j<imWidth;j++)
		{
			offset=i*imWidth+j;
			for(int k=0;k<nChannels;k++)
				pDestData[offset*extNChannels+k]=pData[offset*nChannels+k];
			for(int k=nChannels;k<extNChannels;k++)
				pDestData[offset*extNChannels+k]=pAddData[offset*addImage.nchannels()+k-nChannels];
		}
}

template <class T>
template <class T1,class T2>
void Image<T>::concatenate(Image<T1> &destImage, const Image<T2> &addImage,float ratio) const
{
	if(addImage.width()!=imWidth || addImage.height()!=imHeight)
	{
		destImage.copy(*this);
		return;
	}
	int extNChannels=nChannels+addImage.nchannels();
	if(destImage.width()!=imWidth || destImage.height()!=imHeight || destImage.nchannels()!=extNChannels)
		destImage.allocate(imWidth,imHeight,extNChannels);
	int offset;
	T1*& pDestData=destImage.data();
	const T2*& pAddData=addImage.data();
	for(int i=0;i<imHeight;i++)
		for(int j=0;j<imWidth;j++)
		{
			offset=i*imWidth+j;
			for(int k=0;k<nChannels;k++)
				pDestData[offset*extNChannels+k]=pData[offset*nChannels+k];
			for(int k=nChannels;k<extNChannels;k++)
				pDestData[offset*extNChannels+k]=pAddData[offset*addImage.nchannels()+k-nChannels]*ratio;
		}
}


template <class T>
template <class T1>
Image<T> Image<T>::concatenate(const Image<T1> &addImage) const
{
	Image<T> destImage;
	concatenate(destImage,addImage);
	return destImage;
}

//------------------------------------------------------------------------------------------
// function to separate the image into two
//------------------------------------------------------------------------------------------
template <class T>
template <class T1,class T2>
void Image<T>::separate(unsigned int firstNChannels, Image<T1> &image1, Image<T2> &image2) const
{
	image1.IsDerivativeImage=IsDerivativeImage;
	image2.IsDerivativeImage=IsDerivativeImage;

	if(firstNChannels>=nChannels)
	{
		image1=*this;
		image2.allocate(imWidth,imHeight,0);
		return;
	}
	if(firstNChannels==0)
	{
		image1.allocate(imWidth,imHeight,0);
		image2=*this;
		return;
	}
	int secondNChannels=nChannels-firstNChannels;
	if(image1.width()!=imWidth || image1.height()!=imHeight || image1.nchannels()!=firstNChannels)
		image1.allocate(imWidth,imHeight,firstNChannels);
	if(image2.width()!=imWidth || image2.height()!=imHeight || image2.nchannels()!=secondNChannels)
		image2.allocate(imWidth,imHeight,secondNChannels);

	for(int i=0;i<imHeight;i++)
		for(int j=0;j<imWidth;j++)
		{
			int offset=i*imWidth+j;
			for(int k=0;k<firstNChannels;k++)
				image1.pData[offset*firstNChannels+k]=pData[offset*nChannels+k];
			for(int k=firstNChannels;k<nChannels;k++)
				image2.pData[offset*secondNChannels+k-firstNChannels]=pData[offset*nChannels+k];
		}
}

template <class T>
void Image<T>::AddBorder(Image<T>& outImg, int borderWidth) const
{
	int extW = imWidth + 2*borderWidth;
	int extH = imHeight + 2*borderWidth;
	int ch = nChannels;
	if (!outImg.matchDimension(extW, extH, ch))
		outImg.allocate(extW, extH, ch);
#if 1
	for(int i=0;i<extH;i++){
		for(int j=0;j<extW;j++){
			int x = j-borderWidth;
			int y = i-borderWidth;
			x = ImageProcessing::EnforceRange(x, imWidth);
			y = ImageProcessing::EnforceRange(y, imHeight);
			T* dst = outImg.pData + (i*extW + j)*ch;
			T* src = pData + (y*imWidth + x)*ch;
			std::memcpy(dst, src, sizeof(T)*ch);
		}
	}
#else
	std::memset(outImg.pData, 0, outImg.nElements*sizeof(T));
	for (int i = 0; i < imHeight; i++){
		T* dstPtr = outImg.rowPtr(i + borderWidth) + borderWidth*ch;
		std::memcpy(dstPtr, pData + i*imWidth*ch, imWidth*ch*sizeof(T));
	}
#endif
}

//------------------------------------------------------------------------------------------
// function to separate the image into two
//------------------------------------------------------------------------------------------
template <class T>
template <class T1>
void Image<T>::getPatch(Image<T1>& patch,float x,float y,int wsize) const
{
	int wlength=wsize*2+1;
	if(patch.width()!=wlength || patch.height()!=wlength || patch.nchannels()!=nChannels)
		patch.allocate(wlength,wlength,nChannels);
	else
		patch.reset();
	ImageProcessing::getPatch(pData,patch.data(),imWidth,imHeight,nChannels,x,y,wsize);
}

//------------------------------------------------------------------------------------------
// function to crop an image
//------------------------------------------------------------------------------------------
template <class T>
template <class T1>
void Image<T>::crop(Image<T1>& patch,int Left,int Top,int Width,int Height) const
{
	if(patch.width()!=Width || patch.height()!=Height || patch.nchannels()!=nChannels)
		patch.allocate(Width,Height,nChannels);
	// make sure that the cropping is valid
	if(Left<0 || Top<0 || Left>=imWidth || Top>=imHeight)
	{
		std::cerr<<"The cropping coordinate is outside the image boundary!"<<std::endl;
		return;
	}
	if(Width<0 || Height<0 || Width+Left>imWidth || Height+Top>imHeight)
	{
		std::cerr<<"The patch to crop is invalid!"<<std::endl;
		return;
	}
	ImageProcessing::cropImage(pData,imWidth,imHeight,nChannels,patch.data(),Left,Top,Width,Height);
}

template <class T>
void Image<T>::flip_horizontal(Image<T>& image)
{
	if(!image.matchDimension(*this))
		image.allocate(*this);
	for(int i = 0;i<imHeight;i++)
		for(int j = 0;j<imWidth;j++)
		{
			int offset1 = (i*imWidth+j)*nChannels;
			int offset2 = (i*imWidth+imWidth-1-j)*nChannels;
			for(int k = 0;k<nChannels;k++)
				image.pData[offset2+k] = pData[offset1+k];
		}
}

template <class T>
void Image<T>::flip_horizontal()
{
	Image<T> temp(*this);
	flip_horizontal(temp);
	copyData(temp);
}

//------------------------------------------------------------------------------------------
// function to multiply image1, image2 and image3 to the current image
//------------------------------------------------------------------------------------------
template <class T>
template <class T1,class T2,class T3>
void Image<T>::Multiply(const Image<T1>& image1,const Image<T2>& image2,const Image<T3>& image3)
{
	if(image1.matchDimension(image2)==false || image2.matchDimension(image3)==false)
	{
		std::cerr<<"Error in image dimensions--function Image<T>::Multiply()!"<<std::endl;
		return;
	}
	if(matchDimension(image1)==false)
		allocate(image1);

	const T1*& pData1=image1.data();
	const T2*& pData2=image2.data();
	const T3*& pData3=image3.data();

#ifdef WITH_SSE
	if (typeid(T) == typeid(float) && typeid(T1) == typeid(float)
		&& typeid(T2) == typeid(float) && typeid(T3) == typeid(float)){
		__m128 *p0 = (__m128*)pData, *p1 = (__m128*)pData1, 
			*p2 = (__m128*)pData2, *p3 = (__m128*)pData3;
		__m128 tmp;
		for (int i = 0; i < nElements / 4; i++){
			tmp = _mm_mul_ps(*p1, *p2);
			*p0 = _mm_mul_ps(tmp, *p3);
			p0++; p1++; p2++; p3++;
		}
	}else
#endif
	{
		for(int i=0;i<nElements;i++)
			pData[i] = pData1[i] * pData2[i] * pData3[i];
	}
}

template <class T>
template <class T1,class T2>
void Image<T>::Multiply(const Image<T1>& image1,const Image<T2>& image2)
{
	if(image1.matchDimension(image2)==false)
	{
		std::cerr<<"Error in image dimensions--function Image<T>::Multiply()!"<<std::endl;
		return;
	}
	if(matchDimension(image1)==false)
		allocate(image1);

	const T1*& pData1=image1.data();
	const T2*& pData2=image2.data();

	for(int i=0;i<nElements;i++)
		pData[i]=pData1[i]*pData2[i];
}

template <class T>
template <class T1,class T2>
void Image<T>::MultiplyAcross(const Image<T1>& image1,const Image<T2>& image2)
{
	if(image1.width() != image2.width() || image1.height()!=image2.height() || image2.nchannels()!=1)
	{
		std::cerr<<"Error in image dimensions--function Image<T>::Multiply()!"<<std::endl;
		return;
	}
	if(matchDimension(image1)==false)
		allocate(image1);

	const T1*& pData1=image1.data();
	const T2*& pData2=image2.data();

	for(int i = 0;i<nPixels;i++)
		for(int j=0;j<nChannels;j++)
			pData[i*nChannels+j] = pData1[i*nChannels+j]*pData2[i];
}

template <class T>
template <class T1>
void Image<T>::Multiplywith(const Image<T1> &image1)
{
	if(matchDimension(image1)==false)
	{
		std::cerr<<"Error in image dimensions--function Image<T>::Multiplywith()!"<<std::endl;
		return;
	}
	const T1*& pData1=image1.data();
	for(int i=0;i<nElements;i++)
		pData[i]*=pData1[i];
}

template <class T>
template <class T1>
void Image<T>::MultiplywithAcross(const Image<T1> &image1)
{
	if(imWidth!=image1.width() || imHeight!=image1.height() || image1.nchannels()!=1)
	{
		std::cerr<<"Error in image dimensions--function Image<T>::MultiplywithAcross()!"<<std::endl;
		return;
	}
	const T1*& pData1=image1.data();
	for(int i=0;i<nPixels;i++)
		for(int j = 0;j<nChannels;j++)
			pData[i*nChannels+j]*=pData1[i];
}


template <class T>
void Image<T>::Multiplywith(float value)
{
#ifdef WITH_SSE
	if (typeid(T) == typeid(float)){
		__m128 *p0 = (__m128*)pData;
		__m128 _val;
		_val = _mm_set_ps1(value);
		for (int i = 0; i < nElements / 4; i++){
			*p0 = _mm_mul_ps(*p0, _val);
			p0++;
		}
	}
	else
#endif
	{
		for (int i = 0; i < nElements; i++)
			pData[i] *= value;
	}
}

//------------------------------------------------------------------------------------------
// function to add image2 to image1 to the current image
//------------------------------------------------------------------------------------------
template <class T>
template <class T1,class T2>
void Image<T>::Add(const Image<T1>& image1,const Image<T2>& image2)
{
	if(image1.matchDimension(image2)==false)
	{
		std::cerr<<"Error in image dimensions--function Image<T>::Add()!"<<std::endl;
		return;
	}
	if(matchDimension(image1)==false)
		allocate(image1);

	const T1*& pData1=image1.data();
	const T2*& pData2=image2.data();

#ifdef WITH_SSE
	if (typeid(T) == typeid(float) && typeid(T1) == typeid(float) && typeid(T2) == typeid(float)){
		__m128 *p0 = (__m128*)pData, *p1 = (__m128*)pData1, *p2 = (__m128*)pData2;
		for (int i = 0; i < nElements / 4; i++){
			*p0 = _mm_add_ps(*p1, *p2);
			p0++; p1++; p2++;
		}
	}
	else
#endif
	{
		for (int i = 0; i < nElements; i++)
			pData[i] = pData1[i] + pData2[i];
	}
}

template <class T>
template <class T1,class T2>
void Image<T>::Add(const Image<T1>& image1,const Image<T2>& image2,float ratio)
{
	if(image1.matchDimension(image2)==false)
	{
		std::cerr<<"Error in image dimensions--function Image<T>::Add()!"<<std::endl;
		return;
	}
	if(matchDimension(image1)==false)
		allocate(image1);

	const T1*& pData1=image1.data();
	const T2*& pData2=image2.data();

#ifdef WITH_SSE
	if (typeid(T) == typeid(float) && typeid(T1) == typeid(float) && typeid(T2) == typeid(float)){
		__m128 *p0 = (__m128*)pData, *p1 = (__m128*)pData1, *p2 = (__m128*)pData2;
		__m128 _rat, _tmp;
		_rat = _mm_set_ps1(ratio);
		for (int i = 0; i < nElements / 4; i++){
			_tmp = _mm_mul_ps(*p2, _rat);
			*p0 = _mm_add_ps(*p1, _tmp);
			p0++; p1++; p2++;
		}
	}
	else
#endif
	{
		for (int i = 0; i < nElements; i++)
			pData[i] = pData1[i] + pData2[i] * ratio;
	}
}

template <class T>
template <class T1>
void Image<T>::Add(const Image<T1>& image1,const float ratio)
{
	if(matchDimension(image1)==false)
	{
		std::cerr<<"Error in image dimensions--function Image<T>::Add()!"<<std::endl;
		return;
	}
	const T1*& pData1=image1.data();

#ifdef WITH_SSE
	if (typeid(T) == typeid(float) && typeid(T1) == typeid(float)){
		__m128 *p0 = (__m128*)pData, *p1 = (__m128*)pData1;
		__m128 _rat, _tmp;
		_rat = _mm_set_ps1(ratio);
		for (int i = 0; i < nElements / 4; i++){
			_tmp = _mm_mul_ps(*p1, _rat);
			*p0 = _mm_add_ps(*p0, _tmp);
			p0++; p1++;
		}
	}
	else
#endif
	{
		for (int i = 0; i < nElements; i++)
			pData[i] += pData1[i] * ratio;
	}
}

template <class T>
template <class T1>
void Image<T>::Add(const Image<T1>& image1)
{
	if(matchDimension(image1)==false)
	{
		std::cerr<<"Error in image dimensions--function Image<T>::Add()!"<<std::endl;
		return;
	}
	const T1*& pData1=image1.data();

#ifdef WITH_SSE
	if (typeid(T) == typeid(float) && typeid(T1) == typeid(float)){
		__m128 *p0 = (__m128*)pData, *p1 = (__m128*)pData1;
		for (int i = 0; i < nElements / 4; i++){
			*p0 = _mm_add_ps(*p0, *p1);
			p0++; p1++;
		}
	}
	else
#endif
	{
		for (int i = 0; i < nElements; i++)
			pData[i] += pData1[i];
	}
}

template <class T>
void Image<T>::Add(const T value)
{
#ifdef WITH_SSE
	if (typeid(T) == typeid(float)){
		__m128 *p0 = (__m128*)pData;
		__m128 _val;
		_val = _mm_set_ps1(value);
		for (int i = 0; i < nElements / 4; i++){
			*p0 = _mm_add_ps(*p0, _val);
			p0++;
		}
	}
	else
#endif
	{
		for (int i = 0; i < nElements; i++)
			pData[i] += value;
	}
}

//------------------------------------------------------------------------------------------
// function to subtract image2 from image1
//------------------------------------------------------------------------------------------
template <class T>
template <class T1,class T2>
void Image<T>::Subtract(const Image<T1> &image1, const Image<T2> &image2)
{
	if(image1.matchDimension(image2)==false)
	{
		std::cerr<<"Error in image dimensions--function Image<T>::Subtract()!"<<std::endl;
		return;
	}
	if(matchDimension(image1)==false)
		allocate(image1);

	const T1*& pData1=image1.data();
	const T2*& pData2=image2.data();

#ifdef WITH_SSE
	if (typeid(T) == typeid(float) && typeid(T1) == typeid(float) && typeid(T2) == typeid(float)){
		__m128 *p0 = (__m128*)pData, *p1 = (__m128*)pData1, *p2 = (__m128*)pData2;
		for (int i = 0; i < nElements / 4; i++){
			*p0 = _mm_sub_ps(*p1, *p2);
			p0++; p1++; p2++;
		}
	}
	else
#endif
	{
		for (int i = 0; i < nElements; i++)
			pData[i] = (T)pData1[i] - pData2[i];
	}
}

template <class T>
void Image<T>::Subtract(const T value)
{
#ifdef WITH_SSE
	if (typeid(T) == typeid(float)){
		__m128 *p0 = (__m128*)pData;
		__m128 _val;
		_val = _mm_set_ps1(value);
		for (int i = 0; i < nElements / 4; i++){
			*p0 = _mm_sub_ps(*p0, _val);
			p0++;
		}
	}
	else
#endif
	{
		for (int i = 0; i < nElements; i++)
			pData[i] -= value;
	}
}

//------------------------------------------------------------------------------------------
// square an image
//------------------------------------------------------------------------------------------
template <class T>
void Image<T>::square()
{
	for(int i = 0;i<nElements;i++)
		pData[i] = pData[i]*pData[i];
}


//------------------------------------------------------------------------------------------
// exp an image
//------------------------------------------------------------------------------------------
template <class T>
void Image<T>::Exp(float sigma)
{
	for(int i = 0;i<nElements;i++)
		pData[i] = exp(-pData[i]/sigma);
}

//------------------------------------------------------------------------------------------
// normalize an image
//------------------------------------------------------------------------------------------
template <class T>
template <class T1>
void Image<T>::normalize(Image<T1>& image, float minV/* = -1*/, float maxV/* = -1*/) const
{
	if(image.width()!=imWidth || image.height()!=imHeight || image.nchannels()!=nChannels)
		image.allocate(imWidth,imHeight,nChannels);

	T realMin = immin();
	T realMax = immax();

	if (minV < 0){
		minV = 0;
	}
	if (maxV < 0){
		if (image.IsFloat())
			maxV = 1;
		else
			maxV = 255;
	}

	if (realMax == realMin){
		image.copy(*this);
		return;
	}
	float ratio = (maxV - minV) / (realMax - realMin);
	T1* data=image.data();
	for (int i = 0; i < nElements; i++)
		data[i] = (float)(pData[i] - realMin)*ratio + minV;
}

template <class T>
void Image<T>::normalize(float minV/* = -1*/, float maxV/* = -1*/)
{
	Image<T> temp;
	normalize(temp, minV, maxV);
	copyData(temp);
}

template <class T>
float Image<T>::norm2() const
{
	float result=0;
	for(int i=0;i<nElements;i++)
		result+=pData[i]*pData[i];
	return result;
}

template <class T>
void Image<T>::threshold(float minV/* = std::numeric_limits<float>::min()*/, float maxV/* = std::numeric_limits<float>::max()*/)
{
	if (minV == std::numeric_limits<float>::min()){
		minV = 0;
	}
	if (maxV == std::numeric_limits<float>::max()){
		if (IsFloat())
			maxV = 1;
		else
			maxV = 255;
	}
	for(int i = 0;i<nPixels*nChannels;i++)
		pData[i] = std::min(std::max(pData[i], minV), maxV);
}

template <class T>
float Image<T>::sum() const
{
	float result=0;
	for(int i=0;i<nElements;i++)
		result+=pData[i];
	return result;
}

template <class T>
template <class T1>
float Image<T>::innerproduct(Image<T1> &image) const
{
	float result=0;
	const T1* pData1=image.data();
	for(int i=0;i<nElements;i++)
		result+=pData[i]*pData1[i];
	return result;
}

#ifndef _MATLAB
//template <class T>
//bool Image<T>::writeImage(QFile &file) const
//{
//	file.write((char *)&imWidth,sizeof(int));
//	file.write((char *)&imHeight,sizeof(int));
//	file.write((char *)&nChannels,sizeof(int));
//	file.write((char *)&IsDerivativeImage,sizeof(bool));
//	if(file.write((char *)pData,sizeof(T)*nElements)!=sizeof(T)*nElements)
//		return false;
//	return true;
//}
//
//template <class T>
//bool Image<T>::readImage(QFile& file)
//{
//	clear();
//	file.read((char *)&imWidth,sizeof(int));
//	file.read((char *)&imHeight,sizeof(int));
//	file.read((char *)&nChannels,sizeof(int));
//	file.read((char *)&IsDerivativeImage,sizeof(bool));
//	if(imWidth<0 ||imWidth>100000 || imHeight<0 || imHeight>100000 || nChannels<0 || nChannels>10000)
//		return false;
//	allocate(imWidth,imHeight,nChannels);
//	if(file.read((char *)pData,sizeof(T)*nElements)!=sizeof(T)*nElements)
//		return false;
//	return true;
//}
//
//template <class T>
//bool Image<T>::writeImage(const QString &filename) const
//{
//	QFile file(filename);
//	if(file.open(QIODevice::WriteOnly)==false)
//		return false;
//	if(!writeImage(file))
//		return false;
//	return true;
//}
//
//template <class T>
//bool Image<T>::readImage(const QString &filename)
//{
//	QFile file(filename);
//	if(file.open(QIODevice::ReadOnly)==false)
//		return false;
//	if(!readImage(file))
//		return false;
//	return true;
//}
#endif

template <class T>
template <class T1>
void Image<T>::Integral(Image<T1>& image) const
{
	ImageProcessing::Integral(pData, image.data(), imWidth, imHeight, nChannels);
}

template <class T>
template <class T1>
void Image<T>::BoxFilter(Image<T1>& image, int r, bool norm /*= true*/) const
{
	ImageProcessing::BoxFilter(pData, image.data(), imWidth, imHeight, nChannels, r, norm);
}

template <class T>
template <class T1>
void Image<T>::BilateralFiltering(Image<T1>& other,int fsize,float filter_sigma,float range_sigma) const
{
	float *pBuffer;
	Image<T1> result(other);
	pBuffer=new float[other.nchannels()];

	// set spatial weight to save time
	float *pSpatialWeight;
	int flength = fsize*2+1;
	pSpatialWeight = new float[flength*flength];
	for(int i=-fsize;i<=fsize;i++)
		for(int j=-fsize;j<=fsize;j++)
			pSpatialWeight[(i+fsize)*flength+j+fsize]=exp(-(float)(i*i+j*j)/(2*filter_sigma*filter_sigma));

	for(int i=0;i<imHeight;i++)
		for(int j=0;j<imWidth;j++)
		{
			float totalWeight=0;
			for(int k=0;k<other.nchannels();k++)
				pBuffer[k]=0;
			for(int ii=-fsize;ii<=fsize;ii++)
				for(int jj=-fsize;jj<=fsize;jj++)
				{
					int x=j+jj;
					int y=i+ii;
					if(x<0 || x>=imWidth || y<0 || y>=imHeight)
						continue;

					// compute weight
					int offset=(y*imWidth+x)*nChannels;
					float temp=0;
					for(int k=0;k<nChannels;k++)
					{
						float diff=pData[offset+k]-pData[(i*imWidth+j)*nChannels+k];
						temp+=diff*diff;
					}
					float weight=exp(-temp/(2*range_sigma*range_sigma));
					weight *= pSpatialWeight[(ii+fsize)*flength+jj+fsize];
					//weight*=exp(-(float)(ii*ii+jj*jj)/(2*filter_sigma*filter_sigma));
					totalWeight+=weight;
					for(int k=0;k<other.nchannels();k++)
						pBuffer[k]+=other.data()[(y*imWidth+x)*other.nchannels()+k]*weight;
				}
			for(int k=0;k<other.nchannels();k++)
				result.data()[(i*imWidth+j)*other.nchannels()]=pBuffer[k]/totalWeight;
		}
	other.copyData(result);
	delete []pBuffer;
	delete []pSpatialWeight;
}


template <class T>
//Image<T>  Image<T>::BilateralFiltering(int fsize,float filter_sigma,float range_sigma)
void  Image<T>::imBilateralFiltering(Image<T>& result,int fsize,float filter_sigma,float range_sigma) const
{
	//Image<T> result(*this);
	result.allocate(*this);

	float *pBuffer;
	pBuffer=new float[nChannels];

	// set spatial weight to save time
	float *pSpatialWeight;
	int flength = fsize*2+1;
	pSpatialWeight = new float[flength*flength];
	for(int i=-fsize;i<=fsize;i++)
		for(int j=-fsize;j<=fsize;j++)
			pSpatialWeight[(i+fsize)*flength+j+fsize]=exp(-(float)(i*i+j*j)/(2*filter_sigma*filter_sigma));

	for(int i=0;i<imHeight;i++)
		for(int j=0;j<imWidth;j++)
		{
			float totalWeight=0;
			for(int k=0;k<nChannels;k++)
				pBuffer[k]=0;
			int offset0 = (i*imWidth+j)*nChannels;
			for(int ii=-fsize;ii<=fsize;ii++)
				for(int jj=-fsize;jj<=fsize;jj++)
				{
					int x=j+jj;
					int y=i+ii;
					if(x<0 || x>=imWidth || y<0 || y>=imHeight)
						continue;

					// compute weight
					int offset=(y*imWidth+x)*nChannels;
					float temp=0;
					for(int k=0;k<nChannels;k++)
					{
						float diff=pData[offset+k]-pData[offset0+k];
						temp+=diff*diff;
					}
					float weight=exp(-temp/(2*range_sigma*range_sigma));
					weight *= pSpatialWeight[(ii+fsize)*flength+jj+fsize];

					//weight*=exp(-(float)(ii*ii+jj*jj)/(2*filter_sigma*filter_sigma));
					totalWeight+=weight;
					for(int k=0;k<nChannels;k++)
						pBuffer[k]+=pData[offset+k]*weight;
				}
			for(int k=0;k<nChannels;k++)
				result.data()[offset0+k]=pBuffer[k]/totalWeight;

		}
	delete []pBuffer;
	delete []pSpatialWeight;
	//return result;
}

template <class T>
template <class T1,class T2>
int Image<T>::kmeansIndex(int pixelIndex,T1& MinDistance,const T2* pDictionary,int nVocabulary,int nDim)
{
	int offset1 = pixelIndex*nChannels;
	T1 Distance = 0;


	int index;
	for(int j = 0;j<nVocabulary;j++)
	{
		int offset2 = j*nDim;
		Distance = 0;
		for(int k = 0;k<nDim;k++)
			Distance += (T1)(pData[offset1+k] - pDictionary[offset2+k])*(pData[offset1+k] - pDictionary[offset2+k]);
		if(j==0)
		{
			MinDistance = Distance;
			index = 0;
		}
		else if(Distance < MinDistance)
		{
			MinDistance = Distance;
			index = j;
		}
	}
	return index;
}

// function to convert an image to visual words based on the vocabulary
template <class T>
template <class T1,class T2>
void Image<T>::ConvertToVisualWords(Image<T1> &result, const T2 *pDictionary, int nDim, int nVocabulary)
{
	if(nChannels !=nDim)
	{
		std::cerr<<"The dimension of the vocabulary must match to the nChannels of the image"<<std::endl;
		return;
	}
	if(result.matchDimension(imWidth,imHeight,1))
		result.allocate(imWidth,imHeight);
	
	bool isFloat = IsFloat();
	for(int i = 0;i<nPixels;i++)
	{
		if(isFloat)
		{
			float minDistance;
			result[i] = kmeansIndex(i,minDistance,pDictionary,nVocabulary);
		}
		else
		{
			int minDistance;
			result[i] = kmeansIndex(i,minDistance,pDictionary,nVocabulary);
		}
	}
}

// function to count the histogram of a specified region
// notice that the range of the coordinate is [0,imWidth] (x) and [0,imHeight] (y)
template <class T>
template <class T1>
Vector<T1> Image<T>::histogramRegion(int nBins,float left,float top,float right,float bottom) const
{
	Vector<T1> histogram(nBins);
	int Left = left,Top = top,Right = right,Bottom = bottom;
	float dLeft,dTop,dRight,dBottom;
	dLeft = 1-(left-Left);
	dTop = 1-(top-Top);
	if(right > Right)
	{
		dRight = right - Right;
	}
	else
	{
		Right --;
		dRight = 1;
	}
	if(bottom > Bottom)
	{
		dBottom = bottom - Bottom;
	}
	else
	{
		Bottom --;
		dBottom = 1;
	}

	for(int i = Top; i <= Bottom; i++)
		for(int j = Left; j <= Right; j++)
		{
			int offset = (i*imWidth+j)*nChannels;
			float weight=1;

			if(Top==Bottom)
				weight *= (dTop+dBottom-1);
			else
			{
				if(i==Top)
					weight *= dTop;
				if(i==Bottom)
					weight *= dBottom;
			}			
			
			if(Left==Right)
				weight *= (dLeft+dRight-1);
			else
			{
				if(j==Left)
					weight *= dLeft;
				if(j==Right)
					weight *= dRight;
			}

			for(int k = 0;k<nChannels;k++)
				histogram[pData[offset+k]] += weight;
		}
	return histogram;
}

//-----------------------------------------------------------------------------------------------------------------------
// functions for bicubic interpolation
//-----------------------------------------------------------------------------------------------------------------------
template <class T>
template <class T1>
void Image<T>::warpImageBicubic(Image<T>& output,const Image<T1>& vx,const Image<T1>& vy) const
{
	float dfilter[3] = {-0.5,0,0.5};
	FImage imdx,imdy,imdxdy;
	imfilter_h(imdx,dfilter,1);
	imfilter_v(imdy,dfilter,1);
	imdx.imfilter_v(imdxdy,dfilter,1);
	warpImageBicubic(output,imdx,imdy,imdxdy,vx,vy);
}

template <class T>
template <class T1,class T2>
void Image<T>::warpImageBicubic(Image<T>& output,const Image<T1>& imdx,const Image<T1>& imdy,const Image<T1>& imdxdy,
																		const Image<T2>& vx,const Image<T2>& vy) const
{
	T* pIm = pData;
	const T1* pImDx = imdx.data();
	const T1* pImDy = imdy.data();
	const T1* pImDxDy = imdxdy.data();
	int width = vx.width();
	int height = vx.height();
	if(!output.matchDimension(width,height,nChannels))
		output.allocate(width,height,nChannels);
	float a[4][4];
	int offsets[2][2];

	T ImgMax;
	if(IsFloat())
		ImgMax = 1;
	else
		ImgMax = 255;

	for(int i  = 0; i<height; i++)
		for(int j = 0;j<width;j++)
		{
			int offset = i*width+j;
			float x = j + vx.pData[offset];
			float y = i + vy.pData[offset];
			int x0 = x;
			int y0 = y;
			int x1 = x0+1;
			int y1 = y0+1;
			x0 = std::min(std::max(x0,0),imWidth-1);
			x1 = std::min(std::max(x1,0),imWidth-1);
			y0 = std::min(std::max(y0,0),imHeight-1);
			y1 = std::min(std::max(y1,0),imHeight-1);

			float dx = x - x0;
			float dy = y- y0;
			float dx2 = dx*dx;
			float dy2 = dy*dy;
			float dx3 = dx*dx2;
			float dy3 = dy*dy2;


			for(int k = 0;k<nChannels;k++)
			{
				offsets[0][0] = (y0*imWidth+x0)*nChannels + k;
				offsets[1][0] = (y0*imWidth+x1)*nChannels + k;
				offsets[0][1] = (y1*imWidth+x0)*nChannels + k;
				offsets[1][1] = (y1*imWidth+x1)*nChannels + k;

				// set the sampling coefficients
				BicubicCoeff(a,pIm,pImDx,pImDy,pImDxDy,offsets);

				// now use the coefficients for interpolation
				output.pData[offset*nChannels+k] = a[0][0] +          a[0][1]*dy +          a[0][2]*dy2 +           a[0][3]*dy3+
					                                                                    a[1][0]*dx +   a[1][1]*dx*dy   + a[1][2]*dx*dy2   + a[1][3]*dx*dy3 + 
																						a[2][0]*dx2 + a[2][1]*dx2*dy + a[2][2]*dx2*dy2 + a[2][3]*dx2*dy3+
																						a[3][0]*dx3 + a[3][1]*dx3*dy + a[3][2]*dx3*dy2 + a[3][3]*dx3*dy3;
				//output.pData[offset*nChannels+k] = std::max(std::min(output.pData[offset*nChannels+k],ImgMax),0);

			}
		}
}


template <class T>
template <class T1,class T2>
void Image<T>::warpImageBicubic(Image<T>& output,const Image<T1>& coeff,const Image<T2>& vx,const Image<T2>& vy) const
{
	T* pIm = pData;
	int width = vx.width();
	int height = vx.height();
	if(!output.matchDimension(width,height,nChannels))
		output.allocate(width,height,nChannels);
	float a[4][4];

	T ImgMax;
	if(IsFloat())
		ImgMax = 1;
	else
		ImgMax = 255;

	for(int i  = 0; i<height; i++)
		for(int j = 0;j<width;j++)
		{
			int offset = i*width+j;
			float x = j + vx.pData[offset];
			float y = i + vy.pData[offset];
			int x0 = x;
			int y0 = y;
			int x1 = x0+1;
			int y1 = y0+1;
			x0 = std::min(std::max(x0,0),imWidth-1);
			x1 = std::min(std::max(x1,0),imWidth-1);
			y0 = std::min(std::max(y0,0),imHeight-1);
			y1 = std::min(std::max(y1,0),imHeight-1);

			float dx = x - x0;
			float dy = y- y0;
			float dx2 = dx*dx;
			float dy2 = dy*dy;
			float dx3 = dx*dx2;
			float dy3 = dy*dy2;


			for(int k = 0;k<nChannels;k++)
			{
				// save the coefficients
				for(int ii = 0;ii<4;ii++)
					for(int jj=0;jj<4;jj++)
						a[ii][jj] = coeff.pData[(offset*nChannels+k)*16+ii*4+jj];


				// set the sampling coefficients

				// now use the coefficients for interpolation
				output.pData[offset*nChannels+k] = a[0][0] +          a[0][1]*dy +          a[0][2]*dy2 +           a[0][3]*dy3+
					                                                                    a[1][0]*dx +   a[1][1]*dx*dy   + a[1][2]*dx*dy2   + a[1][3]*dx*dy3 + 
																						a[2][0]*dx2 + a[2][1]*dx2*dy + a[2][2]*dx2*dy2 + a[2][3]*dx2*dy3+
																						a[3][0]*dx3 + a[3][1]*dx3*dy + a[3][2]*dx3*dy2 + a[3][3]*dx3*dy3;
				//output.pData[offset*nChannels+k] = std::max(std::min(output.pData[offset*nChannels+k],ImgMax),0);

			}
		}
}

template <class T>
template <class T1>
void Image<T>::BicubicCoeff(float a[][4],const T* pIm,const T1* pImDx,const T1* pImDy,const T1* pImDxDy,const int offsets[][2]) const
{
		a[0][0] = pIm[offsets[0][0]];
		a[1][0] = pImDx[offsets[0][0]];
		a[2][0] = -3*pIm[offsets[0][0]] + 3*pIm[offsets[1][0]] -2*pImDx[offsets[0][0]] - pImDx[offsets[1][0]];
		a[3][0] =   2*pIm[offsets[0][0]] -  2*pIm[offsets[1][0]] +   pImDx[offsets[0][0]] +pImDx[offsets[1][0]];

		a[0][1] = pImDy[offsets[0][0]];
		a[1][1] = pImDxDy[offsets[0][0]];
		a[2][1] = -3*pImDy[offsets[0][0]] + 3*pImDy[offsets[1][0]] - 2*pImDxDy[offsets[0][0]] - pImDxDy[offsets[1][0]];
		a[3][1] = 2*pImDy[offsets[0][0]] - 2*pImDy[offsets[1][0]] + pImDxDy[offsets[0][0]] + pImDxDy[offsets[1][0]];

		a[0][2] =      -3*pIm[offsets[0][0]]      + 3*pIm[offsets[0][1]]       -2*pImDy[offsets[0][0]]        - pImDy[offsets[0][1]];
		a[1][2] = -3*pImDx[offsets[0][0]] + 3*pImDx[offsets[0][1]] -2*pImDxDy[offsets[0][0]] - pImDxDy[offsets[0][1]];
		a[2][2] =		     9*pIm[offsets[0][0]]      -        9*pIm[offsets[1][0]]     -        9*pIm[offsets[0][1]]     +    9*pIm[offsets[1][1]] +
								6*pImDx[offsets[0][0]]   +    3*pImDx[offsets[1][0]]   -     6*pImDx[offsets[0][1]] -    3*pImDx[offsets[1][1]] +
								6*pImDy[offsets[0][0]]   -     6*pImDy[offsets[1][0]] +      3*pImDy[offsets[0][1]] -    3*pImDy[offsets[1][1]] +
							4*pImDxDy[offsets[0][0]] + 2*pImDxDy[offsets[1][0]] + 2*pImDxDy[offsets[0][1]] + pImDxDy[offsets[1][1]];
		a[3][2] =		    -6*pIm[offsets[0][0]]      +      6*pIm[offsets[1][0]]     +       6*pIm[offsets[0][1]]     -     6*pIm[offsets[1][1]] +
							(-3)*pImDx[offsets[0][0]]   -     3*pImDx[offsets[1][0]]   +    3*pImDx[offsets[0][1]] +   3*pImDx[offsets[1][1]] +
							(-4)*pImDy[offsets[0][0]]   +    4*pImDy[offsets[1][0]]    -    2*pImDy[offsets[0][1]] +   2*pImDy[offsets[1][1]] +
						(-2)*pImDxDy[offsets[0][0]]  - 2*pImDxDy[offsets[1][0]]   -    pImDxDy[offsets[0][1]]   -  pImDxDy[offsets[1][1]];

		a[0][3] =      2*pIm[offsets[0][0]]        - 2*pIm[offsets[0][1]]       + pImDy[offsets[0][0]]        + pImDy[offsets[0][1]];
		a[1][3] = 2*pImDx[offsets[0][0]]  - 2*pImDx[offsets[0][1]]  + pImDxDy[offsets[0][0]] + pImDxDy[offsets[0][1]];
		a[2][3] =		    -6*pIm[offsets[0][0]]      +      6*pIm[offsets[1][0]]     +       6*pIm[offsets[0][1]]     -     6*pIm[offsets[1][1]] +
							(-4)*pImDx[offsets[0][0]]   -     2*pImDx[offsets[1][0]]   +    4*pImDx[offsets[0][1]] +   2*pImDx[offsets[1][1]] +
							(-3)*pImDy[offsets[0][0]]   +    3*pImDy[offsets[1][0]]    -    3*pImDy[offsets[0][1]] +   3*pImDy[offsets[1][1]] +
						(-2)*pImDxDy[offsets[0][0]]  -     pImDxDy[offsets[1][0]] -  2*pImDxDy[offsets[0][1]]   -  pImDxDy[offsets[1][1]];
		a[3][3] =		     4*pIm[offsets[0][0]]      -        4*pIm[offsets[1][0]]     -        4*pIm[offsets[0][1]]     +    4*pIm[offsets[1][1]] +
								2*pImDx[offsets[0][0]]   +    2*pImDx[offsets[1][0]]   -     2*pImDx[offsets[0][1]] -    2*pImDx[offsets[1][1]] +
								2*pImDy[offsets[0][0]]   -     2*pImDy[offsets[1][0]] +      2*pImDy[offsets[0][1]] -    2*pImDy[offsets[1][1]] +
								pImDxDy[offsets[0][0]] +     pImDxDy[offsets[1][0]] +      pImDxDy[offsets[0][1]] + pImDxDy[offsets[1][1]];
}

template <class T>
template <class T1>
void Image<T>::warpImageBicubicCoeff(Image<T1>& output) const
{
	// generate derivatie filters
	Image<float> imdx,imdy,imdxdy;
	float dfilter[3] = {-0.5,0,0.5};
	imfilter_h(imdx,dfilter,1);
	imfilter_v(imdy,dfilter,1);
	imdx.imfilter_v(imdxdy,dfilter,1);

	T* pIm = pData;
	const T1* pImDx = imdx.data();
	const T1* pImDy = imdy.data();
	const T1* pImDxDy = imdxdy.data();

	if(!output.matchDimension(imWidth,imHeight,nChannels*16))
		output.allocate(imWidth,imHeight,nChannels*16);
	float a[4][4];
	int offsets[2][2];

	for(int i  = 0; i<imHeight; i++)
		for(int j = 0;j<imWidth;j++)
		{
			int offset = i*imWidth+j;
			int x0 = j;
			int y0 = i;
			int x1 = x0+1;
			int y1 = y0+1;
			x0 = std::min(std::max(x0,0),imWidth-1);
			x1 = std::min(std::max(x1,0),imWidth-1);
			y0 = std::min(std::max(y0,0),imHeight-1);
			y1 = std::min(std::max(y1,0),imHeight-1);

			for(int k = 0;k<nChannels;k++)
			{
				offsets[0][0] = (y0*imWidth+x0)*nChannels + k;
				offsets[1][0] = (y0*imWidth+x1)*nChannels + k;
				offsets[0][1] = (y1*imWidth+x0)*nChannels + k;
				offsets[1][1] = (y1*imWidth+x1)*nChannels + k;

				// set the sampling coefficients
				BicubicCoeff(a,pIm,pImDx,pImDy,pImDxDy,offsets);

				// save the coefficients
				for(int ii = 0;ii<4;ii++)
					for(int jj=0;jj<4;jj++)
						output.pData[(offset*nChannels+k)*16+ii*4+jj] = a[ii][jj];
			}
		}
}

template <class T>
template <class T1>
void Image<T>::warpImageBicubicRef(const Image<T>& ref,Image<T>& output,const Image<T1>& vx,const Image<T1>& vy) const
{
	float dfilter[3] = {-0.5,0,0.5};
	FImage imdx,imdy,imdxdy;
	imfilter_h(imdx,dfilter,1);
	imfilter_v(imdy,dfilter,1);
	imdx.imfilter_v(imdxdy,dfilter,1);
	warpImageBicubicRef(ref,output,imdx,imdy,imdxdy,vx,vy);
}

template <class T>
template <class T1>
void Image<T>::warpImageBicubicRef(const Image<T>& ref,Image<T>& output,const Image<T1>& flow) const
{
	FImage vx,vy;
	flow.DissembleFlow(vx,vy);
	warpImageBicubicRef(ref,output,vx,vy);
}

template <class T>
template <class T1>
void Image<T>::DissembleFlow(Image<T1>& vx,Image<T1>& vy) const
{
	if(!vx.matchDimension(imWidth,imHeight,1))
		vx.allocate(imWidth,imHeight,1);
	if(!vy.matchDimension(imWidth,imHeight,1))
		vy.allocate(imWidth,imHeight,1);
	for(int i =0;i<vx.npixels();i++)
	{
		vx.data()[i] = pData[i*2];
		vy.data()[i] = pData[i*2+1];
	}
}

template <class T>
template <class T1,class T2>
void Image<T>::warpImageBicubicRef(const Image<T>& ref,Image<T>& output,const Image<T1>& imdx,const Image<T1>& imdy,const Image<T1>& imdxdy,
																		const Image<T2>& vx,const Image<T2>& vy) const
{
	T* pIm = pData;
	const T1* pImDx = imdx.data();
	const T1* pImDy = imdy.data();
	const T1* pImDxDy = imdxdy.data();
	int width = vx.width();
	int height = vx.height();
	if(!output.matchDimension(width,height,nChannels))
		output.allocate(width,height,nChannels);
	float a[4][4];
	int offsets[2][2];

	T ImgMax;
	if(IsFloat())
		ImgMax = 1;
	else
		ImgMax = 255;

	for(int i  = 0; i<height; i++)
		for(int j = 0;j<width;j++)
		{
			int offset = i*width+j;
			float x = j + vx.pData[offset];
			float y = i + vy.pData[offset];
			if(x<0 || x>imWidth-1 || y<0 || y>imHeight-1)
			{
				for(int k = 0; k<nChannels;k++)
					output.pData[offset*nChannels+k] = ref.pData[offset*nChannels+k];
				continue;
			}
			int x0 = x;
			int y0 = y;
			int x1 = x0+1;
			int y1 = y0+1;
			x0 = std::min(std::max(x0,0),imWidth-1);
			x1 = std::min(std::max(x1,0),imWidth-1);
			y0 = std::min(std::max(y0,0),imHeight-1);
			y1 = std::min(std::max(y1,0),imHeight-1);

			float dx = x - x0;
			float dy = y- y0;
			float dx2 = dx*dx;
			float dy2 = dy*dy;
			float dx3 = dx*dx2;
			float dy3 = dy*dy2;


			for(int k = 0;k<nChannels;k++)
			{
				offsets[0][0] = (y0*imWidth+x0)*nChannels + k;
				offsets[1][0] = (y0*imWidth+x1)*nChannels + k;
				offsets[0][1] = (y1*imWidth+x0)*nChannels + k;
				offsets[1][1] = (y1*imWidth+x1)*nChannels + k;

				// set the sampling coefficients
				BicubicCoeff(a,pIm,pImDx,pImDy,pImDxDy,offsets);

				// now use the coefficients for interpolation
				output.pData[offset*nChannels+k] = a[0][0] +          a[0][1]*dy +          a[0][2]*dy2 +           a[0][3]*dy3+
					                                                                    a[1][0]*dx +   a[1][1]*dx*dy   + a[1][2]*dx*dy2   + a[1][3]*dx*dy3 + 
																						a[2][0]*dx2 + a[2][1]*dx2*dy + a[2][2]*dx2*dy2 + a[2][3]*dx2*dy3+
																						a[3][0]*dx3 + a[3][1]*dx3*dy + a[3][2]*dx3*dy2 + a[3][3]*dx3*dy3;
				//output.pData[offset*nChannels+k] = std::max(std::min(output.pData[offset*nChannels+k],ImgMax),0);
				//if(!(output.pData[offset*nChannels+k]<100000 && output.pData[offset*nChannels+k]>-100000)) // bound the values
				//	output.pData[offset*nChannels+k] = ref.pData[offset*nChannels+k];

			}
		}
}

template <class T>
template <class T1,class T2>
void Image<T>::warpImageBicubicRef(const Image<T>& ref,Image<T>& output,const Image<T1>& coeff,const Image<T2>& vx,const Image<T2>& vy) const
{
	T* pIm = pData;
	int width = vx.width();
	int height = vx.height();
	if(!output.matchDimension(width,height,nChannels))
		output.allocate(width,height,nChannels);
	float a[4][4];

	T ImgMax;
	if(IsFloat())
		ImgMax = 1;
	else
		ImgMax = 255;

	for(int i  = 0; i<height; i++)
		for(int j = 0;j<width;j++)
		{
			int offset = i*width+j;
			float x = j + vx.pData[offset];
			float y = i + vy.pData[offset];
			if(x<0 || x>imWidth-1 || y<0 || y>imHeight-1)
			{
				for(int k = 0; k<nChannels;k++)
					output.pData[offset*nChannels+k] = ref.pData[offset*nChannels+k];
				continue;
			}

			int x0 = x;
			int y0 = y;
			int x1 = x0+1;
			int y1 = y0+1;
			x0 = std::min(std::max(x0,0),imWidth-1);
			x1 = std::min(std::max(x1,0),imWidth-1);
			y0 = std::min(std::max(y0,0),imHeight-1);
			y1 = std::min(std::max(y1,0),imHeight-1);

			float dx = x - x0;
			float dy = y- y0;
			float dx2 = dx*dx;
			float dy2 = dy*dy;
			float dx3 = dx*dx2;
			float dy3 = dy*dy2;


			for(int k = 0;k<nChannels;k++)
			{
				// load the coefficients
				for(int ii = 0;ii<4;ii++)
					for(int jj=0;jj<4;jj++)
						a[ii][jj] = coeff.pData[(offset*nChannels+k)*16+ii*4+jj];


				// now use the coefficients for interpolation
				output.pData[offset*nChannels+k] = a[0][0] +          a[0][1]*dy +          a[0][2]*dy2 +           a[0][3]*dy3+
					                                                                    a[1][0]*dx +   a[1][1]*dx*dy   + a[1][2]*dx*dy2   + a[1][3]*dx*dy3 + 
																						a[2][0]*dx2 + a[2][1]*dx2*dy + a[2][2]*dx2*dy2 + a[2][3]*dx2*dy3+
																						a[3][0]*dx3 + a[3][1]*dx3*dy + a[3][2]*dx3*dy2 + a[3][3]*dx3*dy3;
				//output.pData[offset*nChannels+k] = std::max(std::min(output.pData[offset*nChannels+k],ImgMax),0);
				//if(!(output.pData[offset*nChannels+k]<100000 && output.pData[offset*nChannels+k]>-100000)) // bound the values
				//	output.pData[offset*nChannels+k] = ref.pData[offset*nChannels+k];

			}
		}
}





template <class T>
template <class T1>
void Image<T>::warpImage(Image<T> &output, const Image<T1>& vx, const Image<T1>& vy) const
{
	if(!output.matchDimension(*this))
		output.allocate(*this);
	ImageProcessing::warpImage(output.data(),pData,vx.data(),vy.data(),imWidth,imHeight,nChannels);
}

template <class T>
template <class T1>
void Image<T>::warpImage_transpose(Image<T> &output, const Image<T1>& vx, const Image<T1>& vy) const
{
	if(!output.matchDimension(*this))
		output.allocate(*this);
	ImageProcessing::warpImage_transpose(output.data(),pData,vx.data(),vy.data(),imWidth,imHeight,nChannels);
}

template <class T>
template <class T1>
void Image<T>::warpImage(Image<T> &output, const Image<T1>& flow) const
{
	if(!output.matchDimension(*this))
		output.allocate(*this);
	ImageProcessing::warpImage(output.data(),pData,flow.data(),imWidth,imHeight,nChannels);
}

template <class T>
template <class T1>
void Image<T>::warpImage_transpose(Image<T> &output, const Image<T1>& flow) const
{
	if(!output.matchDimension(*this))
		output.allocate(*this);
	ImageProcessing::warpImage_transpose(output.data(),pData,flow.data(),imWidth,imHeight,nChannels);
}

template <class T>
T Image<T>::maximum() const
{
	T Max = pData[0];
	for(int i = 0;i<nElements; i++)
		Max = std::max(Max,pData[i]);
	return Max;
}

template <class T>
T Image<T>::minimum() const
{
	T Min = pData[0];
	for(int i = 0;i<nElements;i++)
		Min = std::min(Min,pData[i]);
	return Min;
}

#ifdef _MATLAB

template <class T>
template <class T1>
void Image<T>::LoadMatlabImageCore(const mxArray *image,bool IsImageScaleCovnersion)
{
	int nDim = mxGetNumberOfDimensions(image);
	const int* imDim = mxGetDimensions(image);
	if(nDim==2)
		allocate(imDim[1],imDim[0]);
	else if(nDim==3)
		allocate(imDim[1],imDim[0],imDim[2]);
	else
		mexErrMsgTxt("The image doesn't have the appropriate dimension!");
	T1* pMatlabPlane=(T1*)mxGetData(image);
	bool IsMatlabFloat;
	if(typeid(T1)==typeid(float) || typeid(T1)==typeid(float) || typeid(T1)==typeid(long float))
		IsMatlabFloat=true;
	else
		IsMatlabFloat=false;
	bool isfloat=IsFloat();
	if(isfloat==IsMatlabFloat || IsImageScaleCovnersion==false)
	{
		ConvertFromMatlab<T1>(pMatlabPlane,imWidth,imHeight,nChannels);
		return;
	}
	int offset=0;
	if(isfloat==true)
		for(int i=0;i<imHeight;i++)
			for(int j=0;j<imWidth;j++)
				for(int k=0;k<nChannels;k++)
					pData[offset++]=(float)pMatlabPlane[k*nPixels+j*imHeight+i]/255;
	else
		for(int i=0;i<imHeight;i++)
			for(int j=0;j<imWidth;j++)
				for(int k=0;k<nChannels;k++)
					pData[offset++]=(float)pMatlabPlane[k*nPixels+j*imHeight+i]*255;
}

template <class T>
bool Image<T>::LoadMatlabImage(const mxArray* matrix,bool IsImageScaleCovnersion)
{
	colorType = RGB; // the color is RGB when we use matlab to load the image
	if(mxIsClass(matrix,"uint8"))
	{
		LoadMatlabImageCore<unsigned char>(matrix,IsImageScaleCovnersion);
		return true;
	}
	if(mxIsClass(matrix,"int8"))
	{
		LoadMatlabImageCore<char>(matrix,IsImageScaleCovnersion);
		return true;
	}
	if(mxIsClass(matrix,"int32"))
	{
		LoadMatlabImageCore<int>(matrix,IsImageScaleCovnersion);
		return true;
	}
	if(mxIsClass(matrix,"uint32"))
	{
		LoadMatlabImageCore<unsigned int>(matrix,IsImageScaleCovnersion);
		return true;
	}
	if(mxIsClass(matrix,"int16"))
	{
		LoadMatlabImageCore<short int>(matrix,IsImageScaleCovnersion);
		return true;
	}
	if(mxIsClass(matrix,"uint16"))
	{
		LoadMatlabImageCore<unsigned short int>(matrix,IsImageScaleCovnersion);
		return true;
	}
	if(mxIsClass(matrix,"single"))
	{
		LoadMatlabImageCore<float>(matrix,IsImageScaleCovnersion);
		return true;
	}
	if(mxIsClass(matrix,"float"))
	{
		LoadMatlabImageCore<float>(matrix,IsImageScaleCovnersion);
		return true;
	}
	mexErrMsgTxt("Unknown type of the image!");
	return false;
}


template <class T>
template <class T1>
void Image<T>::ConvertFromMatlab(const T1 *pMatlabPlane, int _width, int _height, int _nchannels)
{
	if(imWidth!=_width || imHeight!=_height || nChannels!=_nchannels)
		allocate(_width,_height,_nchannels);
	int offset=0;
	for(int i=0;i<imHeight;i++)
		for(int j=0;j<imWidth;j++)
			for(int k=0;k<nChannels;k++)
				pData[offset++]=pMatlabPlane[k*nPixels+j*imHeight+i];
}

// convert image data to matlab matrix
template <class T>
template <class T1>
void Image<T>::ConvertToMatlab(T1 *pMatlabPlane) const
{
	int offset=0;
	for(int i=0;i<imHeight;i++)
		for(int j=0;j<imWidth;j++)
			for(int k=0;k<nChannels;k++)
				pMatlabPlane[k*nPixels+j*imHeight+i]=pData[offset++];
}

template <class T>
void Image<T>::OutputToMatlab(mxArray *&matrix) const
{
	int dims[3];
	dims[0]=imHeight;
	dims[1]=imWidth;
	dims[2]=nChannels;
	int nDims;
	nDims = (nChannels ==1)? 2:3;
	if(typeid(T) == typeid(unsigned char))
		matrix=mxCreateNumericArray(nDims, dims,mxUINT8_CLASS, mxREAL);
	if(typeid(T) == typeid(char))
		matrix=mxCreateNumericArray(nDims, dims,mxINT8_CLASS, mxREAL);
	if(typeid(T) == typeid(short int))
		matrix=mxCreateNumericArray(nDims, dims,mxINT16_CLASS, mxREAL);
	if(typeid(T) == typeid(unsigned short int))
		matrix=mxCreateNumericArray(nDims, dims,mxUINT16_CLASS, mxREAL);
	if(typeid(T) == typeid(int))
		matrix=mxCreateNumericArray(nDims, dims,mxINT32_CLASS, mxREAL);
	if(typeid(T) == typeid(unsigned int))
		matrix=mxCreateNumericArray(nDims, dims,mxUINT32_CLASS, mxREAL);
	if(typeid(T) == typeid(float))
		matrix=mxCreateNumericArray(nDims, dims,mxSINGLE_CLASS, mxREAL);
	if(typeid(T) == typeid(float))
		matrix=mxCreateNumericArray(nDims, dims,mxfloat_CLASS, mxREAL);

	ConvertToMatlab<T>((T*)mxGetData(matrix));
}

#endif

}  // namespace cpm

#endif // _IMAGE_HPP


