#ifndef STOCHASTIC_H
#define STOCHASTIC_H




namespace cpm
{


enum SortType{SortAscending,SortDescending};

class CStochastic
{
public:
	CStochastic(void);
	~CStochastic(void);
	static void ConvertInt2String(int x,char* string,int BitNumber=3);
	static double UniformSampling();
	static int UniformSampling(int R);
	static double GaussianSampling();
	template <class T> static void GetMeanVar(T* signal,int length,double* mean,double* var);
	static int Sampling(double* Density,int NumSamples);
	static double GetMean(double *signal,int length);
	static void Generate1DGaussian(double* pGaussian,int size,double sigma=0);
	static void Generate2DGaussian(double* pGaussian,int size,double sigma=0);
	static double entropy(double* pDensity,int n);

	template <class T> static T sum(int NumData,T* pData);
	template <class T> static void Normalize(int NumData,T* pData);
	template <class T> static T mean(int NumData, T* pData);
	template <class T> static void sort(int number, T* pData,int *pIndex,SortType m_SortType=SortDescending);
	template <class T> static T Min(int NumData, T* pData);
	template <class T> static T Min(int NumData, T* pData1,T* pData2);
	template <class T> static T Max(int NumData ,T* pData);
	template <class T> static int FindMax(int NumData,T* pData);
	template <class T1,class T2> static void ComputeVectorMean(int Dim,int NumData,T1* pData,T2* pMean,double* pWeight=NULL);
	template <class T1,class T2> static void ComputeMeanCovariance(int Dim,int NumData,T1* pData,T2* pMean,T2* pCovarance,double* pWeight=NULL);
	template <class T1,class T2> static double VectorSquareDistance(int Dim,T1* pVector1,T2* pVector2);
	template <class T1> static void KMeanClustering(int Dim,int NumData,int NumClusters,T1* pData,int *pPartition,double** pClusterMean=NULL,int MaxIterationNum=10,int MinClusterSampleNumber=2);
	template <class T> static double norm(T* X,int Dim);
	template <class T1,class T2> static int FindClosestPoint(T1* pPointSet,int NumPoints,int nDim,T2* QueryPoint);
	template <class T1,class T2> static void GaussianFiltering(T1* pSrcArray,T2* pDstArray,int NumPoints,int nChannels,int size,double sigma);
};


}  // namespace cpm


#include "Stochastic.hpp"

#endif  // STOCHASTIC_H
