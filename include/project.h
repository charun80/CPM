#ifndef _PROJECT_H
#define _PROJECT_H

#include <cstdio>
#include <vector>

namespace cpm
{


template <class T>
void _Release1DBuffer(T* pBuffer)
{
	if(pBuffer!=NULL)
		delete []pBuffer;
	pBuffer=NULL;
}

template <class T>
void _Rlease2DBuffer(T** pBuffer,size_t nElements)
{
	for(size_t i=0;i<nElements;i++)
		delete [](pBuffer[i]);
	delete []pBuffer;
	pBuffer=NULL;
}

}  // namespace cpm

#ifdef _MATLAB
#include "mex.h"
#endif


#ifndef WIN32
    #define strcmpi strcasecmp
#endif

#endif  // _PROJECT_H
