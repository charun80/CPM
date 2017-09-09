#ifndef _MEMORY_h
#define _MEMORY_h

#include <memory.h>
#include <malloc.h>


#ifdef WITH_SSE
#include <xmmintrin.h>
#endif

template <class T>
inline void* xmalloc(T size){
#ifdef WITH_SSE
#ifdef WIN32
	return _aligned_malloc(size, 32);
#else
	return memalign(32, size);
#endif
#else
	return malloc(size);
#endif
}

template <class T>
inline void xfree(T* ptr){
#if defined(WITH_SSE) && defined(WIN32)
	_aligned_free(ptr);
#else
	free(ptr);
#endif
}

#ifdef WITH_SSE

// for windows and linux
typedef union _m128{
	__m128 m;
	__m128i mi;
	float m128_f32[4];
	unsigned short m128i_u16[8];
}hu_m128;

#endif


#endif
