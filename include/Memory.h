#ifndef _MEMORY_h
#define _MEMORY_h

#include <memory.h>
#include <malloc.h>


#ifdef USE_SIMD
#include <xmmintrin.h>
#endif

namespace cpm
{


#ifdef USE_SIMD

// for windows and linux
typedef union _m128{
	__m128 m;
	__m128i mi;
	float m128_f32[4];
	unsigned short m128i_u16[8];
}hu_m128;

#endif

}  // namespace cpm

#endif
