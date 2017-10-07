#ifndef __SIMD_H_
#define __SIMD_H_


#ifdef __SSE2__

    #define SIMD_MESSAGE "Using SSE2 instruction set"
    #pragma message( SIMD_MESSAGE )
    
    
    #define WITH_SSE

    #include <xmmintrin.h>
#else
    
    #define SIMD_MESSAGE "Compiling without SSE2 instruction set"
    
#endif // __SSE2__

#undef SIMD_MESSAGE

#endif // __SIMD_H_
