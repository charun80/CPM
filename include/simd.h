#ifndef __SIMD_H_
#define __SIMD_H_

#include <string>
#include <iostream>
#include <cassert>

// see https://software.intel.com/sites/landingpage/IntrinsicsGuide
//     https://gcc.gnu.org/onlinedocs/gcc/x86-Built-in-Functions.html


#ifdef __AVX__
    
    #define SIMD_MESSAGE "Using AVX instruction set"
    #define NSIMDBYTES 32   // Number of bytes in SIMD registers
    #define USE_SIMD

    
    #include "simd_avx.h"
    
#else // __AVX__
#ifdef __SSE__

    #define SIMD_MESSAGE "Using SSE instruction set"
    #define NSIMDBYTES 16   // Number of bytes in SIMD registers
    #define USE_SIMD

    #include "simd_sse.h"

#else // __SSE__
#if defined __ARM_NEON__ || defined __ARM_NEON

    #define SIMD_MESSAGE "Using NEON instruction set"
    #define NSIMDBYTES 16  // Number of bytes in SIMD registers
    #define USE_SIMD

    #include "simd_neon.h"

#else // __NEON__

    #define SIMD_MESSAGE "Neither SSE nor AVX instruction set available - no compilation possible"
    #error SIMD_MESSAGE
    

#endif // __NEON__
#endif // __SSE__
#endif // __AVX__



#ifdef USE_SIMD
    // common definitions in case simd is available

    static const size_t NSimdBytes = NSIMDBYTES;
    static const size_t NSimdFloats = NSIMDBYTES / sizeof(float);
    static const size_t NSimdChars  = NSIMDBYTES / sizeof(char);
    
    #undef NSIMDBYTES
    
    
    inline static bool isAligned( const void *fptr )
    {
        size_t ptr_i = size_t(fptr);
        return (0 == (ptr_i % NSimdBytes));
    }
    
    
    inline static simdqi_t* simdqi_ptrcast( unsigned char* fptr ) {  return reinterpret_cast<simdqi_t*>( fptr ); }
    inline static const simdqi_t* simdqi_ptrcast( const unsigned char* fptr ) {  return reinterpret_cast<const simdqi_t*>( fptr ); }
    
    inline static simdsf_t* simdsf_ptrcast( float* fptr ) { return reinterpret_cast<simdsf_t*>( fptr ); }
    inline static const simdsf_t* simdsf_ptrcast( const float* fptr ) { return reinterpret_cast<const simdsf_t*>( fptr ); }


    inline static void* __xmalloc( size_t fSize, const std::string &fErrpos )
    {
        fSize = ((fSize + NSimdBytes - 1) / NSimdBytes) * NSimdBytes;
        assert( 0 == (fSize % NSimdBytes) );
        
        void *lMemPtr = NULL;
        const int lMemAlignError = posix_memalign( &lMemPtr, NSimdBytes, fSize );
        
        if( 0 != lMemAlignError )
        {
            std::cerr << "Error: allocating memory in " << fErrpos << ": error = " << lMemAlignError << std::endl;
            exit(1);
        }
        
        return lMemPtr;
    }

#else  // #ifdef USE_SIMD


    inline static void* __xmalloc( size_t fSize, const std::string &fErrpos )
    {
        void *lMemPtr = malloc( fSize );
        const int lMemAlignError = 
        
        if( NULL != lMemPtr )
        {
            std::cerr << "Error: allocating memory in " << fErrpos << std::endl;
            exit(1);
        }
        
        return lMemPtr;
    }


#endif // USE_SIMD


#define __SIMD_STR_HELPER(x) #x
#define __SIMD_STR(x) __SIMD_STR_HELPER(x)

//#define xmalloc(Size) __xmalloc( Size, " " )
#define xmalloc(Size) __xmalloc( Size,  __FILE__ ":" __SIMD_STR( __LINE__ ) )


#pragma message( SIMD_MESSAGE )
#undef SIMD_MESSAGE


#endif // __SIMD_H_

