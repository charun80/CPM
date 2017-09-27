#ifndef __SIMD_H_
#define __SIMD_H_


#define MAX(a,b)  (((a)>(b)) ? (a) : (b))


// see https://software.intel.com/sites/landingpage/IntrinsicsGuide
//     https://gcc.gnu.org/onlinedocs/gcc/x86-Built-in-Functions.html


#ifdef __AVX__
    
    #define SIMD_MESSAGE "Using AVX instruction set"
    #define USE_SIMD


    #include <immintrin.h>


    // Number of bytes in SIMD registers
    #define NSIMDBYTES 32

    // SIMD vector type
    typedef __v8sf   simdsf_t;
    typedef __v16hi  simdhi_t;
    typedef __v32qi  simdqi_t;
    
    #define simdsf_sqrt(x)   __builtin_ia32_sqrtps256(x)
    #define simdsf_max(x,y)  __builtin_ia32_maxps256(x,y)
    
    inline static simdsf_t simdsf_init( float x ) 
    { 
        simdsf_t sx = {x,x,x,x, x,x,x,x};
        return ( sx ); 
    }
    
    #ifndef __AVX2__
        // No AVX2 available
        
        inline static unsigned int simdqi_sumAbsDiff( const simdqi_t &x, const simdqi_t &y )
        {           
            typedef unsigned int uint;
            
            const __v16qi *X2 = (__v16qi*)&x;
            const __v16qi *Y2 = (__v16qi*)&y;
            
            __v2di z;
            {
                __v2di Z[2];
                Z[0] = __builtin_ia32_psadbw128 ( X2[0], Y2[0] );
                Z[1] = __builtin_ia32_psadbw128 ( X2[1], Y2[1] );
                
                __v2di z = Z[0] + Z[1];
            }
            
            return ( uint(z[0]) + uint(z[1]) );
        }
        
    #else
        // AVX2 available
        
        inline static unsigned int simdqi_sumAbsDiff( simdqi_t x, simdqi_t y )
        {
            typedef unsigned int uint;
            __v8hi z;
            {
                typedef union {
                    __v16hi vd;
                    __v8hi  vs[2];
                } v16hi_2v8hi;
                
                v16hi_2v8hi Z;
                Z.vd = __builtin_ia32_psadbw256( x, y );
                
                // using SSE here
                z = Z.vs[0] + Z.vs[1];
            }
            
            uint zui = z[0];
            for (int i = 1; i < 8; ++i)
                zui += uint(z[i]);
            
            return zui;
        }
        
        
    #endif // __SSE2__

#else
#ifdef __SSE__

    #define SIMD_MESSAGE "Using SSE instruction set"
    #define USE_SIMD

    #include <xmmintrin.h>

    // Number of bytes in SIMD registers
    #define NSIMDBYTES 16
    
    typedef __v4sf  simdsf_t;
    typedef __v16qi simdqi_t;
    
    // Avaliable in SSE1
    #define simdsf_sqrt(x)   __builtin_ia32_sqrtps(x)
    #define simdsf_max(x,y)  __builtin_ia32_maxps(x,y)
    
    
    #ifdef __SSE2__
    
    // Available in SSE2
    #define simd8i_absdiff_2i(x,y) __builtin_ia32_psadbw128(x,y)
    
    
    
    
    inline static simdsf_t simdsf_init( float x ) 
    { 
        simdsf_t sx = {x,x,x,x};
        return ( sx ); 
    }
    
    
    #ifndef __SSE2__
        // No SSE2 available
        
        inline static unsigned int simdqi_sumAbsDiff( const simdqi_t &x, const simdqi_t &y )
        {
            const __v8qi *X2 = (__v8qi*)&x;
            const __v8qi *Y2 = (__v8qi*)&y;
            
            __v1di Z[2];
            
            Z[0] = __builtin_ia32_psadbw ( X2[0], Y2[0] );
            Z[1] = __builtin_ia32_psadbw ( X2[1], Y2[1] );
            return ( (unsigned int)(Z[0]) + (unsigned int)(Z[1]) );
        }
        
    #else
        // SSE2 available
        
        inline static unsigned int simdqi_sumAbsDiff( simdqi_t x, simdqi_t y )
        {
            __v2di z = __builtin_ia32_psadbw128( x, y );
            return (z[0] + z[1]);
        }
        
        
    #endif // __SSE2__

#else
#if defined __ARM_NEON__ || defined __ARM_NEON__

    #define SIMD_MESSAGE "Using NEON instruction set"
    #define USE_SIMD

    #include <arm_neon.h>
    
    // Number of bytes in SIMD registers
    #define NSIMDBYTES 16
    
    typedef float32x4 simdsf_t;
    
    static inline simdsf_t simdsf_sqrt( simdsf_t x  )
    {
        return (1.f / vrsqrteq_f32(x));
    }
    
    
    static inline simdsf_t simdsf_max( simdsf_t x, simdsf_t y )
    {
        simdsf_t z;
        
        for (int i = 0; i < 4; ++i)
            z[i] = MAX( x[i], y[i] );
        
        return z;
    }
    
    
    inline static simdsf_t simdsf_init( float x ) 
    { 
        simdsf_t sx = {x,x,x,x};
        return ( sx ); 
    }
    

#else

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


    inline static simdsf_t* simdsf_ptrcast( void* fptr ) 
    {
        return ((simdsf_t*) fptr);
    }


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


#define __STR_HELPER(x) #x
#define __STR(x) __STR_HELPER(x)

#define xmallox(Size) __xmalloc( Size, __FILE__ ":" __STR(__LINE__) )

#undef __STR_HELPER
#undef __STR


#pragma message( SIMD_MESSAGE )
#undef SIMD_MESSAGE

#undef NSIMDBYTES

#endif // __SIMD_H_
