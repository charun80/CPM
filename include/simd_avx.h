#ifndef __SIMD_AVX_H_
#define __SIMD_AVX_H_

#include <immintrin.h>


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
    
    
#endif // __AVX2__


#endif // __SIMD_AVX_H_

