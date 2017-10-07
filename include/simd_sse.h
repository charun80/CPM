#ifndef __SIMD_SSE_H_
#define __SIMD_SSE_H_


#include <xmmintrin.h>


typedef __v4sf  simdsf_t;
typedef __v16qi simdqi_t;

// Avaliable in SSE1
#define simdsf_sqrt(x)   __builtin_ia32_sqrtps(x)
#define simdsf_max(x,y)  __builtin_ia32_maxps(x,y)

    

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


#endif // __SIMD_SSE_H_

