#ifndef __SIMD_SSE_H_
#define __SIMD_SSE_H_


#include <xmmintrin.h>
#include <stdint.h>


typedef __v4sf  simdsf_t;
typedef __v16qi simdqi_t;


inline static simdsf_t simdsf_init( float x );

// Avaliable in SSE1
#define simdsf_sqrt(x)   __builtin_ia32_sqrtps(x)
#define simdsf_max(x,y)  __builtin_ia32_maxps(x,y)

inline static uint32_t simdqi_sumAbsDiff( const simdqi_t &x, const simdqi_t &y );


// Implementation


inline static simdsf_t simdsf_init( float x ) 
{ 
    simdsf_t sx = {x,x,x,x};
    return ( sx ); 
}



#ifndef __SSE2__
    // No SSE2 available
    
    inline static uint32_t simdqi_sumAbsDiff( const simdqi_t &x, const simdqi_t &y )
    {
        const __v8qi *X2 = reinterpret_cast<const __v8qi*>(&x);
        const __v8qi *Y2 = reinterpret_cast<const __v8qi*>(&y);
        
        __v1di z = __builtin_ia32_psadbw ( X2[0], Y2[0] ) + __builtin_ia32_psadbw ( X2[1], Y2[1] );
        
        return z[0];
    }
    
    
#else
    // SSE2 available
    
    inline static uint32_t simdqi_sumAbsDiff( const simdqi_t &x, const simdqi_t &y )
    {
        __v2di z = __builtin_ia32_psadbw128( x, y );
        return (z[0] + z[1]);
    }
    
    
    
#endif // __SSE2__


#endif // __SIMD_SSE_H_

