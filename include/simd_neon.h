#ifndef __SIMD_NEON_H_
#define __SIMD_NEON_H_


#include <arm_neon.h>


typedef float32x4 simdsf_t;


static inline simdsf_t simdsf_sqrt( simdsf_t x  )
{
    return (1.f / vrsqrteq_f32(x));
}


static inline simdsf_t simdsf_max( simdsf_t x, simdsf_t y )
{
    simdsf_t z;
    
    #define SIMDSF_MAX(a,b)  (((a)>(b)) ? (a) : (b))
    
    for (int i = 0; i < 4; ++i)
        z[i] = MAX( x[i], y[i] );
    
    #undef SIMDSF_MAX
    
    return z;
}


inline static simdsf_t simdsf_init( float x ) 
{ 
    simdsf_t sx = {x,x,x,x};
    return ( sx ); 
}


#endif // __SIMD_NEON_H_

