#ifndef __SIMD_NEON_H_
#define __SIMD_NEON_H_


#include <arm_neon.h>

// see https://gcc.gnu.org/onlinedocs/gcc-4.9.1/gcc/ARM-NEON-Intrinsics.html
//     https://github.com/thenifty/neon-guide
//     http://infocenter.arm.com/help/topic/com.arm.doc.dui0491i/Badcdfad.html



typedef float32x4  simdsf_t;
typedef uint8x16_t simdqi_t;


inline static simdsf_t simdsf_init( float x );
static inline simdsf_t simdsf_sqrt( simdsf_t x  );
static inline simdsf_t simdsf_max( simdsf_t x, simdsf_t y );

inline static unsigned int simdqi_sumAbsDiff( const simdqi_t &x, const simdqi_t &y );



// Implementation


inline static unsigned int simdqi_sumAbsDiff( const simdqi_t &x, const simdqi_t &y )
{
    uint8x16_t z16 = vabdq_u8( x, y );
    
    uint8x8_t z16high = vget_high_u8( z16 );
    uint8x8_t z16low  = vget_low_u8( z16 );
    
    // 8 values left
    uint16x8_t z8 = vaddl_u8( z16high, z16low );
    
    uint16x4_t z8high = vget_high_u16( z8 );
    uint16x4_t z8low  = vget_low_u16( z8 );
    
    // 4 values left
    uint32x4_t z4 = vmovl_u16( vadd_u16( z8high, z8low ) );
    
    uint32x2_t z4high = vget_high_u32( z4 );
    uint32x2_t z4low  = vget_low_u32( z4 );
    
    // 2 values left
    uint32x2_t z2 = vadd_u32( z4high, z4low );
    
    return (z2[0] + z2[1]);
}



static inline simdsf_t simdsf_sqrt( simdsf_t x  )
{
    return (1.f / vrsqrteq_f32(x));
}


static inline simdsf_t simdsf_max( const simdsf_t &x, const simdsf_t &y )
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

