#ifndef __SIMD_AVX_H_
#define __SIMD_AVX_H_

#include <immintrin.h>
#include <stdint.h>


// SIMD vector type
typedef __v8sf   simdsf_t;
typedef __v16hi  simdhi_t;
typedef __v32qi  simdqi_t;


inline static simdsf_t simdsf_init( float x );

#define simdsf_sqrt(x)   __builtin_ia32_sqrtps256(x)
#define simdsf_max(x,y)  __builtin_ia32_maxps256(x,y)

inline static uint32_t simdqi_sumAbsDiff( const simdqi_t &x, const simdqi_t &y );



// Implementation


inline static simdsf_t simdsf_init( float x ) 
{ 
    simdsf_t sx = {x,x,x,x, x,x,x,x};
    return ( sx ); 
}

#ifndef __AVX2__
    // No AVX2 available
    
    inline static uint32_t simdqi_sumAbsDiff( const simdqi_t &x, const simdqi_t &y )
    {
        const __v16qi *X2 = (__v16qi*)&x;
        const __v16qi *Y2 = (__v16qi*)&y;
        
        __v2di z;
        {
            __v2di Z[2];
            Z[0] = __builtin_ia32_psadbw128 ( X2[0], Y2[0] );
            Z[1] = __builtin_ia32_psadbw128 ( X2[1], Y2[1] );
            
            z = Z[0] + Z[1];
        }
        
        return ( z[0] + z[1] );
    }
    
    
#else
    // AVX2 available
    
    inline static uint32_t simdqi_sumAbsDiff( const simdqi_t &x, const simdqi_t &y )
    {
        // see: https://software.intel.com/sites/landingpage/IntrinsicsGuide/#text=sad&techs=AVX2&expand=4559
         
        union {
            __v16hi vhi;
            __v8hi  v2hi[2];
        } z4;
        
        z4.vhi = __builtin_ia32_psadbw256( x, y ); // only 4 values are used
        
        union {
            __v8hi vhi;
            __v4si vdi;
        } z2;
        // using SSE here -
        z2.vhi = z4.v2hi[0] + z4.v2hi[1];
        
        // ... the low 32 bits of 64-bit elements
        return (z2.vdi[0] + z2.vdi[2]);
    }

    
    
#endif // __AVX2__


#endif // __SIMD_AVX_H_

