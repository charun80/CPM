#include "ctypesCPM.h"

#include "CPM.h"
#include <cassert>


extern "C"
int cpmFlowFromGray(
    const float* const f_inImg1_data_pf,
    const float* const f_inImg2_data_pf,
    const int f_nRows_i, const int f_nCols_i,
    const int f_nSteps,
    float **f_outMatching_ppf )
{
    const size_t l_numPixels = f_nRows_i * f_nCols_i;
    const size_t l_imgSize = sizeof(float) * l_numPixels;

    // obtain input
    FImage l_cpmImg1, l_cpmImg2;

    l_cpmImg1.allocate( f_nCols_i, f_nRows_i, 1);
    memcpy( l_cpmImg1.pData, f_inImg1_data_pf, l_imgSize );

    l_cpmImg2.allocate( f_nCols_i, f_nRows_i, 1);
    memcpy( l_cpmImg2.pData, f_inImg2_data_pf, l_imgSize );

    // compute flow
    FImage l_cpmMatches;
    {
        CPM cpm;
        cpm.SetStep( f_nSteps );
        cpm.Matching( l_cpmImg1, l_cpmImg1, l_cpmMatches );
    }

    const int l_numMatches = l_cpmMatches.height();

    if ( 0 < l_numMatches ) {
        assert( 4 == l_cpmMatches.width() );
        assert( 1 == l_cpmMatches.nchannels() );

        // Generate output
        *f_outMatching_ppf = l_cpmMatches.pData;
        l_cpmMatches.pData = NULL;
    }
    else {
        *f_outMatching_ppf = NULL;
    }

    return l_numMatches;
}

