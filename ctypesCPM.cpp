#include "ctypesCPM.h"

#include "CPM.cpp"
#include <cassert>


extern "C"
DLL_PUBLIC sFlowResult cpmFlowFromBGR(
    const float* const f_inImg1_data_pf,
    const float* const f_inImg2_data_pf,
    const int f_nRows_i, const int f_nCols_i,
    const int f_nSteps )
{
    sFlowResult l_Result = {NULL,0};

    const size_t l_numPixels = f_nRows_i * f_nCols_i * 3;
    const size_t l_imgSize = sizeof(float) * l_numPixels;

    // obtain input
    FImage l_cpmImg1, l_cpmImg2;

    l_cpmImg1.allocate( f_nCols_i, f_nRows_i, 3);
    //l_cpmImg1.allocate( f_nRows_i, f_nCols_i, 3);
    memcpy( l_cpmImg1.pData, f_inImg1_data_pf, l_imgSize );
    l_cpmImg1.setColorType ( BGR );

    l_cpmImg2.allocate( f_nCols_i, f_nRows_i, 3);
    //l_cpmImg2.allocate( f_nRows_i, f_nCols_i, 3);
    memcpy( l_cpmImg2.pData, f_inImg2_data_pf, l_imgSize );
    l_cpmImg2.setColorType ( BGR );

    // compute flow
    FImage l_cpmMatches;
    {
        CPM cpm;
        cpm.SetStep( f_nSteps );
        std::cout << "Start Matching" << std::endl;
        cpm.Matching( l_cpmImg1, l_cpmImg2, l_cpmMatches );
        std::cout << "Done Matching" << std::endl;
    }

    const int l_numMatches = l_cpmMatches.height();

    if ( 0 < l_numMatches ) {
        assert( 4 == l_cpmMatches.width() );
        assert( 1 == l_cpmMatches.nchannels() );

        // Generate output
        l_Result.m_outMatching_pf = l_cpmMatches.pData;
        l_cpmMatches.pData = NULL;

        l_Result.m_nMatches_i = l_numMatches;
    }

    return l_Result;
}

