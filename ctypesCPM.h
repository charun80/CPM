#ifndef PYCPM_H_INCLUDED
#define PYCPM_H_INCLUDED


// example: see https://github.com/nlgranger/CPPyModule/blob/master/example-distutils/flipping_wrapper.h


#define DLL_PUBLIC __attribute__ ((visibility ("default")))

extern "C"
{

struct sFlowResult
{
    float* m_outMatching_pf;
    int    m_nMatches_i;
};



DLL_PUBLIC sFlowResult cpmFlowFromBGR(
    const float* const f_inImg1_data_pf,
    const float* const f_inImg2_data_pf,
    const int f_nRows_i, const int f_nCols_i, const int f_nChannels,
    const int f_nSteps );


DLL_PUBLIC void clearFlowResult( sFlowResult* const f_flowRes );



} // extern C

#endif // PYCPM_H_INCLUDED
