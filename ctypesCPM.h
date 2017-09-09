#ifndef PYCPM_H_INCLUDED
#define PYCPM_H_INCLUDED


// example: see https://github.com/nlgranger/CPPyModule/blob/master/example-distutils/flipping_wrapper.h


extern "C"
int cpmFlowFromGray(
    const float* const f_inImg1_data_pf,
    const float* const f_inImg2_data_pf,
    const int f_nRows_i, const int f_nCols_i,
    const int f_nSteps,
    float **f_outMatching_ppf );

#endif // PYCPM_H_INCLUDED
