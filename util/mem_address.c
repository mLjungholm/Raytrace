#include "mex.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *ptr;
    int i;
    mwSize n;
      /* Function accepts only one input argument */
      if(nrhs > 1)
      {
          mexErrMsgTxt("Only one input argument accepted\n\n");
      }   
      /* This function works only with row vectors */
      if( mxGetM(prhs[0]) > 1 )
      {
          mexErrMsgTxt("Only row vectors are accepted\n\n");
      }
      /* Get the Pointer to the Real Part of the mxArray */
      ptr = mxGetPr(prhs[0]);
      n = mxGetN(prhs[0]);
      mexPrintf("Pointer Address of mxArray = %x\n\n", prhs[0]);
      for(i = 0; i < n; i++)
      {
          mexPrintf("Value of Real Argument at array index %d = %f\n", i, ptr[i]);
          mexPrintf("Pointer Address = %x\n\n", &ptr[i]);
      }
  }