#include <string.h>

#include <mex.h>

const size_t ptr_size = sizeof(void*);
const mxClassID ptr_class = sizeof(void*) == 8 ? mxUINT64_CLASS : mxUINT32_CLASS;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mxArray *ret_ptr;
  void *ret_ptr_data;
  mxDouble *rhs_pr;

  /* Validate argument count. */
  if (nlhs > 1 || nrhs < 1 || nrhs > 1) {
    mexErrMsgTxt("Incorrect number of input and/or output arguments.");
  }

  /* Create scalar of appropriate class for return value. */
  ret_ptr = mxCreateNumericMatrix(1, 1, ptr_class, mxREAL);
  ret_ptr_data = mxGetData(ret_ptr);

  /* Get pointer to real part of array of doubles. */
  rhs_pr = mxGetPr(prhs[0]);

  /* Copy the pointer value to return value. */
  memcpy(ret_ptr_data, &rhs_pr, ptr_size);

  /* Assign return value. */
  plhs[0] = ret_ptr;
}
