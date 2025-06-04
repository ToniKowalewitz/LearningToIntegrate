/*----------------------------------------------------------------------*/
/*
  MATLAB (R) is a trademark of The Mathworks (R) Corporation

  Function:    complex2real
  Filename:    complex2real.c
  Programmer:  James Tursa
  Version:     3.00
  Date:        March 13, 2020
  Copyright:   (c) 2018, 2019, 2020 by James Tursa, All Rights Reserved

  Revision History: 1.00 , 2018-10-19
                    * Initial Release
                    2.00 , 2020-02-27
                    * Updated for R2019a CrossLink field change
                    3.00 , 2020-03-13
                    * Updated for R2019b mxArray header changes (ndim)

  This code uses the BSD License:

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

     * Redistributions of source code must retain the above copyright
       notice, this list of conditions and the following disclaimer.
     * Redistributions in binary form must reproduce the above copyright
       notice, this list of conditions and the following disclaimer in
       the documentation and/or other materials provided with the distribution

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
  POSSIBILITY OF SUCH DAMAGE.

 -------------------------------------------------------------------------

 Building:

 COMPLEX2REAL requires that a mex routine be built (one time only). This
 process is typically self-building the first time you call the function
 as long as you have the files complex2real.m and complex2real.c in the
 same directory somewhere on the MATLAB path. If you need to manually
 build the mex function, here are the commands:

 >> cd ______ <-- change directory to where complex2real.c is

 and then

 >> mex complex2real.c -R2018a

 Note that you cannot use COMPLEX2REAL on versions earlier than R2018a.

 -------------------------------------------------------------------------

 *** DISCLAIMER ***

 This mex routine uses undocumented API functions and hacks into the mxArray
 itself in order to create the converted shared data copy.  I have tried to
 code it in such a way that the routine will not crash MATLAB, but caveat
 emptor.

 -------------------------------------------------------------------------

 Usage:

 Y = complex2real(X)

 Where X = any full numeric variable
       Y = X reinterpreted inplace as a real variable (shared data copy)

 If X is complex, then the first dimension of Y is 2x the first dimension of X.

*/
/* -------------------------------------------------------------------- */
/* -------------------------------------------------------------------- */
/* -------------------------------------------------------------------- */

/* Includes ----------------------------------------------------------- */

#include "mex.h"

#include "matlab_version.h"
#include "matlab_version.c"

#include <stdio.h>  /* sprintf , sscanf */

/* Undocumented API function ------------------------------------------ */

#undef mxCreateSharedDataCopy /* Will have to replace this hack at some point */
mxArray *mxCreateSharedDataCopy(const mxArray *);

/* Top level version needed for mxArray HACK stuff */

int version = 0;

/* mxArray header hack ------------------------------------------------ */
/* Need this to set the pr & pi data pointers directly, and also to use */
/* for CrossLink shared data copy checking.                             */

#include "mxArray_header.h"
/* flags bits 17-12 values for the number of data bytes
 *                      0                000000
 *                      1 -    32        000001
 *                     33 -    64        000011
 *                     65 -    96        000101
 *                     97 -   128        000111
 *                    129 -   160        001001
 *                    161 -   192        001011
 *                    193 -   224        001101
 *                    225 -   256        001111
 *                    257 -   352        010001
 *                    353 -   960        010011
 *                    961 -  2048        010101
 *                   2049 -  4096        010111
 *                   4097 -  6144        011001
 *                   6145 -  8192        011011
 *                   8193 -  9216        011101
 *                   9217 - 10240        011111
 *                  10241 - 13312        100001
 *                  13313 - 16384        100011
 *                  16385 -              000000
 */

#define NBYTES 18
size_t flags_nbytes[NBYTES] = {32,64,96,128,160,192,224,256,352,960,2048,
                               4096,6144,8192,9216,10240,13312,16384};

#define NBYTES_SHIFT 12
#define NBYTES_BITS ( 31u << NBYTES_SHIFT )

#define ISCOMPLEX_SHIFT 11
#define ISCOMPLEX_BIT ( 1u << ISCOMPLEX_SHIFT )

#define ISFULLREALDOUBLESCALAR_BIT 1u

/* -------------------------------------------------------------------- */
/* Calculate the flags bit pattern for the number of bytes in object    */
/* -------------------------------------------------------------------- */

unsigned int flags_nbytes_bits(size_t n)
{
    unsigned int f;

    if( n ) {
        for( f = 0; f < NBYTES; f++ ) {
            if( n <= flags_nbytes[f] ) {
                return ((f<<1)|1u) << NBYTES_SHIFT;
            }
        }
    }
    return 0;
}

/* -------------------------------------------------------------------- */
/* Fixes the number of bytes bits field                                 */
/* -------------------------------------------------------------------- */

void fix_flags_nbytes(mxArray *mx)
{
    MXH
    unsigned int flags;

    if( version >= 0x2019b ) {
        HACK(mx);
        flags = HACKGET(flags);
        flags &= ~NBYTES_BITS; /* reset the number of data bytes bits to start with */
        if( (mxIsNumeric(mx) || mxIsLogical(mx)) && !(flags & ISCOMPLEX_BIT) ) { /* if variable is real or logical */
            flags |= flags_nbytes_bits(mxGetNumberOfElements(mx) * mxGetElementSize(mx));
        }
        HACKSET(flags, flags );
    }
}

/* -------------------------------------------------------------------- */
/* Turn a real mxArray into a complex mxArray inplace                   */
/* Assumes the first dimension has already been checked for being even  */
/* -------------------------------------------------------------------- */

void real2complex(mxArray *mx)
{
    MXH

    HACK(mx);
    HACKSET(flags, HACKGET(flags) | ISCOMPLEX_BIT ); /* Set the iscomplex bit */
    mxSetM(mx,mxGetM(mx)/2); /* Halve the first dimension */
    fix_flags_nbytes(mx); /* Fix up the number of bytes bits */
}

/* -------------------------------------------------------------------- */
/* Turn a complex mxArray into a real mxArray inplace                   */
/* -------------------------------------------------------------------- */

void complex2real(mxArray *mx)
{
    MXH

    HACK(mx);
    HACKSET(flags, HACKGET(flags) & ~ISCOMPLEX_BIT ); /* Reset the iscomplex bit */
    mxSetM(mx,mxGetM(mx)*2); /* Double the first dimension */
    fix_flags_nbytes(mx); /* Fix up the number of bytes bits */
}

/* -------------------------------------------------------------------- */
/* Gateway                                                              */
/* -------------------------------------------------------------------- */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

/* Version check */

    if( TARGET_API_VERSION == R2017b ) {
        mexErrMsgTxt("This code cannot be used with R2017b memory model");
    }
    if( !version ) {
        version = matlab_version();
    }

/* Input/Output argument checks */

    if( nlhs > 1 ) {
        mexErrMsgTxt("Too many outputs.");
    }
    if( nrhs != 1 || !mxIsNumeric(prhs[0]) || mxIsSparse(prhs[0]) ) {
        mexErrMsgTxt("Need exactly one full numeric input.");
    }

/* Create the shared data copy, then change it inplace */

    plhs[0] = mxCreateSharedDataCopy(prhs[0]); /* Create shared data copy */
    if( mxIsComplex(plhs[0]) ) {
        complex2real(plhs[0]);
    }
}
