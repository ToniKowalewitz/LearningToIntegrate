/* mxArray_header.h */
/*************************************************************************************
 *
 * MATLAB (R) is a trademark of The Mathworks (R) Corporation
 *
 * Filename:    mxArray_header.h
 * Programmer:  James Tursa
 * Version:     1.00
 * Date:        March 13, 2020
 * Copyright:   (c) 2020 by James Tursa, All Rights Reserved
 *
 * Change Log:
 * 2020/Mar/13 --> 1.00, Initial Release
 *
 *  This code uses the BSD License:
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are
 *  met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in
 *       the documentation and/or other materials provided with the distribution
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 * mxArray_header.h is an include file that is used to hack into the mxArray
 * variable struct for various versions of MATLAB. Also included is the
 * hacking code to get and set various fields.
 *
 */

/* mxArray header hack ------------------------------------------------ */
/* Need this to set the pr & pi data pointers directly, and also to use */
/* for CrossLink shared data copy checking.                             */

/**************************************************************************/
struct mxArray_header_2017b {
	void *RevCrossLink; /* Name of variable (R2008b-),
                           NULL (R2009a - R2010b), or
                           Reverse CrossLink pointer (R2011a+) */
	mxClassID ClassID; /*  0 = (unknown)
					       1 = cell
						   2 = struct
						   3 = logical
						   4 = char
						   5 = void
						   6 = double
						   7 = single
						   8 = int8
						   9 = uint8
					      10 = int16
						  11 = uint16
						  12 = int32
						  13 = uint32
						  14 = int64
						  15 = uint64
						  16 = function_handle
						  17 = opaque (User Defined Class indicator new classdef style)
						  18 = object (User Defined Class indicator old @directory style)
						  19 = index (deprecated)
						  20 = sparse (deprecated)
					   */
	int VariableType;  /*  0 = normal
					       1 = persistent
						   2 = global
						   3 = sub-element (field or cell)
						   4 = temporary
						   5 = (unknown)
						   6 = property of opaque class object
					   */
	mxArray *CrossLink;  /* Address of next shared-data variable in linked list */
	size_t ndim;
	unsigned int RefCount; /* Number of additional sub-elements identical to this one */
                           /* Sub-element means a struct field element or a cell element */
	unsigned int flags;
/*          R2017b and earlier:
 *              bit  0 = is scalar double full
 *              bit  2 = is empty double full (R2015a and earlier)
 *              bit  4 = is temporary
 *              bit  5 = is sparse
 *              bit  9 = is numeric
 *              bits 24 - 31 = User Bits
 */
	union {
		size_t M;        /* Row size for 2D matrices, or            */
		size_t *dims;    /* Pointer to dims array for nD > 2 arrays */
	} Mdims;
	size_t N;            /* Column size for 2D matrices = prod(dims(2:end)) */
	void *pr;            /* Pointer to real data (or pointer to cell or field elements */
	void *pi;            /* Pointer to imag data (or pointer to cell or field information */
	union {
		mwIndex *ir;        /* Pointer to row values for sparse arrays, or           */
		mxClassID ClassID;  /* User Defined Class ID (opaque new style), or          */
		char *ClassName;    /* Pointer to User Defined Class Name (object old style) */
	} irClassNameID;
	union {
		mwIndex *jc;        /* Pointer to column values for sparse arrays, or */
		mxClassID ClassID;  /* User Defined Class ID (object old style)       */
	} jcClassID;
	size_t nzmax;           /* Number of elements allocated for sparse matrix */
};

/**************************************************************************/
/* Pi data pointer removed, flags bits definition changes */
struct mxArray_header_2018a {
	void *RevCrossLink; /* Name of variable (R2008b-),
                           NULL (R2009a - R2010b), or
                           Reverse CrossLink pointer (R2011a+) */
	mxClassID ClassID; /*  0 = (unknown)
					       1 = cell
						   2 = struct
						   3 = logical
						   4 = char
						   5 = void
						   6 = double
						   7 = single
						   8 = int8
						   9 = uint8
					      10 = int16
						  11 = uint16
						  12 = int32
						  13 = uint32
						  14 = int64
						  15 = uint64
						  16 = function_handle
						  17 = opaque (User Defined Class indicator new classdef style)
						  18 = object (User Defined Class indicator old @directory style)
						  19 = index (deprecated)
						  20 = sparse (deprecated)
					   */
	int VariableType;  /*  0 = normal
					       1 = persistent
						   2 = global
						   3 = sub-element (field or cell)
						   4 = temporary
						   5 = (unknown)
						   6 = property of opaque class object
					   */
	mxArray *CrossLink;  /* Address of next shared-data variable in linked list */
	size_t ndim;
	unsigned int RefCount; /* Number of additional sub-elements identical to this one */
                           /* Sub-element means a struct field element or a cell element */
	unsigned int flags;
/*          R2018a and later:
 *	            bit   0 = is scalar double full
 *              bit   4 = is sparse
 *              bit   7 = is numeric
 *              bit  11 = is complex
 *              bits 12 - 17 = Various values related to byte usage of non-complex variables
 *              Number of Data Bytes      Bits
 *              flags bits 17-12 values for the number of data bytes
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
 *              bits 24 - 31 = User Bits
 *          So it used to be that you only had to check if Pi is not NULL to
 *          see if a numeric variable is complex. Now you have to check bit
 *          flag 11 for 1 or 0 to see if a numeric variable is complex.
 */
	union {
		size_t M;        /* Row size for 2D matrices, or            */
		size_t *dims;    /* Pointer to dims array for nD > 2 arrays */
	} Mdims;
	size_t N;            /* Column size for 2D matrices = prod(dims(2:end)) */
	void *pr;            /* Pointer to real data (or pointer to cell or field elements */
	union {
		mwIndex *ir;        /* Pointer to row values for sparse arrays, or           */
		mxClassID ClassID;  /* User Defined Class ID (opaque new style), or          */
		char *ClassName;    /* Pointer to User Defined Class Name (object old style) */
	} irClassNameID;
	union {
		mwIndex *jc;        /* Pointer to column values for sparse arrays, or */
		mxClassID ClassID;  /* User Defined Class ID (object old style)       */
	} jcClassID;
	size_t nzmax;           /* Number of elements allocated for sparse matrix */
};

/**************************************************************************/
/* CrossLink field moved */
struct mxArray_header_2019a {
	void *RevCrossLink; /* Name of variable (R2008b-),
                           NULL (R2009a - R2010b), or
                           Reverse CrossLink pointer (R2011a+) */
	mxArray *CrossLink;  /* Address of next shared-data variable in linked list */
	mxClassID ClassID; /*  0 = (unknown)
					       1 = cell
						   2 = struct
						   3 = logical
						   4 = char
						   5 = void
						   6 = double
						   7 = single
						   8 = int8
						   9 = uint8
					      10 = int16
						  11 = uint16
						  12 = int32
						  13 = uint32
						  14 = int64
						  15 = uint64
						  16 = function_handle
						  17 = opaque (User Defined Class indicator new classdef style)
						  18 = object (User Defined Class indicator old @directory style)
						  19 = index (deprecated)
						  20 = sparse (deprecated)
					   */
	int VariableType;  /*  0 = normal
					       1 = persistent
						   2 = global
						   3 = sub-element (field or cell)
						   4 = temporary
						   5 = (unknown)
						   6 = property of opaque class object
					   */
	size_t ndim;
	unsigned int RefCount; /* Number of additional sub-elements identical to this one */
                           /* Sub-element means a struct field element or a cell element */
	unsigned int flags;
/*          R2018a and later:
 *	            bit   0 = is scalar double full
 *              bit   4 = is sparse
 *              bit   7 = is numeric
 *              bit  11 = is complex
 *              bits 12 - 17 = Various values related to byte usage of non-complex variables
 *              Number of Data Bytes      Bits
 *              flags bits 17-12 values for the number of data bytes
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
 *              bits 24 - 31 = User Bits
 *          So it used to be that you only had to check if Pi is not NULL to
 *          see if a numeric variable is complex. Now you have to check bit
 *          flag 11 for 1 or 0 to see if a numeric variable is complex.
 */
	union {
		size_t M;        /* Row size for 2D matrices, or            */
		size_t *dims;    /* Pointer to dims array for nD > 2 arrays */
	} Mdims;
	size_t N;            /* Column size for 2D matrices = prod(dims(2:end)) */
	void *pr;            /* Pointer to real data (or pointer to cell or field elements */
	union {
		mwIndex *ir;        /* Pointer to row values for sparse arrays, or           */
		mxClassID ClassID;  /* User Defined Class ID (opaque new style), or          */
		char *ClassName;    /* Pointer to User Defined Class Name (object old style) */
	} irClassNameID;
	union {
		mwIndex *jc;        /* Pointer to column values for sparse arrays, or */
		mxClassID ClassID;  /* User Defined Class ID (object old style)       */
	} jcClassID;
	size_t nzmax;           /* Number of elements allocated for sparse matrix */
};

/**************************************************************************/
/* ndim field moved, flags 2D bit definition added */
struct mxArray_header_2019b {
	void *RevCrossLink; /* Name of variable (R2008b-),
                           NULL (R2009a - R2010b), or
                           Reverse CrossLink pointer (R2011a+) */
	mxArray *CrossLink;  /* Address of next shared-data variable in linked list */
	mxClassID ClassID; /*  0 = (unknown)
					       1 = cell
						   2 = struct
						   3 = logical
						   4 = char
						   5 = void
						   6 = double
						   7 = single
						   8 = int8
						   9 = uint8
					      10 = int16
						  11 = uint16
						  12 = int32
						  13 = uint32
						  14 = int64
						  15 = uint64
						  16 = function_handle
						  17 = opaque (User Defined Class indicator new classdef style)
						  18 = object (User Defined Class indicator old @directory style)
						  19 = index (deprecated)
						  20 = sparse (deprecated)
					   */
	int VariableType;  /*  0 = normal
					       1 = persistent
						   2 = global
						   3 = sub-element (field or cell)
						   4 = temporary
						   5 = (unknown)
						   6 = property of opaque class object
					   */
	unsigned int RefCount; /* Number of additional sub-elements identical to this one */
                           /* Sub-element means a struct field element or a cell element */
	unsigned int flags;
/*          R2019b and later:
 *	            bit   0 = is scalar double full
 *              bit   4 = is sparse
 *              bit   7 = is numeric
 *              bit  11 = is complex
 *              bits 12 - 17 = Various values related to byte usage of non-complex variables
 *              Number of Data Bytes      Bits
 *              flags bits 17-12 values for the number of data bytes
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
 *              bit  18 = is 2D variable
 *              bits 24 - 31 = User Bits
 *          So it used to be that you only had to check if Pi is not NULL to
 *          see if a numeric variable is complex. Now you have to check bit
 *          flag 11 for 1 or 0 to see if a numeric variable is complex.
 */
	union {
		size_t M;        /* Row size for 2D matrices, or            */
		size_t *dims;    /* Pointer to dims array for nD > 2 arrays */
	} Mdims;
	union {
		size_t N;        /* Column size for 2D matrices, or         */
		size_t ndim;     /* Number of dimensions for nD > 2 arrays  */
	} Nndim;
    size_t unknown;      /* As a result of ndim being combined above */
	void *pr;            /* Pointer to real data (or pointer to cell or field elements */
	union {
		mwIndex *ir;        /* Pointer to row values for sparse arrays, or           */
		mxClassID ClassID;  /* User Defined Class ID (opaque new style), or          */
		char *ClassName;    /* Pointer to User Defined Class Name (object old style) */
	} irClassNameID;
	union {
		mwIndex *jc;        /* Pointer to column values for sparse arrays, or */
		mxClassID ClassID;  /* User Defined Class ID (object old style)       */
	} jcClassID;
	size_t nzmax;           /* Number of elements allocated for sparse matrix */
};

/**************************************************************************/

/* Related macros */
/* Every time MATLAB changes the mxArray, will need to add a new struct definition
 * above and then modify all of the macros below to include the new definitions.
 */

/* Used to define the hacking structs */
#define MXH \
        struct mxArray_header_2017b *mxh_2017b; \
        struct mxArray_header_2018a *mxh_2018a; \
        struct mxArray_header_2019a *mxh_2019a; \
        struct mxArray_header_2019b *mxh_2019b;

/* General hacking macro for pointing into a mxArray */
#define HACK_2019b(mx) (                    (void *)(mxh_2019b = (struct mxArray_header_2019b *)(mx))                  )
#define HACK_2019a(mx) ( version<=0x2019a ? (void *)(mxh_2019a = (struct mxArray_header_2019a *)(mx)) : HACK_2019b(mx) )
#define HACK_2018a(mx) ( version<=0x2018b ? (void *)(mxh_2018a = (struct mxArray_header_2018a *)(mx)) : HACK_2019a(mx) )
#define HACK(mx)       ( version<=0x2017b ? (void *)(mxh_2017b = (struct mxArray_header_2017b *)(mx)) : HACK_2018a(mx) )

/* General hacking macro for getting a field value */
#define HACKGET_2019b(f) (                    mxh_2019b->f                    )
#define HACKGET_2019a(f) ( version<=0x2019a ? mxh_2019a->f : HACKGET_2019b(f) )
#define HACKGET_2018a(f) ( version<=0x2018b ? mxh_2018a->f : HACKGET_2019a(f) )
#define HACKGET(f)       ( version<=0x2017b ? mxh_2017b->f : HACKGET_2018a(f) )

/* Special hacking macro for getting N */
#define HACKGETN_2019b (                    mxh_2019b->Nndim.N            )
#define HACKGETN_2019a ( version<=0x2019a ? mxh_2019a->N : HACKGETN_2019b )
#define HACKGETN_2018a ( version<=0x2018b ? mxh_2018a->N : HACKGETN_2019a )
#define HACKGETN       ( version<=0x2017b ? mxh_2017b->N : HACKGETN_2018a )

/* Special hacking macro for getting ndim */
#define HACKGETndim_2019b (                    mxh_2019b->Nndim.ndim               )
#define HACKGETndim_2019a ( version<=0x2019a ? mxh_2019a->ndim : HACKGETndim_2019b )
#define HACKGETndim_2018a ( version<=0x2018b ? mxh_2018a->ndim : HACKGETndim_2019a )
#define HACKGETndim       ( version<=0x2017b ? mxh_2017b->ndim : HACKGETndim_2018a )

/* General hacking macro for setting a field value */
#define HACKSET_2019b(f,x) (                    (mxh_2019b->f = (x))                      )
#define HACKSET_2019a(f,x) ( version<=0x2019a ? (mxh_2019a->f = (x)) : HACKSET_2019b(f,x) )
#define HACKSET_2018a(f,x) ( version<=0x2018b ? (mxh_2018a->f = (x)) : HACKSET_2019a(f,x) )
#define HACKSET(f,x)       ( version<=0x2017b ? (mxh_2017b->f = (x)) : HACKSET_2018a(f,x) )
