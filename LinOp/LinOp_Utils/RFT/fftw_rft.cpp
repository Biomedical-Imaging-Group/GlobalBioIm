/* computes an rft using the fftw 3 library
 ************************* fftw_rft ****************************************
 *   Copyright (C) 2018 by Rainer Heintzmann                               *
 *   heintzmann@gmail.com                                                  *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; Version 2 of the License.               *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************
 To compile:
 mex fftw_rft.cpp libfftw3f-3.lib -LC:\Users\pi96doc\Documents\Programming\Lib\libfftw3\ -IC:\Users\pi96doc\Documents\Programming\Lib\libfftw3\
 */

// #define DEBUG

#include "mex.h"
#include "matrix.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "fftw3.h"  // locally under C:\Users\pi96doc\Documents\Programming\Lib, mex compile using:
// mex fftw_rft.cpp libfftw3.a libfftw3_threads.a -LC:\Users\pi96doc\Documents\Programming\Lib\libfftw3\ -IC:\Users\pi96doc\Documents\Programming\Lib\libfftw3\ -lm

#define BaseType float

#define MAXDIM 10

// #define NOCOPY   // if defined, no copies are made and the original data is used. For this the FFTW Planning mechanism currently entirely fails.
// #define COPYCOMPACT    // if defined, also the compact dimensions will be copied to a different memory location.

static int wisdomDone=0, didImport=0;

#ifndef DEBUG
#define dbgprintf dummy
#else
#define dbgprintf printf
#endif

void dummy(char* d, ...) {return;}

fftwf_iodim dims[MAXDIM],howmany_dims[MAXDIM];
BaseType * fre=0, * falloc=0;;

// N always refers to the input data dimensions
fftwf_plan CreateRFTPlan(int NumDims, int * dirYes, int *N, BaseType * inRe, BaseType * inIm, BaseType * outRe, BaseType * outIm, int Direction, unsigned FFTWFlag) {
  fftwf_plan myPlan;
  int MaxDim=0, InputStride=1,OutputStride=1,InputSize=1,OutputSize=1,BigSize=1;
  int t,h,TDims=0,HDims=0,cutDir=0;
  for(int k = 0; k < NumDims; k++)  // remove empty dimensions
      if (N[k] > 1)
      { N[MaxDim] = N[k]; dirYes[MaxDim] = dirYes[k]; MaxDim++; TDims += (dirYes[k] > 0);}
  for(int k = NumDims; k < MAXDIM; k++)  // reset the rest
      { N[MaxDim] = 1; dirYes[MaxDim] = 0;}
  for(cutDir = 0; cutDir < MAXDIM; cutDir++)  // reset the rest
       if (dirYes[cutDir])
           break;
  NumDims=MaxDim;  // change the Maximum number of dimensions
   
  HDims = NumDims - TDims;
  t=TDims-1;
  h=HDims-1;
  InputStride=1;OutputStride=1;
  for(int k = 0; k < NumDims; k++)
  {
    if (k!=cutDir) {InputSize = N[k];OutputSize = N[k];BigSize=N[k];}
    else if (Direction > 0)
            {InputSize = N[k];BigSize=InputSize;OutputSize = (N[k]/2+1);}
        else
            {InputSize = N[k];OutputSize = ((N[k]-1)*2);BigSize=OutputSize;}

    if (dirYes[k] > 0) {
        dims[t].n = BigSize;  // has to be the large size
        dims[t].is = InputStride;
        dims[t].os = OutputStride;
        dbgprintf("for %d dims[%d].n =%d .is=%d .os=%d \n",k,t,dims[t].n,dims[t].is,dims[t].os);
        t--;
    } else {
        howmany_dims[h].n = InputSize;  // number of repititions
        howmany_dims[h].is = InputStride;
        howmany_dims[h].os = OutputStride;
        dbgprintf("for %d howmany_dims[%d].n =%d .is=%d .os=%d \n",k,h,howmany_dims[h].n,howmany_dims[h].is,howmany_dims[h].os);
        h--;
    }
    InputStride *= InputSize;OutputStride *= OutputSize;
    dbgprintf("Direction %d, cutDir %d, N[%d]=%d, InputSize %d, OutputSize %d, InputStride %d, OutputStride %d\n",Direction,cutDir,k,N[k],InputSize,OutputSize,InputStride,OutputStride);
  }
  
  if (Direction>0)
      if (HDims <= 0)
        myPlan=fftwf_plan_guru_split_dft_r2c(TDims, dims, 0, NULL, inRe, outRe, outIm, FFTWFlag); // FFTW_MEASURE, FFTW_ESTIMATE
      else
        myPlan=fftwf_plan_guru_split_dft_r2c(TDims, dims, HDims, howmany_dims, inRe, outRe, outIm, FFTWFlag); // FFTW_MEASURE
  else
      if (HDims <= 0)
        myPlan=fftwf_plan_guru_split_dft_c2r(TDims, dims, 0, NULL, inRe, inIm, outRe, FFTW_ESTIMATE);  // howmany can be 1 and iodim size of that direction
      else
        myPlan=fftwf_plan_guru_split_dft_c2r(TDims, dims, HDims, howmany_dims, inRe, inIm, outRe, FFTW_ESTIMATE);  // howmany can be 1 and iodim size of that direction
              
  return myPlan;
}


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {  // arguments are:  real-valued array, direction, numthreads

#define B_OUT     plhs[0]

  int k, numCPU;
  size_t NumElIn = 1, NumElOut=1;
  const mwSize *N;
  BaseType *pr, *pi, *pr2, *pi2,*prB,*piB, * fim=0, * InRe=0, * InIm=0;
  static long MatLeng = 0;
  fftwf_plan myPlan,myPlan2;
  int NumDims=1, Direction = 1, cutDir=1, status;
  int InDimensions[MAXDIM], dirYes[MAXDIM], YesDims,doWisdom=0;
  mwSize OutDimensions[MAXDIM];
  double * YesData=0;
  const int FNLgth=10000;
  char WisdomFileName[FNLgth];  
  if (nrhs != 4 && nrhs != 5) {
      mexErrMsgTxt("Four or five input argument required (data, direction, transformDirVector, number of threads, WisdomFileName).");
  }

  if (!mxIsSingle(prhs[0])) {
      mexErrMsgTxt( "Array must be single");
  }

  if (!mxIsDouble(prhs[2])) {
      mexErrMsgTxt( "Transform dimension vector must be double");
  }

  if (!mxIsDouble(prhs[1])) {
      mexErrMsgTxt( "Direction Flag must be double");
  }
  Direction = (int) mxGetScalar(prhs[1]);
  if (Direction!=1 && Direction!=-1)
      if (Direction==2)
          {doWisdom=1;Direction=1;}
      else if (Direction==-2)
          {doWisdom=1;Direction=-1;}
      else
      mexErrMsgTxt( "Direction needs to be 1 or -1 (or 2 or -2 for FFTW_MEASURE).");
          
  if (Direction==1)
  { if (mxIsComplex(prhs[0]))
        mexErrMsgTxt( "Input array for rft must not be complex");}
#ifdef NOCOPY  // in this case currently no empty dummy array is generated.
   else
     if (!mxIsComplex(prhs[0]))
         mexErrMsgTxt( "Input array irft must be complex");
#endif

  if (!mxIsDouble(prhs[2])) {
      mexErrMsgTxt("NumThreads must be double");
  }

  numCPU = (int) mxGetScalar(prhs[3]);
  if (numCPU > 8) {
      mexErrMsgTxt("NumThreads > 8 was requested but only 8 threads are supported");
  }
  
  WisdomFileName[0]=0;
  if (nrhs==5) {
      status=mxGetString(prhs[4], WisdomFileName, FNLgth-1);
      if (status != 0)
          mexErrMsgTxt("Wisdomfilename has to be a string.");
      }

  NumDims = (int) mxGetNumberOfDimensions(prhs[0]);
  if (NumDims >= MAXDIM) {
      mexErrMsgTxt("The input array has more than maximally allowed number of dimensions");
  }

  if (NumDims >= MAXDIM)
     mexErrMsgTxt("Number of dimensions above maximally allowed dimensions.");

  N = mxGetDimensions(prhs[0]);  // Input Dimensions

  YesDims = (int) mxGetNumberOfElements(prhs[2]);
  if (YesDims >= MAXDIM)
      mexErrMsgTxt("The TransformDim vector has more than maximally allowed number of dimensions");

  if (YesDims != NumDims)
      mexErrMsgTxt( "The TransformDim vector must agree to number of dimensions of data");

  YesData = mxGetPr(prhs[2]);  // Input Dimensions
  for(k=0;k<MAXDIM;k++) {
      if (k < NumDims)
        dirYes[k]=(int) YesData[k];
      else
        dirYes[k]= 0;
  }
  
//  OutDimensions = (int*) mxMalloc( sizeof(int) * NumDims);
  NumElOut=1;
  cutDir=-1; 
  for(k=0;k<NumDims;k++) {  // find the cut direction
      if (N[k]>1 && dirYes[k] == 1)
        {cutDir=k;break;}
  } 
  dbgprintf("Cut direction is %d\n",cutDir);
  if (cutDir < 0)
     mexErrMsgTxt( "No transform direction found or transform direction is singleton in size.");
  
  for(k=0;k<NumDims;k++) {
    if (k!=cutDir)
        OutDimensions[k] = (int) N[k];
    else
        if (Direction>0)
            OutDimensions[k] = (((int) N[k])/2)+1;
        else
            OutDimensions[k] = (((int) N[k])-1)*2;
    InDimensions[k] = (int) N[k];
    NumElIn *= InDimensions[k];
    NumElOut *= OutDimensions[k];
    dbgprintf("Direction=%d, Dimension N[%d]=%d was %d is %d\n",Direction,k,(int) N[k],InDimensions[k],OutDimensions[k]);
  }

  //B_OUT = mxCreateNumericArray(NumDims, N, mxDOUBLE_CLASS, mxCOMPLEX);
  if (Direction > 0) {    
    BaseType *xre = (BaseType* ) mxMalloc( sizeof(BaseType) * NumElOut);
    BaseType *xim = (BaseType* ) mxMalloc( sizeof(BaseType) * NumElOut);
    B_OUT  = mxCreateNumericMatrix(0, 0, mxSINGLE_CLASS, mxCOMPLEX);  // make the output array
    mxSetDimensions(B_OUT , OutDimensions, NumDims);
    
    mxSetData(B_OUT , xre);  
    mxSetImagData(B_OUT , xim);

    pr = (BaseType *) mxGetPr(prhs[0]);
    pr2 = (BaseType* ) mxGetPr(B_OUT);
    pi2 = (BaseType* ) mxGetPi(B_OUT);

#ifndef NOCOPY
#ifdef COPYCOMPACT    
    falloc = (BaseType *) fftwf_malloc(sizeof(BaseType) * (NumElOut * 2 + NumElIn));  // (BaseType* ) mxMalloc( sizeof(BaseType) * NumElOut)
    fim = falloc + NumElOut; // (BaseType* ) mxMalloc( sizeof(BaseType) * NumElOut)
    InRe = falloc + 2*NumElOut;
    fre=falloc;
    memcpy(InRe,pr,sizeof(BaseType) * NumElIn);   // restore the input data    
#else
    falloc = (BaseType *) fftwf_malloc(sizeof(BaseType) * (NumElOut * 2));  // (BaseType* ) mxMalloc( sizeof(BaseType) * NumElOut)
    fre=falloc;
    fim = fre + NumElOut; // (BaseType* ) mxMalloc( sizeof(BaseType) * NumElOut)
    InRe=pr;
#endif    
#else
    fre=pr2;
    fim=pi2;
    InRe=pr;
#endif
  }
  else {
    // BaseType *xre = (BaseType *) fftwf_malloc(sizeof(BaseType) * NumElOut);  // (BaseType* ) mxMalloc( sizeof(BaseType) * NumElOut)
    BaseType *xre = (BaseType* ) mxMalloc( sizeof(BaseType) * NumElOut);
    B_OUT  = mxCreateNumericMatrix(0, 0, mxSINGLE_CLASS, mxREAL);  // make the output array
    mxSetDimensions(B_OUT , OutDimensions, NumDims);
    mxSetData(B_OUT , xre);  
    pr = (BaseType *) mxGetPr(prhs[0]);
    if (mxIsComplex(prhs[0]))
        pi = (BaseType *) mxGetPi(prhs[0]);
    pr2 = (BaseType* ) mxGetPr(B_OUT);
#ifndef NOCOPY
#ifdef COPYCOMPACT    
    falloc = (BaseType *) fftwf_malloc(sizeof(BaseType) * (NumElOut + 2* NumElIn));  // (BaseType* ) mxMalloc( sizeof(BaseType) * NumElOut)
    fim=0;
    InRe = falloc + NumElOut;
    InIm = falloc + NumElOut + NumElIn;
    memcpy(InRe,pr,sizeof(BaseType) * NumElIn);   // copy the input data    
    if (mxIsComplex(prhs[0]))
        memcpy(InIm,pi,sizeof(BaseType) * NumElIn);   // copy the input data
    else
        memset(InIm,0,sizeof(BaseType) * NumElIn);   // clear the imaginary part  (CAREFUL! This could go wrong, if the double zero does not correspond to char 0)
    fre=falloc;
#else
    falloc = (BaseType *) fftwf_malloc(sizeof(BaseType) * (2* NumElIn));  // (BaseType* ) mxMalloc( sizeof(BaseType) * NumElOut)
    fim=0;
    InRe = falloc;
    InIm = falloc + NumElIn;
    memcpy(InRe,pr,sizeof(BaseType) * NumElIn);   // copy the input data    
    if (mxIsComplex(prhs[0]))
        memcpy(InIm,pi,sizeof(BaseType) * NumElIn);   // copy the input data
    else
        memset(InIm,0,sizeof(BaseType) * NumElIn);   // clear the imaginary part  (CAREFUL! This could go wrong, if the double zero does not correspond to char 0)
    fre=pr2;
#endif
#else
    fre=pr2;
    fim=0;
    InRe=pr;
    InIm=pi;
#endif
  }

  fftwf_plan_with_nthreads(numCPU);

  if (!wisdomDone || doWisdom)  { 
    mxArray * myCpy=0;
    if (!didImport)
        { fftwf_import_wisdom_from_filename(WisdomFileName); didImport=1;printf("WARNING: FFTW-Wisdom was not yet imported. Importing the file %s as defined by the global FFTW_WisdomFilename\n",WisdomFileName);
          fftwf_init_threads();}
    myPlan = CreateRFTPlan(NumDims, dirYes, InDimensions, InRe, InIm, fre, fim, Direction, FFTW_WISDOM_ONLY);  // FFTW_MEASURE, FFTW_EXHAUSTIVE
    if (myPlan == 0) {  // Here the data first has to be saved and then restored.
        int doBackup=0;
#ifdef NOCOPY        
        doBackup=1;
#else
#ifndef COPYCOMPACT
        if (Direction > 0)
            doBackup=1;
#endif  
#endif
        if (doBackup) {
            myCpy=mxDuplicateArray(prhs[0]);
            prB = (BaseType *) mxGetPr(myCpy);
        }
        printf("WARNING: No FFTW-Wisdom exists for this plan size. Estimating with FFTW_PATIENT and saving to global FFTW_WisdomFilename=%s.\n",WisdomFileName);
        myPlan = CreateRFTPlan(NumDims, dirYes, InDimensions, InRe, InIm, fre, fim, Direction, FFTW_MEASURE);  // FFTW_MEASURE, FFTW_PATIENT, FFTW_EXHAUSTIVE
        wisdomDone=1; doWisdom=0;
        if (Direction > 0) {
            if (doBackup) {
                memcpy(pr,prB,sizeof(BaseType) * NumElIn);// restore the input data
                if (mxIsComplex(prhs[0]))
                    piB = (BaseType *) mxGetPi(myCpy);
            }
            else
                memcpy(InRe,pr,sizeof(BaseType) * NumElIn);// restore the input data
        } else {
            if (doBackup) {
                memcpy(pr,prB,sizeof(BaseType) * NumElIn); // restore the input data
                if (mxIsComplex(prhs[0]))
                    {piB = (BaseType *) mxGetPi(myCpy);
                    memcpy(pi,piB,sizeof(BaseType) * NumElIn);}
            } else {
                memcpy(InRe,pr,sizeof(BaseType) * NumElIn);// restore the input data
                if (mxIsComplex(prhs[0]))
                    memcpy(InIm,pi,sizeof(BaseType) * NumElIn);   // copy the input data
                else
                    memset(InIm,0,sizeof(BaseType) * NumElIn);   // clear the imaginary part  (CAREFUL! This could go wrong, if the double zero does not correspond to char 0)
            }
        }
        fftwf_export_wisdom_to_filename(WisdomFileName);  // save the new wisdom
        
//        fftwf_cleanup_threads();
//        fftwf_init_threads();        
        myPlan2 = CreateRFTPlan(NumDims, dirYes, InDimensions, InRe, InIm, fre, fim, Direction, FFTW_WISDOM_ONLY);  // FFTW_MEASURE, FFTW_EXHAUSTIVE
        if (myPlan2 == 0) 
            mexErrMsgTxt("Even the same plan does not work again.");
        
        if (doBackup)
            mxDestroyArray(myCpy);
    }  // if (myPlan==0)
    // mexErrMsgTxt("The first call to rfftw was used to generate wisdom. Please call again.");
  }  else  // no need to use wisdom
    myPlan = CreateRFTPlan(NumDims, dirYes, InDimensions, InRe, InIm, fre, fim, Direction, FFTW_ESTIMATE);  // FFTW_MEASURE, 

  if(!myPlan)
    { mexErrMsgTxt("Real to half complex FFT using FFTW failed to create a plan.");return;}

  // fftwf_execute_split_dft(PlanForward, pr, pi, pr2, pi2);
  fftwf_execute(myPlan);
  
   if (Direction > 0)
     {memcpy(pr2,fre,sizeof(BaseType) * NumElOut);memcpy(pi2,fim,sizeof(BaseType) * NumElOut);}  // fftwf_free(fim); copy results to Matlab array
   else
     {memcpy(pr2,fre,sizeof(BaseType) * NumElOut);} // fftwf_free(fre);
    
  fftwf_destroy_plan(myPlan);
//  fftwf_cleanup_threads();   // breaks the WISDOM accumulation
#ifndef NOCOPY
    if (falloc) fftwf_free(falloc);
#endif
//  if (fim) fftwf_free(fim);// This should NOT be called, since there was only one allocation operation
  return;
}
