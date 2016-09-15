#include <string.h> /* needed for memcpy() */
#include <mex.h>
#include <math.h>
#include <omp.h>
#include "matrix.h"
#include "matLib2D.h"

/*mex -v mexFunctions/source/projectSpMat2x2.c 
 * CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -largeArrayDims
 * -ImexFunctions/headers/ -outdir mexFunctions/source/ */

/* In a mxArray to access the element X[i][j][z] you can do it by referring
   to the element X[i+j*dims[0]+z*dims[0]*dims[1]]
 */

/* =========================================================================
%
%  Author: stamatis.lefkimmiatis@epfl.ch
%
% =========================================================================*/


void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
  
  if (nrhs < 3 || nrhs > 4)
    mexErrMsgTxt("At least 3 and at most 4 input arguments are expected.\n");
  /*set up input arguments */
  
  int i;
  int steps[2];
  double c[1];
  
  double* X=(double *)mxGetPr(prhs[0]); // matrix input of size nx x ny x...x 3
  int  number_of_dims=mxGetNumberOfDimensions(prhs[0]);
  const mwSize *dims=mxGetDimensions(prhs[0]); 
  
  if (dims[number_of_dims-1]!=3)
    mexErrMsgTxt("The last dimension of the input should be equal to 3.\n");
  
  size_t numel_X=mxGetNumberOfElements(prhs[0]);
  int num_of_mat=numel_X/3;//number of input matrices.
  
  double p = mxGetScalar(prhs[1]);
  if (p < 1)
    mexErrMsgTxt("The order of the norm should be greater or equal to 1.\n");
  
  double* rho =(double *)mxGetPr(prhs[2]);
  size_t num_of_rho=mxGetNumberOfElements(prhs[2]);
  if (!(num_of_rho==num_of_mat || num_of_rho==1))
    mexErrMsgTxt("rho should be either a scalar or equal to the number of the input matrices.\n");
  
  double c0;
  if (nrhs < 4)
    c0=0.0;
  else
    c0=mxGetScalar(prhs[3]);
  
  //Create output arguments
  plhs[0]= mxCreateNumericArray(number_of_dims, dims, mxDOUBLE_CLASS, mxREAL); //Projected matrices.
  if (plhs[0] == NULL)
    mexErrMsgTxt("Could not create mxArray.\n");
  
  double *Xp=(double *)mxGetPr(plhs[0]); //Projected matrices
  
    
  int k;
  double tmp[3];
  double tmp2[3];
  
  if(p==1){
    #pragma omp parallel for shared(X, Xp) private(i, k, tmp, tmp2)
    for(i=0; i < num_of_mat; i++){
      for (k=0;k<3;k++)
        tmp[k]=X[i+num_of_mat*k];
      
      if (num_of_rho==1)
        projectS1(tmp2, tmp, rho[0]);
      else
        projectS1(tmp2, tmp, rho[i]);
      
      for (k=0;k<3;k++)
        Xp[i+num_of_mat*k]=tmp2[k];
    }
  }
  else if (p==2){
    #pragma omp parallel for shared(X, Xp) private(i, k, tmp, tmp2)
    for(i=0; i < num_of_mat; i++){
      for (k=0;k<3;k++)
        tmp[k]=X[i+num_of_mat*k];
      
      if (num_of_rho==1)
        projectS2(tmp2, tmp, rho[0]);
      else
        projectS2(tmp2, tmp, rho[i]);
      
      for (k=0;k<3;k++)
        Xp[i+num_of_mat*k]=tmp2[k];
    }
  }
  else if(mxIsInf(p)){
    #pragma omp parallel for shared(X, Xp) private(i, k, tmp, tmp2)
    for(i=0; i < num_of_mat; i++){
      for (k=0;k<3;k++)
        tmp[k]=X[i+num_of_mat*k];
      
      if (num_of_rho==1)
        projectSinf(tmp2, tmp, rho[0]);
      else
        projectSinf(tmp2, tmp, rho[i]);
      
      for (k=0;k<3;k++)
        Xp[i+num_of_mat*k]=tmp2[k];
    }    
  }
  
  else{
    #pragma omp parallel for shared(X, Xp) private(i, k, tmp, tmp2)
    for(i=0; i < num_of_mat; i++){
      for (k=0;k<3;k++)
        tmp[k]=X[i+num_of_mat*k];
      
      if (num_of_rho==1)
        projectSp(tmp2, tmp, rho[0], p, c, c0);
      else
        projectSp(tmp2, tmp, rho[i], p, c, c0);
      
      for (k=0;k<3;k++)
        Xp[i+num_of_mat*k]=tmp2[k];
    }    
  }
}
