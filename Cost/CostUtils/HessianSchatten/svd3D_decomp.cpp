#include <mex.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "matrix.h"
#include "matLib3D.h"

/***************************************************************************
  Let X be a NxMxKx6 matrix such that:
  
  P_mn = [X(n,m,k,1) X(n,m,k,2) X(n,m,k,3)
          X(n,m,k,2) X(n,m,k,4) X(n,m,k,5)
          X(n,m,k,3) X(n,m,k,5) X(n,m,k,6)] 
          
  is a symmetric matrix. Then the present function computes the eigenvalues
  E(n,m,k,1) E(n,m,k,2) E(n,m,k,3) and the eigenvector 
          V1 = [V(n,m,k,1) V(n,m,k,2)  V(n,m,k,3)] 
          V2 = [V(n,m,k,4) V(n,m,k,5)  V(n,m,k,6)] 
          V2 = [V(n,m,k,5) V(n,m,k,8)  V(n,m,k,9)]  
  Hence the function outputs two matrices E of size NxMxKx3 and V of size NxMxKx9.
  
  Compilation:
     -linux: mex svd2D_decomp.cpp CFLAGS="\$CFLAGS -openmp" LDFLAGS="\$LDFLAGS -openmp" -largeArrayDims
     -mac  : mex svd3D_decomp.cpp -DUSE_BLAS_LIB -DNEW_MATLAB_BLAS -DINT_64BITS -largeArrayDims CXX=/usr/local/Cellar/gcc/6.3.0_1/bin/g++-6 CXXOPTIMFLAGS="-O3
                   -mtune=native -fomit-frame-pointer -fopenmp" LDOPTIMFLAGS=" -O " LINKLIBS="$LINKLIBS -lmwblas -lmwlapack -L"/usr/local/Cellar/gcc/6.3.0_1/lib/gcc/6" -L/ -fopenmp"
  
  Copyright (C) 2017 E. Soubies emmanuel.soubies@epfl.ch

****************************************************************************/


void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	double* X=(double *)mxGetPr(prhs[0]);                  // matrix input
	int  number_of_dims=mxGetNumberOfDimensions(prhs[0]);  // number of dimensions of the input matrix
    const mwSize *dims=mxGetDimensions(prhs[0]);           // dimension vector
    
  	if (dims[number_of_dims-1]!=6)
    	mexErrMsgTxt("The last dimension of the input should be equal to 6.\n");
    
    size_t numel_X=mxGetNumberOfElements(prhs[0]);         // number of elements in the input matrix
  	int num_of_mat=numel_X/6;                              // number lateral entries (i.e. number of svd to compute)
  		
  	//Create output arguments
  	mwSize  dimsout[4];
  	dimsout[0]=dims[0];dimsout[1]=dims[1];dimsout[2]=dims[2];dimsout[3]=3;
  	plhs[0]= mxCreateNumericArray(number_of_dims, dimsout, mxDOUBLE_CLASS, mxREAL); 
  	if (plhs[0] == NULL)
    	mexErrMsgTxt("Could not create mxArray.\n"); 
    dimsout[3]=9;
    plhs[1]= mxCreateNumericArray(number_of_dims, dimsout, mxDOUBLE_CLASS, mxREAL); 
  	if (plhs[1] == NULL)
    	mexErrMsgTxt("Could not create mxArray.\n"); 	
    double *Ye=(double *)mxGetPr(plhs[0]);    // eigenvalues
	double *Yv=(double *)mxGetPr(plhs[1]);    // eigenvectors 
    
	int k,i;	
	double E[3];
	double D[3];   
    double V[9];
	
    #pragma omp parallel for shared(X, Ye,Yv) private(i,k,V,E,D)
    for(i=0; i < num_of_mat; i++){
        for (k=0;k<3;k++)   // get the matrix value [X(1,1) X(2,1)=X(1,2), X(2,2)]
            V[k]=X[i+num_of_mat*k];
        for (k = 4; k < 6; k++)
            V[k] = X[i+num_of_mat*(k-1)];
        for (k = 7; k < 9; k++)
            V[k] = X[i+num_of_mat*(k-3)];
        V[3]=X[i+num_of_mat];
        V[6]=X[i+num_of_mat*2];
        
        tred2(V, D, E);
        tql2(V, D, E);
        
  		for (k=0;k<3;k++)  // set result
        	Ye[i+num_of_mat*k]=D[k];
        
        for (k=0;k<9;k++){
            Yv[i+num_of_mat*k]=V[k];
        }       
    }   	
}   	
    	

