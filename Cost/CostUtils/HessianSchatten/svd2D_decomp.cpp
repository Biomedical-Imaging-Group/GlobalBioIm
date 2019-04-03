#include <mex.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "matrix.h"

/***************************************************************************
  Let X be a NxMx3 matrix such that:
  
  P_mn = [X(n,m,1) X(n,m,2)
          X(n,m,2) X(n,m,3)]
          
  is a symmetric matrix. Then the present function computes the eigenvalues
  E(n,m,1) E(n,m,2) and the first eigenvector [V(n,m,1) V(n,m,2)] (the second 
  one being [V(n,m,2) -V(n,m,1)]). Hence the function outputs two matrices E
  and V of size NxMx2.
  
  Compilation:
     -linux: mex -v svd2D_decomp.cpp CFLAGS="\$CFLAGS -openmp" LDFLAGS="\$LDFLAGS -openmp" -largeArrayDims
     -mac  : mex svd2D_decomp.cpp -DUSE_BLAS_LIB -DNEW_MATLAB_BLAS -DINT_64BITS -largeArrayDims CXX=/usr/local/Cellar/gcc/6.3.0_1/bin/g++-6 CXXOPTIMFLAGS="-O3
                   -mtune=native -fomit-frame-pointer -fopenmp" LDOPTIMFLAGS=" -O " LINKLIBS="$LINKLIBS -lmwblas -lmwlapack -L"/usr/local/Cellar/gcc/6.3.0_1/lib/gcc/6" -L/ -fopenmp"
     -mac (another option): "brew install llvm" on terminal, then use the command
            mex svd2D_decomp.cpp -DUSE_BLAS_LIB -DNEW_MATLAB_BLAS -DINT_64BITS -largeArrayDims CXX="/usr/local/opt/llvm/bin/clang++"
 
  Copyright (C) 2017 
  E. Soubies emmanuel.soubies@epfl.ch
  F. Soulez

****************************************************************************/

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {

	double* X=(double *)mxGetPr(prhs[0]);                  // matrix input
	int  number_of_dims=mxGetNumberOfDimensions(prhs[0]);  // number of dimensions of the input matrix
    const mwSize *dims=mxGetDimensions(prhs[0]);           // dimension vector
    
  	if (dims[number_of_dims-1]!=3)
    	mexErrMsgTxt("The last dimension of the input should be equal to 3.\n");
    
    size_t numel_X=mxGetNumberOfElements(prhs[0]);         // number of elements in the input matrix
  	int num_of_mat=numel_X/3;                              // number lateral entries (i.e. number of svd to compute)
  		
  	//Create output arguments
  	mwSize  dimsout[3];
  	dimsout[0]=dims[0];dimsout[1]=dims[1];dimsout[2]=2;
  	plhs[0]= mxCreateNumericArray(number_of_dims, dimsout, mxDOUBLE_CLASS, mxREAL); 
  	if (plhs[0] == NULL)
    	mexErrMsgTxt("Could not create mxArray.\n"); 	
    plhs[1]= mxCreateNumericArray(number_of_dims, dimsout, mxDOUBLE_CLASS, mxREAL); 
  	if (plhs[1] == NULL)
    	mexErrMsgTxt("Could not create mxArray.\n"); 	
    double *Ye=(double *)mxGetPr(plhs[0]);    // eigenvalues
	double *Yv=(double *)mxGetPr(plhs[1]);    // eigenvectors
	
	int k,i;
	double n;
	double tmp[3];
	double E[2];
	double U[2];    
	double trace;
  	double delta;
	
    #pragma omp parallel for shared(X,Ye,Yv) private(i, k, trace, delta, n, tmp, E, U)
    for(i=0; i < num_of_mat; i++){
    	for (k=0;k<3;k++)   // get the matrix value [X(1,1) X(2,1)=X(1,2), X(2,2)]
        	tmp[k]=X[i+num_of_mat*k];
        	
        if (fabs(tmp[1]) < 1e-15){
    		E[0]=tmp[0];
    		E[1]=tmp[2];
    		U[0]=1.0;
    		U[1]=0.0;
    	}
    	else{
    		trace=tmp[0]+tmp[2];
    		delta=(tmp[0]-tmp[2])*(tmp[0]-tmp[2])+4*tmp[1]*tmp[1];
    		E[0]=0.5*(trace+sqrt(delta));
    		E[1]=0.5*(trace-sqrt(delta));
    		n=sqrt((E[0]-tmp[0])*(E[0]-tmp[0])+tmp[1]*tmp[1]);
    		U[0]=tmp[1]/n;
    		U[1]=(E[0]-tmp[0])/n;
  		}
  		
  		for (k=0;k<2;k++){  // set result
        	Ye[i+num_of_mat*k]=E[k];
  			Yv[i+num_of_mat*k]=U[k];
  		}
    }
    	
}   	
    	

