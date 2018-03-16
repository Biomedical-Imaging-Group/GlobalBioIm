#include <mex.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "matrix.h"
#include "matLib3D.h"

/***************************************************************************

  Reconstruct X from E and V obtained by svd3D_decomp
  
  Compilation:
     -linux: mex svd3D_recomp.cpp CFLAGS="\$CFLAGS -openmp" LDFLAGS="\$LDFLAGS -openmp" -largeArrayDims
     -mac  : mex svd3D_recomp.cpp -DUSE_BLAS_LIB -DNEW_MATLAB_BLAS -DINT_64BITS -largeArrayDims CXX=/usr/local/Cellar/gcc/6.3.0_1/bin/g++-6 CXXOPTIMFLAGS="-O3
                   -mtune=native -fomit-frame-pointer -fopenmp" LDOPTIMFLAGS=" -O " LINKLIBS="$LINKLIBS -lmwblas -lmwlapack -L"/usr/local/Cellar/gcc/6.3.0_1/lib/gcc/6" -L/ -fopenmp"
 
  Copyright (C) 2017 E. Soubies emmanuel.soubies@epfl.ch

****************************************************************************/

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {

	double* E=(double *)mxGetPr(prhs[0]);                       // input eigenvalues
	double* V=(double *)mxGetPr(prhs[1]);                       // input eigenvectors
	int  number_of_dimsE=mxGetNumberOfDimensions(prhs[0]);      // number of dimensions of E
	int  number_of_dimsV=mxGetNumberOfDimensions(prhs[1]);      // number of dimensions of V
    const mwSize *dimsE=mxGetDimensions(prhs[0]);               // dimension vector E
    const mwSize *dimsV=mxGetDimensions(prhs[1]);               // dimension vector V
    
    if ((number_of_dimsE!=4) || (number_of_dimsV!=4)){
    	mexErrMsgTxt("The inputs should be 3D matrices.\n");
        }
  	if ((dimsE[number_of_dimsE-1]!=3) || (dimsV[number_of_dimsV-1]!=9)){
    	mexErrMsgTxt("The last dimension of E inputs should be equal to 3 and of V should be equal to 9.\n");
        }
    if ((dimsE[0]!=dimsV[0]) || (dimsE[1]!=dimsV[1]) || (dimsE[2]!=dimsV[2])){
    	mexErrMsgTxt("The inputs should have the same size in the three first dimensions.\n");
        }
    
    size_t numel_X=mxGetNumberOfElements(prhs[0]);         // number of elements in the input matrix
  	int num_of_mat=numel_X/3;                              // number lateral entries (i.e. number of svd to compute)
  		
  	//Create output arguments
  	mwSize  dimsout[4];
  	dimsout[0]=dimsE[0];dimsout[1]=dimsE[1];dimsout[2]=dimsE[2];dimsout[3]=6;
  	plhs[0]= mxCreateNumericArray(number_of_dimsE, dimsout, mxDOUBLE_CLASS, mxREAL); 
  	if (plhs[0] == NULL)
    	mexErrMsgTxt("Could not create mxArray.\n"); 		
    double *Y=(double *)mxGetPr(plhs[0]);    // reconstructed matrix
	
	int k,i;
	double ee[3];
	double vv[9];
	double tmp[6];  
	
    #pragma omp parallel for shared(E,V,Y) private(i, k,ee,vv,tmp)
    for(i=0; i < num_of_mat; i++){
    	for (k=0;k<3;k++)   // get the eigenvalues 
        	ee[k]=E[i+num_of_mat*k];
        for (k=0;k<9;k++)   // get the eigenvectors
        	vv[k]=V[i+num_of_mat*k];

        eigen3x3SymRec(tmp,vv,ee);
  		for (k=0;k<6;k++){  // set result
        	Y[i+num_of_mat*k]=tmp[k];
  		}
    }
    	
}   	
