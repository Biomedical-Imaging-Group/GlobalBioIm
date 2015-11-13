#include "mex.h"
#include "math.h"

//-------------------------------------------------------------------------
void mexFunction(int nlhs, mxArray *plhs[],
				int nrhs, const mxArray *prhs[]) {
	
	double** inVects = mxCalloc(nrhs, sizeof(double*));
	mwSize* lengths = mxCalloc(nrhs, sizeof(mwSize));
	mwSize outSize = 1;
	for (int i = 0; i < nrhs; ++i) {
		lengths[i] = fmax(mxGetN(prhs[i]), mxGetM(prhs[i]));
		inVects[i] = mxGetPr(prhs[i]);
		outSize *= lengths[i];
		//mexPrintf("%d \n", lengths[i]);
	}
	
	mwSize* lengthProds = mxCalloc(nrhs, sizeof(mwSize));
	lengthProds[0] = 1;
	for (int i = 1; i < nrhs; ++i) {
			lengthProds[i] = lengthProds[i-1] * lengths[i-1];
			//mexPrintf("%d \n", lengthProds[i]);
	}
	
	
	// output
	plhs[0] = mxCreateDoubleMatrix(nrhs, outSize, mxREAL);
	double* grid = mxGetPr( plhs[0] );
	
	for (int i = 0; i < outSize; ++i) {
		for (int dim = 0; dim < nrhs; ++dim) {
			grid[nrhs*i + dim] = inVects[dim][i / lengthProds[dim] % lengths[dim]];
		
		}
	}
	
return;	
}