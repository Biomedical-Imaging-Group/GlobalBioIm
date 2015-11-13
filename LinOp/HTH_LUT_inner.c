/*
 *
 *
 *
 *
 * */


#include "mex.h"
#include "math.h"
#include <pthread.h>

#define DIMS 3
#define NUM_THREADS 4

// container to hold thread arguments
struct tDataType {
 mwSize start;
 mwSize end;
 mwSize tInd;
};
struct tDataType tDataArray[NUM_THREADS];

double* gConvPhi;
double* yStart;
double* yStep;
double* ySize;

double* Ps;
mwSize numThetas;


double* xStart;
double* xStep;
double* xSize;

mwSize outSize[DIMS];
mwSize numPoints;
double* HTg ;

void *doWork(void *tDataIn);
				

//-------------------------------------------------------------------------
void mexFunction(int nlhs, mxArray *plhs[],
				int nrhs, const mxArray *prhs[]) {
	
	// input
	gConvPhi = mxGetPr( prhs[0] );
	
	yStart = mxGetPr( prhs[1] );
	yStep = mxGetPr( prhs[2] );
	ySize = mxGetPr( prhs[3] );
	
	Ps = mxGetPr( prhs[4] );
	numThetas = mxGetDimensions( prhs[4] )[2];
	
	xStart = mxGetPr( prhs[5] );
	xStep = mxGetPr( prhs[6] );
	xSize = mxGetPr( prhs[7] );
	
	for (int dimInd=0; dimInd < DIMS; ++dimInd) {
		outSize[dimInd] = (mwSize) xSize[dimInd];
	}
	
	// output
	plhs[0] = mxCreateNumericArray(DIMS, outSize, mxDOUBLE_CLASS, mxREAL);
	numPoints = mxGetNumberOfElements( plhs[0] );
	HTg = mxGetPr( plhs[0] );
	
	 // create threads  ------------------------------------------
  
  pthread_attr_t attr;
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE); // this is default, but more portable to define
  
  mwSize sliceSize = numPoints / NUM_THREADS;
  mexPrintf("Starting %d threads\n", NUM_THREADS);
	
	pthread_t threads[NUM_THREADS];
	void* status;
	
	for(mwSize t=0; t<NUM_THREADS; t++){
		
		tDataArray[t].start = t * sliceSize;
		
		if (t == NUM_THREADS - 1) {
			tDataArray[t].end = numPoints;
		} else {
			tDataArray[t].end =  (t+1) * sliceSize;
		}
		tDataArray[t].tInd = t;
		
		mexPrintf("starting real thread %d from %d to %d \n", t, tDataArray[t].start , tDataArray[t].end );
		pthread_create(&threads[t], &attr, doWork, (void *) &tDataArray[t]);
	}
	pthread_attr_destroy(&attr); // destroy attribute container
	
	// wait for all the threads to finish
	for(long t=0; t<NUM_THREADS; t++) {
		pthread_join(threads[t], &status);
		//mexPrintf("joined thread %d\n", t);
	}
	
	
	return;
}

void *doWork(void *tDataIn) {
	
	struct tDataType *tData;
	tData = (struct tDataType *) tDataIn;
	
	// main loop
	for (mwSize i = tData->start; i < tData->end; ++i){
		mwSize xInd[3];
		xInd[0] = i % outSize[0];
		xInd[1] = ( i / outSize[0] ) % outSize[1];
		xInd[2] = i / outSize[0] / outSize[1];
		
		// find x position
		double pos[DIMS];
		for (int dimInd = 0; dimInd < DIMS; ++dimInd) {
			pos[dimInd] = xStart[dimInd] + xInd[dimInd]*xStep[dimInd];
		}
		
		
		// loop over thetas
		for (mwSize thetaInd = 0; thetaInd < numThetas; ++thetaInd){
			mwSize PsOff = DIMS * (DIMS-1) * thetaInd;
			mwSize gConvPhiOff = 0; // only have one thing to interpolate from
			
			// project x position to a y position using Ps
			double proj[DIMS - 1] = {0, 0};
			for (mwSize row = 0; row < DIMS-1; ++row){
				for (mwSize col = 0; col < DIMS; ++col){
					proj[row] += Ps[row + col*(DIMS-1) + PsOff] * pos[col];
					
				}}
			
			// find the corresponding indices in gConvPhi
			mwSize projInds[DIMS-1];
			double projDist[DIMS-1]; // distance away from projInd, scaled [0, 1];
			for (int dimInd = 0; dimInd < DIMS-1; ++dimInd) {
				projInds[dimInd] = (proj[dimInd] - yStart[dimInd]) / yStep[dimInd];
				projDist[dimInd] = (proj[dimInd] - (yStep[dimInd]*projInds[dimInd] + yStart[dimInd])) / yStep[dimInd];
			}
			
			
			mwSize ulInd = projInds[0] + projInds[1]*((mwSize) ySize[0]) + gConvPhiOff;
			double ul,ur,ll,lr;
			
			// --check boundary conditions
			/*
			if (projInds[0] < 0 || projInds[1] < 0 ||
							projInds[0] >= ySize[0] || projInds[1] >= ySize[1] ) {
				ul = 0;
			} else {
				ul = gConvPhi[ulInd];
			}
			
			if (projInds[0] < -1 || projInds[1] < 0 ||
							projInds[0] >= ySize[0] -1 || projInds[1] >= ySize[1] ) {
				ll = 0;
			} else {
				ll = gConvPhi[ulInd + 1];
			}
			
			if (projInds[0] < 0 || projInds[1] < -1 ||
							projInds[0] >= ySize[0] || projInds[1] >= ySize[1] -1 ) {
				ur = 0;
			} else {
				ur = gConvPhi[ulInd + ((mwSize) ySize[0])];
			}
			
			
			if (projInds[0] < 0 -1 || projInds[1] < -1 ||
							projInds[0] >= ySize[0] -1 || projInds[1] >= ySize[1] -1 ) {
				lr = 0;
			} else {
				lr = gConvPhi[ulInd + ((mwSize) ySize[0]) + 1];
			}*/
			
			if (!(projInds[0] < 0 || projInds[1] < 0 ||
							projInds[0] >= ySize[0]-1 || projInds[1] >= ySize[1]-1 )) {
				ul = gConvPhi[ulInd];
				ll = gConvPhi[ulInd + 1];
				ur = gConvPhi[ulInd + ((mwSize) ySize[0])];
				lr = gConvPhi[ulInd + ((mwSize) ySize[0]) + 1];
				
				// do actual bilinear interp
				double top = ul*(1-projDist[1]) + ur*projDist[1];
				double bottom = ll*(1-projDist[1]) + lr*projDist[1];
				double val = top*(1-projDist[0]) + bottom*(projDist[0]);
				
				HTg[i] += val;
			}
			
			// print debug info
			
			/*
			 * mexPrintf("xInd: %d\n", i);
			 * mexPrintf("xSub: [%d %d %d]\n", xInd[0], xInd[1], xInd[2]);
			 * mexPrintf("xPos: [%g %g %g]\n", pos[0], pos[1], pos[2]);
			 * mexPrintf("yPos: [%g %g]\n", proj[0], proj[1]);
			 * mexPrintf("yInd: [%d %d]\n", projInds[0], projInds[1]);
			 * mexPrintf("yDist: [%g %g]\n", projDist[0], projDist[1]);
			 *
			 * mexPrintf("ul: %g\n", ul);
			 * mexPrintf("ll: %g\n", ll);
			 * mexPrintf("ur: %g\n", ur);
			 * mexPrintf("lr: %g\n", lr);
			 * mexPrintf("val: %g\n", val);
			 * mexPrintf("\n");
			 */
			
			
		} // thetas
	} // points in x
	
	return NULL;

}



