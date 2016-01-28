/* 	
 function m_vmlmb_first
 
  [csave, isave, dsave] = m_vmlmb_first(n, m, fatol, frtol, sftol, sgtol, ...
                                        sxtol, epsilon, costheta);
                                      
  MexFile to call from Matlab the op_vmlmb_first function of the OptimPack library
  See test_optim_vmlmb.m and test_optim_vmlmb_const.m scripts for test examples
 
  Compilation :
        mex m_vmlmb_first.c -loptimpack
  
  September 2013
  IRAP, Observatoire Midi-Pyrénées, Toulouse
  Hervé Carfantan
*/

#include "mex.h"
#include "optimpack.h"

void mexFunction( int nlhs1, mxArray *plhs[],    
                  int nrhs, const mxArray *prhs[] )
{ 
    double fatol, frtol, sftol, sgtol, sxtol, epsilon, costheta;
    int n, m;
    mwSize dims[]={1,12}; /* size of isave */
    double *dsave;
    char *csave; 
    int *isave;

    /* Control of the number of inputs and outputs */ 
    if (nrhs > 9)
        mexErrMsgTxt("9 input argument required."); 
    else if (nlhs1 !=3) 
        mexErrMsgTxt("3 output argument required.");
    
  
    /* Control of the inputs*/
    if (!mxIsNumeric(prhs[0]) || !mxIsScalar(prhs[0]) || !mxIsNumeric(prhs[1]) || !mxIsScalar(prhs[1])
        || !mxIsNumeric(prhs[2]) || !mxIsScalar(prhs[2]) || !mxIsNumeric(prhs[3]) || !mxIsScalar(prhs[3])
        || !mxIsNumeric(prhs[4]) || !mxIsScalar(prhs[4]) || !mxIsNumeric(prhs[5]) || !mxIsScalar(prhs[5])
        || !mxIsNumeric(prhs[6]) || !mxIsScalar(prhs[6]) || !mxIsNumeric(prhs[7]) || !mxIsScalar(prhs[7])
        || !mxIsNumeric(prhs[8]) || !mxIsScalar(prhs[8]) ) 
        mexErrMsgTxt("The 9 input must be scalar."); 

    if ( (double )(int )mxGetScalar(prhs[0]) != mxGetScalar(prhs[0]) )
        mexErrMsgTxt("The first input n must be an integer.");
    if ( (double )(int )mxGetScalar(prhs[1]) != mxGetScalar(prhs[1]) )
        mexErrMsgTxt("The second input m must be an integer.");
    n        = (int )mxGetScalar(prhs[0]);
    m        = (int )mxGetScalar(prhs[1]);
    fatol    = mxGetScalar(prhs[2]);
    frtol    = mxGetScalar(prhs[3]);
    sftol    = mxGetScalar(prhs[4]);
    sgtol    = mxGetScalar(prhs[5]);
    sxtol    = mxGetScalar(prhs[6]);
    epsilon  = mxGetScalar(prhs[7]);
    costheta = mxGetScalar(prhs[8]);
    
    /* Display the inputs */
    /*
    mexPrintf("n = %d, m = %d, fatol = %.3e, frtol = %.3e, sftol = %.3e\n",n ,m, fatol, frtol, sftol);
    mexPrintf("sgtol = %.3e, sxtol = %.3e, epsilon = %.3e, costheta = %f\n",sgtol, sxtol, epsilon, costheta);
    */
    
    /* Create the output arguments */
    
    /* Argument csave */
    csave = (char *)mxCalloc(128, sizeof(char));
       
    /* Argument isave */
    plhs[1] = mxCreateNumericArray(2,dims,mxINT32_CLASS,mxREAL); /* For Linux 64 bit, int codded with 32 bit */
    isave   = (op_integer_t *)mxGetPr(plhs[1]);                 
    
    /* Argument dsave */
    plhs[2] = mxCreateDoubleMatrix(OP_VMLMB_DSAVE_NUMBER(n,m),1, mxREAL);
    dsave   = mxGetPr(plhs[2]);
   
    /* Call to the op_vmlmb_first function */
    op_vmlmb_first(n, m, fatol, frtol, sftol, sgtol, sxtol, epsilon, costheta, csave, isave, dsave); 

    /* Allocation for csave and free memory */
    plhs[0] = mxCreateString(csave);
    mxFree(csave);

    /* Display the ouputs */
    /*
    mexPrintf("csave = |%s|\n", csave);
    mexPrintf("isave[");
       for (k=0;k<12;k++)
           mexPrintf(" %ld,",isave[k]);
    mexPrintf("];\n");
    mexPrintf("dsave[");    
        for (k=0;k<d;k++)
           mexPrintf(" %.3e,",dsave[k]);
    mexPrintf("];\n");
    */
}
