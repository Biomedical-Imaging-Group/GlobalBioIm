/* 	
  function m_vmlmb_next

  [task, csave] =  m_op_vmlmb_next( x, f, g, isave, dsave);

  MexFile to call from Matlab the op_vmlmb_next function of the OptimPack library
  See test_optim_vmlmb.m and test_optim_vmlmb_const.m scripts for test examples

  WARNING : 
       To be coherent with the OptimPack library in terms of limited memory 
       used, but contrarily to the usual Matlab functions, the isave and 
       dsave variables are passed by reference in the m_vmlmb_next function, 
       so they might both be considered input and output parameters 

  Compilation :
        mex m_vmlmb_next.c -loptimpack
  
  September 2013
  IRAP, Observatoire Midi-Pyrénées, Toulouse
  Hervé Carfantan
*/

#include "mex.h"
#include "optimpack.h"

#define INDEX_OF_M 4
#define INDEX_OF_N 5

void mexFunction( int nlhs, mxArray *plhs[],    
                  int nrhs, const mxArray *prhs[] )
{ 
    double *x, *f, *g, *dsave;
    char *csave;
    int *task;
    op_logical_t *isave, *active;
    int n, m, k ;
    mwSize dims[]={1,1};
    unsigned int nl, nc;
    size_t len_csave;

    /* Control of the number of inputs and outputs */ 
    if (nrhs != 6) {
        mexErrMsgTxt("6 input argument are required."); } 
    else if (nlhs != 2) {
        mexErrMsgTxt("2 output argument allowed."); }
    if ( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]))
        mexErrMsgTxt("The first input x must be a real valued vector.");
    x = mxGetPr(prhs[0]);
    /* mexPrintf("x = [");
    for (k=0;k<mxGetM(prhs[0])*mxGetN(prhs[0]);k++)
        mexPrintf("%.3f, ",x[k]);
    mexPrintf("]\n");
    */
        
    /* Control the input f */
    if (!mxIsDouble(prhs[1]) || !mxIsScalar(prhs[1])  || mxIsComplex(prhs[1]) )
        mexErrMsgTxt("The second input f must be a real valued scalar.");
    f = mxGetPr(prhs[1]);
    /* mexPrintf("f = %.3f",f[0]); */
    
    /* Control the input g */
    if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) )
        mexErrMsgTxt("The third input g must be a real valued vector");
    g = mxGetPr(prhs[2]);

    /* Control the input active */
    if ( mxGetM(prhs[3])*mxGetN(prhs[3])==0 )
        active = NULL;
    else {
        if ( !mxIsNumeric(prhs[3]) || (mxGetClassID(prhs[3])!=mxINT32_CLASS) || mxIsComplex(prhs[3]) )
            mexErrMsgTxt("The fourth input active must be an integer vector or an empty vector");
        active = (op_integer_t *)mxGetPr(prhs[3]);
    }
    
    /* Control the input isave */
    if ( !mxIsNumeric(prhs[4]) || (mxGetClassID(prhs[4])!=mxINT32_CLASS) || mxIsComplex(prhs[4]) ||
         (mxGetM(prhs[4])*mxGetN(prhs[4])!=OP_VMLMB_ISAVE_NUMBER) )
        mexErrMsgTxt("The fifth input isave must be an integer vector of length 12");
    isave = (op_integer_t *)mxGetPr(prhs[4]);  /* op_integer_t are 16 bits integers */
    n = isave[INDEX_OF_N];
    m = isave[INDEX_OF_M];
    /* mexPrintf("n = %d, m = %d\n",n,m);
    mexPrintf("isave = [");
    for (k=0;k<12;k++)
        mexPrintf("%d, ",isave[k]);
    mexPrintf("]\n");
    */
            
    /* Control the input dsave */
    if ( !mxIsNumeric(prhs[5]) || !mxIsDouble(prhs[5]) || mxIsComplex(prhs[5]) )
        mexErrMsgTxt("The sixth input dsave must be a vector of real valued double");
    dsave = mxGetPr(prhs[5]);  
    /* printf("dsave = [");
    for (k=0;k<OP_VMLMB_DSAVE_NUMBER(n,m);k++)
        mexPrintf("%.3f, ",dsave[k]);
    mexPrintf("]\n");
    */

    /* Control the size of the vectors */
    if (mxGetM(prhs[0])*mxGetN(prhs[0])!=n)
        mexErrMsgTxt("Incorrect dimension for the first input vector x");
    if (mxGetM(prhs[2])*mxGetN(prhs[2])!=n)
        mexErrMsgTxt("Incorrect dimension for the third input vector g");
    if (mxGetM(prhs[3])*mxGetN(prhs[3])!=0 && mxGetM(prhs[3])*mxGetN(prhs[3])!=n)
        mexErrMsgTxt("Incorrect dimension for the fourth input vector active");
    if (mxGetM(prhs[5])*mxGetN(prhs[5])!=OP_VMLMB_DSAVE_NUMBER(n,m))
        mexErrMsgTxt("Incorrect dimension for input vector dsave");

    /* Create the outputs */

    /* Ouput task */
    plhs[0] = mxCreateNumericArray(2,dims,mxINT32_CLASS,mxREAL); /* on Linux 64 bits, int coded with 32 bits */
    task = (op_integer_t *)mxGetPr(plhs[0]);  /* op_integer_t are 16 bits integers */

    /* Ouput csave */
    csave= (char *)mxCalloc(128, sizeof(char));


     
    /* Display the inputs */
    /*
    mexPrintf("x = [%f %f], f = %f, g = [%f %f]\n", x[0], x[1], f, g[0], g[1]); 
    mexPrintf("isave[");
       for (k=0;k<12;k++)
           mexPrintf(" %ld,",isave[k]);
    mexPrintf("];\n");
    mexPrintf("dsave[");    
        for (k=0;k<d;k++)
           mexPrintf(" %.3e,",dsave[k]);
    mexPrintf("];\n");
    */
    
    /* Call to the op_vmlmb_next function */
    task[0] = op_vmlmb_next(x, f, g, active, NULL, csave, isave, dsave);

    /* Allocation for csave and free memory */
    plhs[1] = mxCreateString(csave);
    mxFree(csave);

    /* Display the ouputs */
    /*
    mexPrintf("task = %d\n csave = |%s|\n", task[0],csave); 
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
