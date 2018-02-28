/* 
 *  m_opl_vmlmb_iterate.c
 *
 * function m_opl_vmlmb_iterate
 *
 *	Definitions for optimization routines implemented in OptimPack
 *	library.
 *
 *-----------------------------------------------------------------------------
 *
 *      Copyright (c) 2018, Ferreol SOULEZ.
 *
 *	This file is part of OptimPack.
 *
 *	OptimPack is  free software; you can redistribute  it and/or modify
 *	it under the  terms of the GNU General  Public License as published
 *	by the Free  Software Foundation; either version 2  of the License,
 *	or (at your option) any later version.
 *
 *	OptimPack is  distributed in the hope  that it will  be useful, but
 *	WITHOUT  ANY  WARRANTY;  without   even  the  implied  warranty  of
 *	MERCHANTABILITY or  FITNESS FOR A PARTICULAR PURPOSE.   See the GNU
 *	General Public License for more details.
 *
 *	You should have  received a copy of the  GNU General Public License
 *	along with OptimPack (file  "LICENSE" in the top source directory);
 *	if  not, write  to the  Free Software  Foundation, Inc.,  59 Temple
 *	Place, Suite 330, Boston, MA 02111-1307 USA
ws = opl_vmlmb_create(dims, mem, fmin=fmin,fatol=fatol, frtol=frtol,sftol=sftol, sgtol=sgtol, sxtol=sxtol);
*/

#include <errno.h>

#include "mex.h"
#include "optimpacklegacy.h"

#define TRUE  1
#define FALSE 0

void mexFunction( int nlhs, mxArray *plhs[],    
                  int nrhs, const mxArray *prhs[] )
{
  double *f;
  
  int* isfree;
  double* x;
  double* g;
  double* h;
  int *task;
  opl_vmlmb_workspace_t* ws;
  int n, m;
  
  if (nrhs < 4 || nrhs > 6) {
    mexErrMsgTxt("expecting between 4 and 6 arguments");
  }
  if (nlhs != 1) {
        mexErrMsgTxt("1 output argument allowed."); 
  } 
  
  ws =  (opl_vmlmb_workspace_t*)mxGetPr(prhs[0]); // get the workspace
  n = opl_vmlmb_get_n(ws);   // number of variable 

    /* Control the input x */
    if ( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]))
        mexErrMsgTxt("The second input x must be a real valued vector.");
    if (mxGetM(prhs[1])*mxGetN(prhs[1])!=n)
        mexErrMsgTxt("Incorrect dimension for the second input vector x");
    x = mxGetPr(prhs[1]);

    
    /* Control the input f */
    if (!mxIsDouble(prhs[2]) || !mxIsScalar(prhs[2])  || mxIsComplex(prhs[2]) )
        mexErrMsgTxt("The third input f must be a real valued scalar.");   
     f = mxGetPr(prhs[2]);
    

    /* Control the input g */
    if (!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) )
        mexErrMsgTxt("The fourth input g must be a real valued vector");
      if (mxGetM(prhs[3])*mxGetN(prhs[3])!=n)
        mexErrMsgTxt("Incorrect dimension for the fourth input vector g");
    g = mxGetPr(prhs[3]);

    if(nrhs>4){
    /* Control the input isfree */
    if ( mxGetM(prhs[4])*mxGetN(prhs[4])==0 )
        isfree = NULL;
    else {
        if ( !mxIsNumeric(prhs[4]) || (mxGetClassID(prhs[4])!=mxINT32_CLASS) || mxIsComplex(prhs[4]) )
            mexErrMsgTxt("The fifth input isfree must be an integer vector or an empty vector");
         if (mxGetM(prhs[4])*mxGetN(prhs[4])!=n)
        mexErrMsgTxt("Incorrect dimension for the fourth input vector active");
   isfree = (opl_integer_t *)mxGetPr(prhs[4]);
    }
    }else{
        isfree=NULL;
        }
    if(nrhs>5){
                if  (!mxIsDouble(prhs[5]) || mxIsComplex(prhs[5]) )
            mexErrMsgTxt("The sixth input h must be a real valued vector");
      if (mxGetM(prhs[5])*mxGetN(prhs[5])!=n)
        mexErrMsgTxt("Incorrect dimension for the sixth input vector g");
    h = mxGetPr(prhs[5]);
        }else{
            h =NULL;
            }
    ;
    mwSize dims[]={1,1};
    plhs[0] = mxCreateNumericArray(2,dims,mxINT16_CLASS,mxREAL);
   task = (opl_integer_t *)mxGetPr(plhs[0]);  /* opl_integer_t are 16 bits integers */
       //mexPrintf("task: %d \n", *task);

    task[0] = opl_vmlmb_iterate( ws, x, f, g, isfree, h);
    
           }