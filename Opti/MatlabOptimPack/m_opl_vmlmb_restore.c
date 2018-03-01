/* 
 *  m_opl_vmlmb_restore.c
 *
 * function m_opl_vmlmb_restore
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
 */

#include <errno.h>
#include<string.h>

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
    void* xo;
  double* g;
  opl_task_t  task;
  opl_vmlmb_workspace_t* ws;
  int n, m;
  
  if (nrhs != 4) {
    mexErrMsgTxt("expecting between 4 arguments");
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

    mwSize *dims = mxGetDimensions(prhs[1]);
    
    plhs[0] = mxCreateNumericArray(mxGetNumberOfDimensions(prhs[1]), dims,mxGetClassID(prhs[1]),mxREAL);
   xo = mxGetPr(plhs[0]);

    task = opl_vmlmb_restore( ws, x, f, g);
    memcpy(xo, x, n*sizeof(x[0]));    
}