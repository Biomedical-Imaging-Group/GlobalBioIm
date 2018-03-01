/* 
 *  m_opl_vmlmb_create.c
 *
 * function m_opl_vmlmb_create
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

#include "mex.h"
#include "optimpacklegacy.h"

#define TRUE  1
#define FALSE 0

void mexFunction( int nlhs, mxArray *plhs[],    
                  int nrhs, const mxArray *prhs[] )
{ 
    double fatol, frtol, sftol, sgtol, sxtol, epsilon, delta;
    int n, m,size, fmin=FALSE;
    mwSize dims[]={1,1};
    opl_vmlmb_workspace_t* ws;
      /* Control of the number of inputs and outputs */ 
    if (nrhs > 9)
        mexErrMsgTxt("9 input argument required."); 
    else if (nlhs !=1) 
        mexErrMsgTxt("1 output argument required.");
    
  
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
    delta = mxGetScalar(prhs[8]);
    
    size = opl_vmlmb_monolithic_workspace_size(n, m);
    dims[1] = size;
    plhs[0] = mxCreateNumericArray(2,dims,mxINT8_CLASS,mxREAL);
    
    ws = opl_vmlmb_monolithic_workspace_init((char *)mxGetPr(plhs[0]), n, m);
    
    if (ws == NULL) {
    if (errno == ENOMEM) {
      mexErrMsgTxt("insufficient memory");
    } else {
      mexErrMsgTxt("unknown error during workspace creation");
    }
  }
  /* Configure VMLMB instance  */
# define SET_ATTRIBUTE(name, invalid)                           \
    {double value = name ;                                      \
    if ((invalid) ||                                            \
        opl_vmlmb_set_##name(ws, value) != OPL_SUCCESS) {       \
      mexErrMsgTxt("invalid value for `" #name "`");            \
  }}
  SET_ATTRIBUTE(fmin, FALSE);
  SET_ATTRIBUTE(fatol, value < 0);
  SET_ATTRIBUTE(frtol, value < 0);
  SET_ATTRIBUTE(sftol, value <= 0 || value >= 1);
  SET_ATTRIBUTE(sgtol, value <= 0 || value >= 1);
  SET_ATTRIBUTE(sxtol, value <= 0 || value >= 1);
  SET_ATTRIBUTE(delta, value < 0);
  SET_ATTRIBUTE(epsilon, value < 0);
# undef SET_ATTRIBUTE

}

