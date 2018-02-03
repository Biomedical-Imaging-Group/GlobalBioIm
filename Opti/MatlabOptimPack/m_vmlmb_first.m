% function m_vmlmb_first 
%
% [csave, isave, dsave] = m_vmlmb_first(n, m, fatol, frtol, sftol, ...
%                                       sgtol, sxtol, epsilon, costheta);
%
% Matlab interface to the op_vmlmb_first function of the OptimPack package
% See the OptimPack doc below for more information
%
%  September 2013
%  IRAP, Observatoire Midi-Pyrénées, Toulouse
%  Hervé Carfantan
%
%
% /* VMLM-B computes  a local  minimizer of  a function of  N variables  by a
%  * limited memory variable metric (BFGS) method; optionally, the parameters
%  * may be bounded.  The user must evaluate the function and the gradient.
%  *
%  * VMLM-B   is   implemented   via   two  functions:   op_vmlmb_first   for
%  * initialization   and  op_vmlmb_next   for  further   iterations.   These
%  * functions use  reverse communication.  The  user must choose  an initial
%  * approximation X to the minimizer, evaluate the function and the gradient
%  * at X, and make the initial call  with TASK set to "start".  On exit TASK
%  * indicates the required action.
%  *
%  * The arguments are:
%  *
%  *   N is the number of parameters.
%  *
%  *   M is  the number of correction  pairs to remember in  order to compute
%  *       the  limited memory  variable metric  (BFGS) approximation  of the
%  *       inverse of the Hessian.  For large problems, M = 3 to 5 gives good
%  *       results.  For  small problems, M should  be less or  equal N.  The
%  *       larger is  M (and N)  the more computer  memory will be  needed to
%  *       store the workspaces (see DSAVE).
%  *
%  *   FRTOL   is  the  relative   error  desired   in  the   function  (e.g.
%  *       FRTOL=1e-8).  Convergence  occurs if the estimate  of the relative
%  *       error between F(X)  and F(XSOL), where XSOL is  a local minimizer,
%  *       is less or  equal FRTOL.  FRTOL must have  a non-negative floating
%  *       point value.
%  *
%  *   FATOL is the absolute error  desired in the function (e.g. FATOL=0.0).
%  *       Convergence occurs  if the estimate of the  absolute error between
%  *       F(X)  and F(XSOL), where  XSOL is  a local  minimizer, is  less or
%  *       equal FATOL.  FATOL must have a non-negative floating point value.
%  *
%  *   SFTOL, SGTOL, and SXTOL are  tolerances for the line search subroutine
%  *       (see  op_csrch).    Recommended  values:  SFTOL=0.001,  SGTOL=0.9,
%  *       SXTOL=0.1  (other   values  may   be  more  suitable   for  highly
%  *       non-quadratic penalty function).
%  *
%  *   EPSILON is a small (strictly positive) value used to discard BFGS updates
%  *       that may yield a non positive definite Hessian approximation.
%  *
%  *   COSTHETA is a small value, in the range [0,1), equals to the cosine of
%  *       the maximum angle between the search direction and the anti-gradient.
%  *       The BFGS recursion is restarted, whenever the search direction is not
%  *       sufficiently "descending".
%  *
%  *   CSAVE is  a character workspace array  of length OP_VMLMB_CSAVE_NUMBER
%  *       (same  as   OP_MSG_SIZE)  which   is  used  to   store  additional
%  *       information on exit with convergence, a warning or an error.
%  *
%  *   ISAVE is an integer workspace array of length OP_VMLMB_ISAVE_NUMBER.
%  *
%  *   DSAVE is a floating point workspace array of length equal to the value
%  *       returned by the macro OP_VMLMB_DSAVE_NUMBER(N, M):
%  *           26 + N + 2*M*(N + 1).
%  *
%  *   X is a  double  precision  array of  length  N.   On  entry,  X is  an
%  *       approximation to the solution.   On exit with TASK=OP_TASK_CONV, X
%  *       is the current approximation.
%  *
%  *   F is the  address of a double precision variable.  On  entry, F is the
%  *       value of  the function  at X.   On final exit,  F is  the function
%  *       value at X.
%  *
%  *   G is a double  precision array of length N.  On entry,  G is the value
%  *       of  the gradient  at X.   On final  exit, G  is the  value  of the
%  *       gradient at X.
%  *
%  *   ACTIVE  is an optional  integer array  with length  N provided  by the
%  *       caller if the  values in X has bounds.  If  the parameters have no
%  *       bounds,  ACTIVE  should   be  NULL  (unconstrained  minimization).
%  *       Otherwise,  elements  set to  zero  in  ACTIVE  indicate that  the
%  *       corresponding values  in X has reached  a bound and  should not be
%  *       changed during  the next step  because the gradient has  the wrong
%  *       sign (i.e.  the steepest descent direction would violate the bound
%  *       constraints):
%  *           ACTIVE[i] = 0 if i-th value has a lower bound XLO[i]
%  *                           and X[i]=XLO[i] and G[i]>=0 
%  *                       0 if i-th value has an upper bound XHI[i]
%  *                           and X[i]=XHI[i] and G[i]<=0
%  *                       1 (or any non-zero value) otherwise
%  *
%  *       ACTIVE needs  only to be  computed (and specified) the  first time
%  *       op_vmlmb_next is called and  when TASK=OP_TASK_NEWX (i.e.  after a
%  *       successful   step).    ACTIVE   may   also   be   specified   when
%  *       TASK=OP_TASK_CONV  (i.e.   after  convergence  if caller  wish  to
%  *       continue with  minimization).  If X has (some)  bounds, the caller
%  *       is responsible for applying the  bounds to X before evaluating the
%  *       function value  F and the gradient G  (i.e. when TASK=OP_TASK_FG),
%  *       e.g.:
%  *           if (X[i] < XLO[i]) X[i] = XLO[i];
%  *           if (X[i] > XHI[i]) X[i] = XHI[i];
%  *
%  *       If H is  not specified (i.e. H is  NULL) or if H[i] > 0  for all i
%  *       such that ACTIVE[i] is non-zero, then ACTIVE is left unchanged.
%  *
%  *   H is an optional double precision  array with length N provided by the
%  *       caller and such that diag(H) is an approximation of the inverse of
%  *       the Hessian matrix.  If H is NULL, then the inverse of the Hessian
%  *       is approximated by a simple rescaling using Shanno & Phua formula.
%  *       Otherwise, if ACTIVE  is NULL, all elements of  H must be strictly
%  *       greater than  zero; else  ACTIVE[i] is  set to zero  if H[i]  <= 0
%  *       (this is the only case  where ACTIVE is modified).  As for ACTIVE,
%  *       H needs only to be specifed  the first time op_vmlmb is called and
%  *       when JOB=2.
%  *
%  *   TASK is  the value returned  by op_vmlmb_first and  op_vmlmb_next.  It
%  *       can have one of the following values:
%  *           OP_TASK_FG - caller must evaluate the function and gradient at
%  *               X and call op_vmlm_next.
%  *           OP_TASK_NEWX  -   a  new  iterate  has   been  computed.   The
%  *               approximation X, function F,  and gradient G are available
%  *               for examination.
%  *
%  *           OP_TASK_CONV  -  the  search  is  successful.   The  solution,
%  *               function value and gradient are available in X, F and G.
%  *
%  *           OP_TASK_WARN -  VMLMB is not  able to satisfy  the convergence
%  *               conditions.   The  exit  value  of  X  contains  the  best
%  *               approximation found so  far.  Warning message is available
%  *               in CSAVE.
%  *           OP_TASK_ERROR then  there is an error in  the input arguments.
%  *               Error message is available in CSAVE.
%  *
%  * The caller must  not modify the workspace arrays  CSAVE, ISAVE and DSAVE
%  * between calls to op_vmlmb_first and further calls to op_vmlmb_next.
%  *
%  * A  typical invocation of  VMLMB for  unconstrained minimization  has the
%  * following outline:
%  *
%  *    // Choose a starting vector:
%  *    for (i=0 ; i<n ; ++i) x[i] = ... ;
%  *  
%  *    // Allocate and setup workspaces:
%  *    csave = malloc(OP_VMLMB_CSAVE_NUMBER*sizeof(char));
%  *    isave = malloc(OP_VMLMB_ISAVE_NUMBER*sizeof(op_integer_t));
%  *    dsave = malloc(OP_VMLMB_DSAVE_NUMBER(n, m)*sizeof(double));
%  *    task = op_vmlmb_first(n, m, fmin, fatol, frtol, sftol, sgtol, sxtol,
%  *                          csave, isave, dsave);
%  *    for (;;) {
%  *      if (task == OP_TASK_FG) {
%  *        f = ...;  // evaluate the function at X; store in F
%  *        g = ...;  // evaluate the gradient of F at X; store in G
%  *      } else if (task == OP_TASK_NEWX) {
%  *         // New successful step: the approximation X, function F, and
%  *         // gradient G, are available for inspection.
%  *      } else {
%  *        // Convergence, or error, or warning
%  *        fprintf(stderr, "%s\n", csave);
%  *        break;
%  *      }
%  *      // Computes next step:
%  *      task = op_vmlmb_next(x, &f, g, NULL, NULL, csave, isave, dsave);
%  *    }
%  *
%  * A typical invocation of VMLMB for bound-constrained minimization has the
%  * following outline:
%  *
%  *    // Choose a starting vector:
%  *    for (i=0 ; i<n ; ++i) x[i] = ... ;
%  *  
%  *    // Allocate and setup workspaces:
%  *    csave = malloc(OP_VMLMB_CSAVE_NUMBER*sizeof(char));
%  *    isave = malloc(OP_VMLMB_ISAVE_NUMBER*sizeof(op_integer_t));
%  *    dsave = malloc(OP_VMLMB_DSAVE_NUMBER(n, m)*sizeof(double));
%  *    task = op_vmlmb_first(n, m, fmin, fatol, frtol, sftol, sgtol, sxtol,
%  *                          csave, isave, dsave);
%  *    eval = 0; // number of evaluations
%  *    for (;;) {
%  *      if (task == OP_TASK_FG) {
%  *        op_bounds_apply(n, x, xmin, xmax); // aply bound constraints
%  *        f = ...;  // evaluate the function at X; store in F
%  *        g = ...;  // evaluate the gradient of F at X; store in G
%  *        ++eval;
%  *      } else if (task == OP_TASK_NEWX) {
%  *         // New successful step: the approximation X, function F, and
%  *         // gradient G, are available for inspection.
%  *      } else {
%  *        // Convergence, or error, or warning
%  *        fprintf(stderr, "%s\n", csave);
%  *        break;
%  *      }
%  *      // Computes next step:
%  *      if (eval == 1 || task == OP_TASK_NEWX) {
%  *        // Computes set of active parameters:
%  *        op_bounds_active(n, active, x, g, xmin, xmax);
%  *      }
%  *      task = op_vmlmb_next(x, &f, g, active, NULL, csave, isave, dsave);
%  *    }
%  *
%  *
%  * HISTORY:
%  *   MINPACK-2 Project. April 1995.
%  *   Argonne National Laboratory and University of Minnesota.
%  *   Brett M. Averick, Richard G. Carter, and Jorge J. Moré.
%  *
%  *   C-version and improvements (bound constraints, preconditioning, ...).
%  *   February 2003 - March 2003.
%  *   Observatoire de Lyon.
%  *   Eric Thiebaut.
%  */