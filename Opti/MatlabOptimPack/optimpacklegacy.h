/*
 * optimpacklegacy.h --
 *
 * Definitions for optimization routines implemented in Optimpacklegacy
 * library.
 *
 *-----------------------------------------------------------------------------
 *
 * Copyright (c) 2003, 2016 Éric Thiébaut.
 *
 * This file is part of OptimPack <https://github.com/emmt/OptimPackLegacy>.
 *
 * OptimPack is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 2 of the License, or (at your option)
 * any later version.
 *
 * OptimPack is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * OptimPack (file "LICENSE" in the top source directory); if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307 USA
 *
 *-----------------------------------------------------------------------------
 */

#ifndef _OPTIMPACKLEGACY_H
#define _OPTIMPACKLEGACY_H 1

/**
 * Customizable data types:
 *   OPL_INTEGER = data type used to store array indices
 *   OPL_LOGICAL = data type of the result of a logical test
 */
#ifndef OPL_INTEGER
# define OPL_INTEGER int
#endif
#ifndef OPL_LOGICAL
# define OPL_LOGICAL int
#endif

typedef OPL_INTEGER opl_integer_t;
typedef OPL_LOGICAL opl_logical_t;


/*---------------------------------------------------------------------------*/
/* USEFUL MACROS */

/* C++ needs to know that types and declarations are C, not C++. */
#ifdef  __cplusplus
# define _OPL_BEGIN_DECLS  extern "C" {
# define _OPL_END_DECLS    }
#else
# define _OPL_BEGIN_DECLS  /* empty */
# define _OPL_END_DECLS    /* empty */
#endif
_OPL_BEGIN_DECLS

/**
 * Possible values for a boolean.
 */
typedef enum {
  OPL_FALSE = 0,
  OPL_TRUE = 1
} opl_boolean_t;


/*---------------------------------------------------------------------------*/
/* STATUS AND ERROR REPORTING */

/**
 * Standard result codes
 *
 * @brief These constants are the standard value returned by most dynamic
 *        buffer routines.
 */
typedef enum {
  /* Success. */
  OPL_SUCCESS = 0,

  /* Warnings. */
  OPL_STP_EQ_STPMIN,
  OPL_STP_EQ_STPMAX,
  OPL_XTOL_TEST_SATISFIED,
  OPL_ROUNDING_ERROR,

  /* Errors. */
  OPL_STPMAX_LT_STPMIN,
  OPL_STPMIN_LT_ZERO,
  OPL_XTOL_LT_ZERO,
  OPL_FTOL_LE_ZERO,
  OPL_GTOL_LE_ZERO,
  OPL_NOT_A_DESCENT,
  OPL_STP_GT_STPMAX, /* OPL_OUT_OF_BOUNDS? */
  OPL_STP_LT_STPMIN, /* OPL_OUT_OF_BOUNDS? */
  OPL_F_LE_FMIN,
  OPL_NOT_POSITIVE_DEFINITE,
  OPL_INSUFFICIENT_MEMORY,
  OPL_ILLEGAL_ADDRESS,
  OPL_INVALID_ARGUMENT,
  OPL_OUT_OF_BOUNDS,
  OPL_CORRUPTED,
  OPL_OVERFLOW,
#if 0
  OPL_SYNTAX_ERROR,
  OPL_OPEN_ERROR,
  OPL_CLOSE_ERROR,
  OPL_READ_ERROR,
  OPL_WRITE_ERROR,
  OPL_STAT_ERROR,
  OPL_FILE_EXISTS,
  OPL_ZIP_ERROR,
#endif
  OPL_SYSTEM_ERROR   /* System error, use errno code. */
} opl_status_t;

/**
 * The context opaque structure.
 */
typedef struct _opl_context opl_context_t;

/**
 * Query status from context.
 *
 * @param ctx   Context for error reporting.
 *
 * @return A status value.
 */
extern opl_status_t opl_get_status(opl_context_t* ctx);

/**
 * Query system error code from context.
 *
 * @param ctx   Context for error reporting.
 *
 * @return If current context status is `OPL_SYSTEM_ERROR`, the system error
 *         code (a.k.a. `errno`) when the error occured is returned; otherwise,
 *         0 is returned.
 */
extern int opl_get_errno(opl_context_t* ctx);

/**
 * Query message from context.
 *
 * @param ctx   Context for error reporting.
 *
 * @return A human readable message.
 */
extern const char* opl_get_message(opl_context_t* ctx);

/**
 * Query an explanation about a status code.
 *
 * @param status   Status value.
 *
 * @return A human readable message.
 */
extern const char* opl_get_default_message(opl_status_t status);

/**
 * Intialize a context.
 *
 * This function initializes a context to some consistent contents.  It must
 * be called before using the context.
 *
 * @param ctx     The address of the context.
 */
extern void opl_initialize_context(opl_context_t* ctx);

/**
 * The different strorage types.
 */
typedef enum {
  OPL_VOLATILE = 0, /** Storage is volatile, a copy must be done. */
  OPL_PERMANENT     /** Storage is permanent, the same address can be used. */
} opl_storage_type_t;

/**
 * Set a context contents.
 *
 * This function set the contents, that is the status and corresponding
 * message, of a given context and return the status value.
 *
 * The given message must be a null terminated string.  If the storage type of
 * the message is `OPL_PERMANENT`, the context message will point to `message`;
 * otherwise, the context message will point to the internal buffer of `ctx`
 * after copying `message` in this buffer.  In this latter case, if the message
 * is too long to fit in the buffer, it will be truncated and terminated by the
 * ellipsis text "[...]".
 *
 * @param ctx       The address of the context.
 * @param status    The value of the status.
 * @param message   The message string.
 * @param type      The storage type of the message string.
 *
 * @return The value of `status`.
 */
extern opl_status_t opl_set_context(opl_context_t* ctx,
                                    opl_status_t status,
                                    const char* message,
                                    opl_storage_type_t type);

/**
 * Set a context with a formated message.
 *
 * This function set the contents, that is the status and corresponding
 * message, of a given context and return the status value.  The message is
 * formatted using the same mechanism as `printf`.  If the resulting message is
 * too long to fit in the internal buffer of the context, it will be truncated
 * and terminated by the ellipsis text "[...]".
 *
 * @param ctx      The address of the context.
 * @param status   The value of the status.
 * @param format   The format string.
 * @param ...      The arguments of the format.
 *
 * @return The value of `status`.
 */
extern opl_status_t opl_format_context(opl_context_t* ctx,
                                       opl_status_t status,
                                       const char* format, ...);

/**
 * Report successful operation.
 *
 * This (inlined) function is equivalent to:
 * @code
 *     opl_set_context(ctx, OPL_SUCCESS,
 *                     opl_get_default_message(OPL_SUCCESS),
 *                     OPL_PERMANENT)
 * @endcode
 *
 * @param ctx   Context for error reporting.
 *
 * @return `OPL_SUCCESS`.
 */
extern opl_status_t opl_success(opl_context_t* ctx);

/*---------------------------------------------------------------------------*/
/* LINE SEARCH */

/**
 * @addtogroup MoreThuenteLineSearch
 * @{
 *
 * @page csrch Moré & Thuente line search
 * @tableofcontents
 *
 * Moré & Thuente line search method (@ref morethuente1984) finds a step that
 * satisfies a sufficient decrease condition and a curvature condition.
 *
 * The function @link #opl_csrch_start should be called at the start of a line
 * search to initiate the search.  The function @link #opl_csrch_iterate should
 * then be called for each new tried point along the line search until the
 * convergence of the search.
 *
 * The method updates an interval with endpoints `stx` and `sty`.  The
 * interval is initially chosen so that it contains a minimizer of the
 * modified function:
 *
 * <pre>
 *       psi(stp) = f(stp) - f(0) - ftol*stp*g(0)
 * </pre>
 *
 * where `f(0)` and `g(0) = f'(0)` are the value of the function and its
 * derivative for `stp = 0`.  If `psi(stp) ≤ 0` and `g(stp) ≥ 0` for some step,
 * then the interval is chosen so that it contains a minimizer of `f(stp)`.
 * The algorithm is designed to find a step that satisfies the sufficient
 * decrease condition:
 *
 * <pre>
 *     f(stp) <= f(0) + ftol*stp*g(0),                            (1)
 * </pre>
 *
 * and the curvature condition:
 *
 * <pre>
 *     abs(g(stp)) <= gtol*abs(g(0)).                             (2)
 * </pre>
 *
 * Relations (1) and (2) are called the strong Wolfe conditions.  If `ftol` is
 * less than `gtol` and if, for example, the function is bounded below, then
 * there is always a step which satisfies both conditions.  If no step can be
 * found that satisfies both conditions, then the algorithm stops with a
 * warning.  In this case `stp` only satisfies the sufficient decrease
 * condition.
 *
 *
 * @subsection example Example
 *
 * A typical invocation of the method has the following outline:
 *
 * <pre>
 * opl_csrch_create_workspace();
 * double f = ...;   // function value for STP=0
 * double g = ...;   // derivative value for STP=0
 * double stp = ...; // guess for next STP value to try (STP > 0.0)
 * double ftol = 1e-3;
 * double gtol = 0.9;
 * double xtol = 0.1;
 * double stpmin = 0;
 * double stpmax = 1e20*stp;
 * opl_csrch_start(ws, f, g, stp, ftol, gtol, xtol, stpmin, stpmax);
 * for (;;) {
 *     opl_task_t task = opl_csrch_get_task(ws);
 *     if (task == OPL_TASK_FG) {
 *         // Evaluate the function and the gradient at `stp`.
 *         f = func(stp);
 *         g = grad(stp);
 *     } else if (task == OPL_TASK_CONV) {
 *         // Search has converged.
 *         break;
 *     } else if (task == OPL_TASK_WARN) {
 *         // Some problem prevents further progress.
 *         fprintf(stderr, "warning in %s\n", opl_csrch_get_reason(ws));
 *         exit(1);
 *     } else {
 *         // An error occured.
 *         fprintf(stderr, "error in %s\n", , opl_csrch_get_reason(ws));
 *         exit(1);
 *     }
 * }
 * </pre>
 *
 * @subsection refs References
 *
 * @anchor morethuente1984
 *     Jorge J. Moré and David J. Thuente, "Line search algorithms with
 *     guaranteed sufficient decrease" in ACM Transactions on Mathematical
 *     Software (TOMS) Volume 20, Issue 3, Pages 286-307 (September 1994).
 *
 *
 * @subsection hist History
 *
 * MINPACK-1 Project. June 1983.
 * Argonne National Laboratory.
 * Jorge J. Moré and David J. Thuente.
 *
 * MINPACK-2 Project. November 1993.
 * Argonne National Laboratory and University of Minnesota.
 * Brett M. Averick, Richard G. Carter, and Jorge J. Moré.
 *
 * Yorick translation an improvements.  October 2001.
 * C-version.  February 2003.
 * Observatoire de Lyon (France).
 * Éric Thiébaut.
 *
 * New C version with structured workspace.
 * December 2016.
 * Centre de Recherche Astrophysique de Lyon (France).
 * Éric Thiébaut.
 */

/**
 * Possible values for an optimization task.
 */
typedef enum {
  OPL_TASK_START     = 0, /*< start line search */
  OPL_TASK_FG        = 1, /*< caller has to compute function and gradient */
  OPL_TASK_FREEVARS  = 2, /*< caller has to determine the free variables */
  OPL_TASK_NEWX      = 3, /*< new variables available for inspection */
  OPL_TASK_CONV      = 4, /*< search has converged */
  OPL_TASK_WARN      = 5, /*< search aborted with warning */
  OPL_TASK_ERROR     = 6  /*< search aborted with error */
} opl_task_t;

/**
 * Opaque workspace used for Moré & Thuente line search.
 */
typedef struct _opl_csrch_workspace opl_csrch_workspace_t;

/**
 * Get the pending task of a line search instance.
 */
extern opl_task_t opl_csrch_get_task(opl_csrch_workspace_t* ws);

/**
 * Get the status of the last operation on a line search instance.
 */
extern opl_task_t opl_csrch_get_status(opl_csrch_workspace_t* ws);


/**
 * Get the message explaining the result of the last operation on a line search
 * instance.
 */
extern const char* opl_csrch_get_reason(opl_csrch_workspace_t* ws);

/**
 * Get the size of a Moré & Thuente line search instance.
 *
 * This function can be used to determine the size (in bytes) of a Moré &
 * Thuente line search instance for applications which want to directly manage
 * memory.  A line search instance can be stored into any memory block of
 * sufficient size and aligned on a multiple of the size of a `double`.
 * However, the line search workspace can be assumed to have been properly
 * initilaized only after a first call to @link #opl_csrch_start.
 *
 * @return The size in bytes of a Moré & Thuente line search instance.
 */
extern size_t opl_csrch_get_workspace_size();

/**
 * Create a new instance of Moré & Thuente line search.
 *
 * @return A new instance of Moré & Thuente line search.  When no longer
 *         needed, this instance must be destroyed with @link
 *         #olp_csrch_destroy_workspace. `NULL` is returned in case of error
 *         (i.e. insufficient memory).
 */
extern opl_csrch_workspace_t* opl_csrch_create_workspace();

/**
 * Destroy an instance of Moré & Thuente line search.
 *
 * This function should only be called on an instance created by @link
 * #opl_csrch_create_workspace.
 *
 * @param ws    A line search instance created by @link
 *              #opl_csrch_create_workspace.
 */
extern void opl_csrch_destroy_workspace(opl_csrch_workspace_t* ws);

/**
 * Initiate a line search by Moré & Thuente method.
 *
 * @param ws     The address of the line search workspace.
 *
 * @param f      The value of the function at the start of the line search
 *               (that is for `stp = 0`).
 *
 * @param g      The value of the directional derivative at the start of the
 *               line search (that is for `stp = 0`).
 *
 * @param stp    The length of the first step to take along the search
 *               direction.  Must be a strictly positive value.
 *
 * @param ftol   A nonnegative tolerance for the sufficient decrease condition.
 *               One should take `0 < ftol < 1/2`.
 *
 * @param gtol   A nonnegative tolerance for the curvature condition.  One
 *               should take `ftol < gtol < 1`.
 *
 * @param xtol   A nonnegative relative tolerance for an acceptable step.
 *               The line search will exit with a warning if the relative
 *               difference between `sty` and `stx` is less than `xtol`.
 *
 * @param stpmin A nonnegative lower bound for the step.
 *
 * @param stpmax A nonnegative upper bound for the step.
 *
 * @return A standard status: `OPL_SUCCESS` or any other value to indicate the
 *         error.
 */
extern opl_status_t
opl_csrch_start(opl_csrch_workspace_t* ws, double f, double g, double stp,
                double ftol, double gtol, double xtol,
                double stpmin, double stpmax);

/**
 * Iterate a line search by Moré & Thuente method.
 *
 * Upon return, @link #opl_csrch_get_task, @link #opl_csrch_get_status or @link
 * #opl_csrch_get_reason can be used to retrieve the next pending task, the
 * status of the last operation and the reason of the result of the las
 * operation.
 *
 * @param ws     The address of the line search workspace.
 *
 * @param f      The value of the function for the current step along the line
 *               search.
 *
 * @param g      The value of the derivative of the function for the current
 *               step along the line search.
 *
 * @param stp    The address of the variable with, on entry, the length of the
 *               current step taken along the search direction.  On exit with
 *               task `OPL_TASK_FG`, the variable will be set with the next
 *               step to take.  On exit with task `OPL_TASK_CONV`, the variable
 *               is left unchanged.
 *
 * @return The following status values can be returned:
 *
 * - `OPL_SUCCESS`, then the task returned by @link #opl_csrch_get_task
 *   indicates the required action.  If `task = OPL_TASK_FG`, then the callar
 *   shall evaluate the function and derivative at `*stp` and call
 *   `opl_csrch_iterate` again.  If `task = OPL_TASK_CONV`, then the search is
 *   successful.  If `task = OPL_TASK_WARN` then the subroutine is not able to
 *   satisfy the convergence conditions. The exit value of `*stp` contains the
 *   best point found during the search.  If `task = OPL_TASK_ERROR` then there
 *   is an error in the input arguments.
 *
 * - `OPL_STP_EQ_STPMIN` with `task = OPL_TASK_WARN`: the search is blocked at
 *   the lower bound.
 *
 * - `OPL_STP_EQ_STPMAX` with `task = OPL_TASK_WARN`: the search is blocked at
 *   the upper bound.
 *
 * - `OPL_XTOL_TEST_SATISFIED` with `task = OPL_TASK_WARN`: `XTOL` test is
 *   satisfied.
 *
 * - `OPL_TASK_WARN` with `task = OPL_TASK_WARN`: rounding errors prevent
 *   progress.
 *
 * - `OPL_OUT_OF_BOUNDS` with `task = OPL_TASK_ERROR`: the step `*stp` is
 *   outside bracket `(stx,sty)`.
 *
 * - `OPL_NOT_A_DESCENT` with `task = OPL_TASK_ERROR`: the descent condition
 *   does not hold.
 *
 * - `OPL_STPMAX_LT_STPMIN` with `task = OPL_TASK_ERROR`: the search bounds are
 *   incompatible.
 */
extern opl_status_t
opl_csrch_iterate(opl_csrch_workspace_t* ws, double f, double g,
                  double *stp);

/**
 * Compute a safeguarded cubic step.
 *
 * This function computes a safeguarded step for a search procedure and updates
 * an interval that contains a step that satisfies a sufficient decrease and a
 * curvature condition [1].
 *
 * The parameter `stx` contains the step with the least function value.  If
 * `brackt` is set to true (i.e.  non-zero) then a minimizer has been bracketed
 * in an interval with endpoints `stx` and `sty`.  The parameter `stp` contains
 * the current step.  The subroutine assumes that if `brackt` is true then:
 *
 * <pre>
 *     min(stx,sty) < stp < max(stx,sty),
 * </pre>
 *
 * and that the derivative at `stx` is negative in the direction of the step.
 *
 * On output, all the parameters passed by address are updated appropriately.
 *
 * @param ctx      The address of a context structure for error reporting.
 *
 * @param brackt_ptr The addresses where the value of `brackt` is stored.
 *                  `brackt` is a logical variable.  On entry, `brackt`
 *                  specifies if a minimizer has been bracketed.  Initially
 *                  `brackt` must be set to false.  On exit, `brackt` specifies
 *                  if a minimizer has been bracketed.  When a minimizer is
 *                  bracketed, `brackt` (i.e. the value at address
 *                  `brackt_ptr`) is set to true.
 *
 * @param stpmin    The lower bound for the step.
 *
 * @param stpmax    The upper bound for the step.
 *
 * @param stx_ptr   The address of `stx`, the best step obtained so far.
 *
 * @param fx_ptr    The address of `fx`, the function at `stx`.
 *
 * @param dx_ptr    The addresses of `dx`, the derivative at `stx` negative in
 *                  the direction of the step, that is, `dx` and `stp - stx`
 *                  must have opposite signs.
 *
 * @param sty_ptr   The address of `sty`, the other endpoint of the interval
 *                  of uncertainty.
 *
 * @param fy_ptr    The address of `fy`, the function at `sty`.
 *
 * @param dy_ptr    The address of `dy`, the derivative at `sty`.
 *
 * @param stp_ptr   The address where the value of `stp`, the current step,
 *                  is stored.  If `brackt` is set true then on input `stp`
 *                  must be between `stx` and `sty`.  On output, the value at
 *                  `stp_ptr` is set to the new step.
 *
 * @param fp        The function at the current step.
 *
 * @param dp        The derivative at the current step.
 *
 * @return The returned value indicates whether the operation was successful
 *         (failure can only occur with wrong arguments).  The context `ctx` is
 *         used to report errors if any.
 */
extern opl_status_t
opl_cstep(opl_context_t* ctx, opl_boolean_t *brackt,
          double stpmin, double stpmax,
          double *stx_ptr, double *fx_ptr, double *dx_ptr,
          double *sty_ptr, double *fy_ptr, double *dy_ptr,
          double *stp_ptr, double  fp,     double  dp);
/**
 * @}
 */

/*---------------------------------------------------------------------------*/

/**
 * @addtogroup VMLMB
 * @{
 *
 * VMLMB is a limited memory variable metric method (BFGS) which can take into
 * account bound constraints
 */

/**
 * The opaque workspace for VMLMB.
 */
typedef struct _opl_vmlmb_workspace opl_vmlmb_workspace_t;

/**
 * Return the size of a VMLMB workspace as one block.
 *
 * The returned size is a multiple of the size of a double precision floating
 * point value.  It may be stored in an array of such values.
 *
 * Typical usage:
 * <pre>
 *    size = opl_vmlmb_monolithic_workspace_size(n, m);
 *    buffer = malloc(size);
 *    ws = opl_vmlmb_monolithic_workspace_init(buffer, n, m);
 *    ....;
 *    free(buffer);
 * </pre>
 * As it is guaranteed that the base address does not change, you can also do:
 * </pre>
 *    size = opl_vmlmb_monolithic_workspace_size(n, m);
 *    ws = opl_vmlmb_monolithic_workspace_init(malloc(size), n, m);
 *    ....;
 *    free(ws);
 * <pre>
 *
 * @param n   The number of variables.
 *
 * @param m   The number of memorized steps.
 *
 * @return The size in bytes of a monolithic workspace for VMLMB.  Zero is
 *         returned if the arguments are invalid (not strictly positive).
 */
extern size_t
opl_vmlmb_monolithic_workspace_size(opl_integer_t n, opl_integer_t m);

/**
 * Initialize a workspace allocated as one block.
 *
 * @param buf   The memory allocated by the caller of the workspace, it must
 *              have at least the size given by
 *              `opl_vmlmb_monolithic_workspace_size(n, m)`.
 *
 * @param n     The number of variables.
 *
 * @param m     The number of memorized steps.
 *
 * @return The workspace initialized with defaults.  `NULL` is returned
 *         if `buf` is not a valid address.
 */
extern opl_vmlmb_workspace_t*
opl_vmlmb_monolithic_workspace_init(void* buf,
                                    opl_integer_t n, opl_integer_t m);

extern opl_vmlmb_workspace_t*
opl_vmlmb_set_defaults(opl_vmlmb_workspace_t* ws);

extern opl_vmlmb_workspace_t*
opl_vmlmb_create(opl_integer_t n, opl_integer_t m);

extern void
opl_vmlmb_destroy(opl_vmlmb_workspace_t* ws);

/*
 * VMLM-B computes a local minimizer of a function of N variables by a limited
 * memory variable metric (BFGS) method; optionally, the parameters may be
 * bounded.  The user must evaluate the function and the gradient.
 *
 * VMLM-B is implemented via two functions: opl_vmlmb_setup for initialization
 * and opl_vmlmb_iterate for further iterations.  These functions use reverse
 * communication.  The user must choose an initial approximation X to the
 * minimizer, evaluate the function and the gradient at X, and make the initial
 * call with TASK set to OPL_TASK_FG.  On exit TASK indicates the required
 * action.
 *
 * The arguments are:
 *
 *   N is the number of parameters.
 *
 *   M is the number of correction pairs to remember in order to compute the
 *       limited memory variable metric (BFGS) approximation of the inverse of
 *       the Hessian.  For large problems, M = 3 to 5 gives good results.  For
 *       small problems, M should be less or equal N.  The larger is M (and N)
 *       the more computer memory will be needed to store the workspaces (see
 *       DSAVE).
 *
 *   FRTOL is the relative error desired in the function (e.g.  FRTOL=1e-8).
 *       Convergence occurs if the estimate of the relative error between F(X)
 *       and F(XSOL), where XSOL is a local minimizer, is less or equal FRTOL.
 *       FRTOL must have a non-negative floating point value.
 *
 *   FATOL is the absolute error desired in the function (e.g. FATOL=0.0).
 *       Convergence occurs if the estimate of the absolute error between F(X)
 *       and F(XSOL), where XSOL is a local minimizer, is less or equal FATOL.
 *       FATOL must have a non-negative floating point value.
 *
 *   SFTOL, SGTOL, and SXTOL are tolerances for the line search subroutine (see
 *       opl_csrch).   Recommended  values: SFTOL=0.001,  SGTOL=0.9,  SXTOL=0.1
 *       (other values  may be more  suitable for highly  non-quadratic penalty
 *       function).
 *
 *   DELTA is a small nonegative value used to compute a small initial step.
 *
 *   EPSILON is a small value, in the range [0,1), equals to the cosine of the
 *       maximum angle between the search direction and the anti-gradient.  The
 *       BFGS recursion is restarted, whenever the search direction is not
 *       sufficiently "descending".
 *
 *   CSAVE is a character workspace array of length OPL_VMLMB_CSAVE_NUMBER
 *       (same as OPL_MSG_SIZE) which is used to store additional information
 *       on exit with convergence, a warning or an error.
 *
 *   ISAVE is an integer workspace array of length OPL_VMLMB_ISAVE_NUMBER.
 *
 *   DSAVE is a floating point workspace array of length equal to the value
 *       returned by the macro OPL_VMLMB_DSAVE_NUMBER(N, M):
 *
 *           26 + N + 2*M*(N + 1).
 *
 *   X is a double precision array of length N.  On entry, X is an
 *       approximation to the solution.  On exit with TASK = OPL_TASK_CONV, X
 *       is the current approximation.
 *
 *   F is the address of a double precision variable.  On entry, F is the value
 *       of the function at X.  On final exit, F is the function value at X.
 *
 *   G is a double precision array of length N.  On entry, G is the value of
 *       the gradient at X.  On final exit, G is the value of the gradient at
 *       X.
 *
 *   ISFREE is an optional integer array with length N provided by the caller
 *       if the values in X have bounds.  If the parameters have no bounds,
 *       ISFREE should be NULL (unconstrained minimization).  Otherwise,
 *       elements set to zero in ISFREE indicate that the corresponding values
 *       in X has reached a bound and should not be changed during the next
 *       step because the gradient has the wrong sign (i.e.  the steepest
 *       descent direction would violate the bound constraints):
 *
 *           ISFREE[i] = 0 if i-th value has a lower bound XLO[i]
 *                         and X[i] = XLO[i] and G[i] >= 0
 *                       0 if i-th value has an upper bound XHI[i]
 *                         and X[i] = XHI[i] and G[i] <= 0
 *                       1 otherwise
 *
 *       ISFREE needs only to be computed when TASK = OPL_TASK_FREEVARS.  If X
 *       has (some) bounds, the caller is responsible for applying the bounds
 *       to X before evaluating the function value F and the gradient G (i.e.,
 *       when TASK = OPL_TASK_FG), e.g.:
 *
 *           if (X[i] < XLO[i]) X[i] = XLO[i];
 *           if (X[i] > XHI[i]) X[i] = XHI[i];
 *
 *       If H is not specified (i.e., H is NULL) or if H[i] > 0 for all i such
 *       that ISFREE[i] is non-zero, then ISFREE is left unchanged.
 *
 *   H is an optional double precision array with length N provided by the
 *       caller and such that diag(H) is an approximation of the inverse of the
 *       Hessian matrix.  If H is NULL, then the inverse of the Hessian is
 *       approximated by a simple rescaling using Shanno & Phua formula.
 *       Otherwise, if ISFREE is NULL, all elements of H must be strictly
 *       greater than zero; else ISFREE[i] is set to zero if H[i] <= 0 (this is
 *       the only case where ISFREE is modified).  As for ISFREE, H needs only
 *       to be specifed the first time opl_vmlmb is called and when JOB=2.
 *
 *   TASK is the value returned by opl_vmlmb_setup and opl_vmlmb_iterate.  It
 *       can have one of the following values:
 *
 *           OPL_TASK_FG - caller must evaluate the function and gradient at X
 *               and call opl_vmlm_next.
 *
 *           OPL_TASK_FREEVARS - if variables are bounded, caller must
 *               determine the set of free variables for the current variables
 *               X and update IFREE accordingly.
 *
 *           OPL_TASK_NEWX - a new iterate has been computed.  The
 *               approximation X, function F, and gradient G are available for
 *               examination.
 *
 *           OPL_TASK_CONV - the search is successful.  The solution, function
 *               value and gradient are available in X, F and G.
 *
 *           OPL_TASK_WARN - VMLMB is not able to satisfy the convergence
 *               conditions.  The exit value of X contains the best
 *               approximation found so far.  Warning message is available in
 *               CSAVE.
 *
 *           OPL_TASK_ERROR then there is an error in the input arguments.
 *               Error message is available in CSAVE.
 *
 * The caller must not modify the workspace arrays CSAVE, ISAVE and DSAVE
 * between calls to opl_vmlmb_setup and further calls to opl_vmlmb_iterate.
 *
 * A typical invocation of VMLMB for unconstrained minimization has the
 * following outline:
 *
 *    // Choose a starting vector:
 *    for (i=0 ; i<n ; ++i) x[i] = ... ;
 *
 *    // Allocate and setup workspaces:
 *    csave = malloc(OPL_VMLMB_CSAVE_NUMBER*sizeof(char));
 *    isave = malloc(OPL_VMLMB_ISAVE_NUMBER*sizeof(opl_integer_t));
 *    dsave = malloc(OPL_VMLMB_DSAVE_NUMBER(n, m)*sizeof(double));
 *    task = opl_vmlmb_setup(n, m, fmin, fatol, frtol, sftol, sgtol, sxtol,
 *                          csave, isave, dsave);
 *    for (;;) {
 *      if (task == OPL_TASK_FG) {
 *        f = ...;  // evaluate the function at X; store in F
 *        g = ...;  // evaluate the gradient of F at X; store in G
 *      } else if (task == OPL_TASK_NEWX) {
 *         // New successful step: the approximation X, function F, and
 *         // gradient G, are available for inspection.
 *      } else {
 *        // Convergence, or error, or warning
 *        fprintf(stderr, "%s\n", csave);
 *        break;
 *      }
 *      // Computes next step:
 *      task = opl_vmlmb_iterate(x, &f, g, NULL, NULL, csave, isave, dsave);
 *    }
 *
 * A typical invocation of VMLMB for bound-constrained minimization has the
 * following outline:
 *
 *    // Choose a starting vector:
 *    for (i=0 ; i<n ; ++i) x[i] = ... ;
 *
 *    // Allocate and setup workspaces:
 *    csave = malloc(OPL_VMLMB_CSAVE_NUMBER*sizeof(char));
 *    isave = malloc(OPL_VMLMB_ISAVE_NUMBER*sizeof(opl_integer_t));
 *    dsave = malloc(OPL_VMLMB_DSAVE_NUMBER(n, m)*sizeof(double));
 *    task = opl_vmlmb_setup(n, m, fmin, fatol, frtol, sftol, sgtol, sxtol,
 *                          csave, isave, dsave);
 *    eval = 0; // number of evaluations
 *    for (;;) {
 *      if (task == OPL_TASK_FG) {
 *        opl_bounds_apply(n, x, xmin, xmax); // aply bound constraints
 *        f = ...;  // evaluate the function at X; store in F
 *        g = ...;  // evaluate the gradient of F at X; store in G
 *        ++eval;
 *      } else if (task == OPL_TASK_FREEVARS) {
 *        // Determine the set of free variables:
 *        opl_bounds_free(n, isfree, x, g, xmin, xmax);
 *      } else if (task == OPL_TASK_NEWX) {
 *         // New successful step: the approximation X, function F, and
 *         // gradient G, are available for inspection.
 *      } else {
 *        // Convergence, or error, or warning
 *        fprintf(stderr, "%s\n", csave);
 *        break;
 *      }
 *      task = opl_vmlmb_iterate(x, &f, g, isfree, NULL, csave, isave, dsave);
 *    }
 *
 *
 * HISTORY:
 *   MINPACK-2 Project. April 1995.
 *   Argonne National Laboratory and University of Minnesota.
 *   Brett M. Averick, Richard G. Carter, and Jorge J. Moré.
 *
 *   C-version and improvements (bound constraints, preconditioning, ...).
 *   February 2003 - March 2003.
 *   Observatoire de Lyon.
 *   Éric Thiébaut.
 */

extern opl_task_t
opl_vmlmb_iterate(opl_vmlmb_workspace_t* ws,
                  double x[], double *f, double g[],
                  opl_logical_t isfree[], const double h[]);

/* Set workspace data so that it can be used for a new optimization with
   the same parameters. */
extern opl_task_t
opl_vmlmb_restart(opl_vmlmb_workspace_t* ws);

/* Restore last line search starting point.  Calling this is only effective if
   task is OPL_TASK_FG. */
extern opl_task_t
opl_vmlmb_restore(opl_vmlmb_workspace_t* ws,
                  double x[], double *f, double g[]);

/**
 * Set or unset a strict lower bound for the objective function.
 *
 * VMLMB can use a strict lower bound for the objective function to estimate
 * the initial step length for the first iteration or after a restart.  If
 * `fmin` is a NaN value or a value strictly smaller than `-DBL_MAX`, it is
 * assumed that there is no such bound.  VMLMB exits with a warning if a strict
 * lower bound is set and `f(x) <= fmin`.
 *
 * @param ws      VMLMB workspace.
 *
 * @param fmin    A strict lower bound for the objective function.
 *
 * @return `OPL_SUCCESS` or any other value to indicate the reason of the
 *         failure.
 */
extern opl_status_t
opl_vmlmb_set_fmin(opl_vmlmb_workspace_t* ws, double value);

/**
 * Get the strict lower bound for the objective function.
 *
 * @param ws      VMLMB workspace.
 *
 * @return The current value of the strict lower bound set for the objective
 *         function.  A NaN value is returned if there is no such bound or in
 *         case of error (i.e. `ws` is `NULL`).
 */
extern opl_integer_t
opl_vmlmb_get_n(opl_vmlmb_workspace_t* ws);

extern opl_integer_t                                                                
opl_vmlmb_get_m(opl_vmlmb_workspace_t* ws);                                          

extern double                                
opl_vmlmb_get_fmin(opl_vmlmb_workspace_t* ws);                                         

extern opl_status_t
opl_vmlmb_set_fatol(opl_vmlmb_workspace_t* ws, double value);

extern opl_status_t
opl_vmlmb_set_frtol(opl_vmlmb_workspace_t* ws, double value);

extern opl_status_t
opl_vmlmb_set_delta(opl_vmlmb_workspace_t* ws, double value);

extern opl_status_t
opl_vmlmb_set_epsilon(opl_vmlmb_workspace_t* ws, double value);

extern opl_status_t
opl_vmlmb_set_sxtol(opl_vmlmb_workspace_t* ws, double value);

extern opl_status_t
opl_vmlmb_set_sftol(opl_vmlmb_workspace_t* ws, double value);

extern opl_status_t
opl_vmlmb_set_sgtol(opl_vmlmb_workspace_t* ws, double value);

extern opl_task_t opl_vmlmb_get_task(opl_vmlmb_workspace_t* ws);

extern opl_status_t opl_vmlmb_get_status(opl_vmlmb_workspace_t* ws);

extern const char* opl_vmlmb_get_reason(opl_vmlmb_workspace_t* ws);

extern double opl_vmlmb_get_sftol(opl_vmlmb_workspace_t* ws);
extern double opl_vmlmb_get_sgtol(opl_vmlmb_workspace_t* ws);
extern double opl_vmlmb_get_sxtol(opl_vmlmb_workspace_t* ws);
extern double opl_vmlmb_get_frtol(opl_vmlmb_workspace_t* ws);
extern double opl_vmlmb_get_fatol(opl_vmlmb_workspace_t* ws);
extern double opl_vmlmb_get_epsilon(opl_vmlmb_workspace_t* ws);
extern double opl_vmlmb_get_delta(opl_vmlmb_workspace_t* ws);
extern double opl_vmlmb_get_step(opl_vmlmb_workspace_t* ws);
extern double opl_vmlmb_get_gnorm(opl_vmlmb_workspace_t* ws);

extern opl_integer_t opl_vmlmb_get_evaluations(opl_vmlmb_workspace_t* ws);
extern opl_integer_t opl_vmlmb_get_iterations(opl_vmlmb_workspace_t* ws);
extern opl_integer_t opl_vmlmb_get_restarts(opl_vmlmb_workspace_t* ws);
/*	Query values of current step size along search direction, curent
	iteration number. */

/**
 * @}
 */

/*---------------------------------------------------------------------------*/
/* APPLY BOUND CONSTRAINTS */

extern void opl_bounds_apply(opl_integer_t n, double x[],
                             const double xmin[], const double xmax[]);
/*	Apply bounds constraints to array X.  Input/output array X has N
	elements, XMIN must have N elements (the lower bounds of X) or be
	NULL (no lower bounds), similarly XMAX must have N elements (the
	upper bounds of X) or be NULL (no upper bounds). */

extern void opl_bounds_free(opl_integer_t n, opl_logical_t isfree[],
                            const double x[], const double g[],
                            const double xmin[], const double xmax[]);
/*	Set elements of ISFREE to true or false whether the corresponding
	elements in X belong to the unbinded set of parameters or not.
	Output array ISFREE is an N element array. Input arrays X, XMIN and
	XMAX are the same as in opl_bounds_apply.  Input array G is the
	gradient which is an N element array.  The unbinded set of parameters
	verify the following conditions:

	  (X[i] > XMIN[i] || G[i] < 0) && (X[i] < XMAX[i] || G[i] > 0) */

extern opl_integer_t opl_bounds_check(opl_integer_t n, const double xmin[],
                                      const double xmax[]);
/*	Check correctness of bounds XMIN and XMAX (see opl_bounds_apply for
	the definition of the arguments).  This function returns -1 if the
	bounds are such that XMIN[i] <= XMAX[i] for all i=0,...,N-1;
	otherwise, the function return the value i of the first index (i >=
	0) for which the condition is violated. */

extern void opl_lower_bound_apply(opl_integer_t n, double x[], double xmin);
extern void opl_lower_bound_free(opl_integer_t n, opl_logical_t isfree[],
                                 const double x[], const double g[],
                                 double xmin);
/*	These routines are similar to opl_bounds_apply and opl_bounds_free but
	for a scalar lower bound XMIN that is the same for all parameters X. */

extern void opl_upper_bound_apply(opl_integer_t n, double x[], double xmax);
extern void opl_upper_bound_free(opl_integer_t n, opl_logical_t isfree[],
                                 const double x[], const double g[],
                                 double xmax);
/*	These routines are similar to opl_bounds_apply and opl_bounds_free but
	for a scalar upper bound XMAX that is the same for all parameters X. */

extern void opl_interval_apply(opl_integer_t n, double x[], double a, double b);
extern void opl_interval_free(opl_integer_t n, opl_logical_t isfree[],
                              const double x[], const double g[],
                              double a, double b);
/*	These routines are similar to opl_bounds_apply and opl_bounds_free
	but for a scalar lower bound XMIN=min(A,B) and a scalar upper bound
	XMAX=max(A,B) that are the same for all parameters X. */

/*---------------------------------------------------------------------------*/
/* UTILITIES */

extern double opl_dnrm2(opl_integer_t n, const double x[]);
/* Returns the Euclidian norm of X: sqrt(X'.X), taking care of overflows. */


extern void opl_dcopy(opl_integer_t n, const double x[], double y[]);
extern void opl_dcopy_free(opl_integer_t n, const double x[],
                           double y[], const opl_logical_t isfree[]);
/* Copy elements of X into Y.  Does Y[i] = X[i] for i=0,...,N-1.  If ISFREE
   is non-NULL, only elements for which ISFREE[i] is true (non-zero) are
   taken into account. */

extern void opl_daxpy(opl_integer_t n, double a,
		     const double x[], double y[]);
extern void opl_daxpy_free(opl_integer_t n, double a,
                           const double x[], double y[],
                           const opl_logical_t isfree[]);
/* Does Y[i] += A*X[i] for i=0,...,N-1.  If ISFREE is non-NULL, only
   elements for which ISFREE[i] is true (non-zero) are taken into
   account. */

extern double opl_ddot(opl_integer_t n, const double x[], const double y[]);
extern double opl_ddot_free(opl_integer_t n, const double x[],
                            const double y[], const opl_logical_t isfree[]);
/* Computes dot product of N-element vectors X and Y. If ISFREE is
   non-NULL, only elements for which ISFREE[i] is true (non-zero) are taken
   into account. */

extern void opl_dscal(opl_integer_t n, double a, double x[]);
/* Scales N-element vector X by scalar A. */

/*---------------------------------------------------------------------------*/
/* YORICK-LIKE ROUTINES */

extern int opl_anyof(opl_integer_t n, const double x[]);
/* Returns true (non-zero) if any element of X is non-zero; returns faslse
   (zero) otherwise. N is the number of elements of X. */

extern int opl_noneof(opl_integer_t n, const double x[]);
/* Returns true (non-zero) if all elements of X are zero; returns faslse
   (zero) otherwise. N is the number of elements of X. */

extern int opl_allof(opl_integer_t n, const double x[]);
/* Returns true (non-zero) if all elements of X are non-zero; returns faslse
   (zero) otherwise. N is the number of elements of X. */

/*---------------------------------------------------------------------------*/
_OPL_END_DECLS
#endif /* _OPTIMPACKLEGACY_H */
