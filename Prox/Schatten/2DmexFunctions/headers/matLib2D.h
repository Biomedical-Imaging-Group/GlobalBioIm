#include <string.h> /* needed for memcpy() */
#include <stdlib.h>
#include <stdio.h>
#include <mex.h>
#include <math.h>
#include "matrix.h"
#include "epph.h" // This is the header file for general lp projections
                   

#define MAXITER 100 //Maximum allowed number of iterations.
#define EPS_M 1e-15 //Machine double-point precision.

#ifndef max
#define max(a, b) ((a)>(b)?(a):(b))
#endif

#ifndef min
#define min(a, b) ((a)<(b)?(a):(b))
#endif

#ifndef COND_ABS
#define COND_ABS(a, b) (b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a))
#endif

#ifndef sign
#define sign(a) (a >= 0 ? (a > 0 ? 1 : 0) : -1)
#endif

/* =========================================================================
%
%  Author: stamatis.lefkimmiatis@epfl.ch
%
% =========================================================================*/


/* print a matrix */
void printMat(ptrdiff_t m, ptrdiff_t n,  double *X, char *name)
{
  ptrdiff_t di, dj;
  mexPrintf("Mat %s:\n", name);
  for(di=0; di<m; di++)
  {
    for(dj=0; dj<n; dj++)
      mexPrintf("%7.6e \t", *(X+dj*m+di));
    mexPrintf("\n");
  }
}

/* print a vector */
void printVec(ptrdiff_t m,  double *X, char *name)
{
  ptrdiff_t di;
  mexPrintf("vector %s:\n", name);
  for(di=0; di<m; di++)
    mexPrintf("%3.2e \t", X[di]);
  mexPrintf("\n");
}


void eig2x2(double A[3], double U[2], double E[2]) {
  
  double k;
  double trace_A;
  double delta_A;
  
  if (fabs(A[1]) < 1e-15){
    E[0]=A[0];
    E[1]=A[2];
    U[0]=1.0;
    U[1]=0.0;}
  else
  {
    trace_A=A[0]+A[2];
    delta_A=(A[0]-A[2])*(A[0]-A[2])+4*A[1]*A[1];
    E[0]=0.5*(trace_A+sqrt(delta_A));
    E[1]=0.5*(trace_A-sqrt(delta_A));
    k=sqrt((E[0]-A[0])*(E[0]-A[0])+A[1]*A[1]);
    U[0]=A[1]/k;
    U[1]=(E[0]-A[0])/k;
  }
  
}


double pnorm(double *x, double k, double p) {
//x:vector of size k
//p: the order of the norm
  
  double res=0;
  int i;
  
  if (p < 1)
    mexErrMsgTxt("The order of the norm should be greater or equal to one");
  
  if (mxIsFinite(p)){
    for (i=0; i < k; i++)
      res+=pow(fabs(x[i]), p);
    res=pow(res, 1.0/p);}
  else{
    for (i=0; i < k; i++)
      res=max(res, fabs(x[i]));}
  
  return res;
}

double rootfind(double * x_, double *c_, int *iter_step, double * v_, int k_, double p_, double c0, double tau, double x1, double x2, double tol)
/*Using Brent's method, find the root of the (epp-tau) function known to lie between x1 and x2. The
 * root, returned as rootfind, will be refined until its accuracy is tol.*/
{
  int iter;
  double a=x1, b=x2, c=x2, d, e, min1, min2;
  double fa, fb;
  
  epp(x_, c_, iter_step, v_, k_, a, p_, c0);
  fa=pnorm(x_, k_, p_)-tau;
  
  
  epp(x_, c_, iter_step, v_, k_, b, p_, c0);
  fb=pnorm(x_, k_, p_)-tau;
  
  double fc, p, q, r, s, tol1, xm;
  
  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
    mexErrMsgTxt("The specified interval doesn't contain a root");
  fc=fb;
  for (iter=1;iter<=MAXITER;iter++) {
    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
      c=a; // Rename a, b, c and adjust bounding interval
      fc=fa; //d.
      e=d=b-a;
    }
    if (fabs(fc) < fabs(fb)) {
      a=b;
      b=c;
      c=a;
      fa=fb;
      fb=fc;
      fc=fa;
    }
    tol1=2.0*EPS_M*fabs(b)+0.5*tol; //Convergence check.
    xm=0.5*(c-b);
    if (fabs(xm) <= tol1 || fb == 0.0) return b;
    if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
      s=fb/fa; //Attempt inverse quadratic interpolation.
      if (a == c) {
        p=2.0*xm*s;
        q=1.0-s;
      } else {
        q=fa/fc;
        r=fb/fc;
        p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
        q=(q-1.0)*(r-1.0)*(s-1.0);
      }
      if (p > 0.0) q = -q; //Check whether in bounds.
      p=fabs(p);
      min1=3.0*xm*q-fabs(tol1*q);
      min2=fabs(e*q);
      if (2.0*p < (min1 < min2 ? min1 : min2)) {
        e=d; //Accept interpolation.
        d=p/q;
      } else {
        d=xm; //Interpolation failed, use bisection.
        e=d;
      }
    } else { //Bounds decreasing too slowly, use bisection.
      d=xm;
      e=d;
    }
    a=b; //Move last best guess to a.
    fa=fb;
    if (fabs(d) > tol1) //Evaluate new trial root.
      b += d;
    else
      b += COND_ABS(tol1, xm);
    
    epp(x_, c_, iter_step, v_, k_, b, p_, c0);
    fb=pnorm(x_, k_, p_)-tau;
  }
  mexErrMsgTxt("Maximum number of iterations exceeded in rootfind");
  return 0.0; //Never get here.
}


void projectS2(double *Xp, double *X, double rho){
  int j;
  double normF;
  
  //Frobenius norm of X
  normF=sqrt(X[0]*X[0]+2*X[1]*X[1]+X[2]*X[2]);
  
  if (normF <= rho){
    //matrix reconstruction
    memcpy(Xp, X, 3*sizeof(double));}
  else{
    //matrix projection
    for (j=0; j < 3; j++)
      Xp[j]=(X[j]/normF)*rho;
  }
}


void projectSinf(double *Xp, double *X, double rho){
  
  double norm_linf;
  double E[2]; // Eigenvalues
  double Ep[2];// Projected Eigenvalues
  double V[2]; // Eigenvector
  
  eig2x2(X,V,E);
  
  norm_linf=max(fabs(E[0]),fabs(E[1]));//Infinity norm.
  
  if (norm_linf <= rho){
    //matrix reconstruction
    memcpy(Xp, X, 3*sizeof(double));}
  else{
    //matrix projection
    //Projection of the eigenvalues
    Ep[0]=sign(E[0])*min(fabs(E[0]), rho);
    Ep[1]=sign(E[1])*min(fabs(E[1]), rho);
    //matrix reconstruction
    Xp[0]=V[0]*V[0]*Ep[0]+V[1]*V[1]*Ep[1];
    Xp[1]=V[0]*V[1]*(Ep[0]-Ep[1]);
    Xp[2]=V[0]*V[0]*Ep[1]+V[1]*V[1]*Ep[0];
  }
}


void projectS1(double *Xp, double *X, double rho){
  
  double norm_l1;
  double E[2]; // Eigenvalues
  double Ep[2];// Projected Eigenvalues
  double V[2]; // Eigenvector
  double S[2]; // Singular Values
  
  eig2x2(X,V,E);
  
  S[0]=fabs(E[0]);
  S[1]=fabs(E[1]);
  
  norm_l1=S[0]+S[1];//Nuclear norm.
  
  if (norm_l1 <= rho){
    //matrix reconstruction
    memcpy(Xp, X, 3*sizeof(double));}
  else{
    //matrix projection
    //Projection of the eigenvalues
    double gamma;  
    if (rho < fabs(S[0]-S[1]))
      gamma=max(S[0],S[1])-rho;
    else
      gamma=(norm_l1-rho)/2;
        
    Ep[0]=sign(E[0])*max(S[0]-gamma, 0.0);
    Ep[1]=sign(E[1])*max(S[1]-gamma, 0.0);
    
    //matrix reconstruction
    Xp[0]=V[0]*V[0]*Ep[0]+V[1]*V[1]*Ep[1];
    Xp[1]=V[0]*V[1]*(Ep[0]-Ep[1]);
    Xp[2]=V[0]*V[0]*Ep[1]+V[1]*V[1]*Ep[0];
    
  }
}

void projectSp(double *Xp, double *X, double rho, double p, double * c, double c0){
  
  double norm_lp;
  double E[2]; // Eigenvalues
  double Ep[2];// Projected Eigenvalues
  double V[2]; // Eigenvector
  
  eig2x2(X,V,E);
  
  norm_lp=pow(fabs(E[0]),p)+pow(fabs(E[1]),p);//lp norm.
  norm_lp=pow(norm_lp,1.0/p);
  
  if (norm_lp <= rho){
    //matrix reconstruction
    memcpy(Xp, X, 3*sizeof(double));}
  else{
    //matrix projection
    //Projection of the eigenvalues -- projectlp(Ep,E,2,rho,p,c,c0);
    int steps[2];
    double q=p/(p-1);
    double norm_lq=pow(fabs(E[0]),q)+pow(fabs(E[1]),q);//lq norm.
    norm_lq=pow(norm_lq,1.0/q);
    double rho_opt=rootfind(Ep, c, steps, E, 2, p, c0, rho, 0, norm_lq, 1e-8);
    epp(Ep, c, steps, E, 2, rho_opt, p, c0);
    
    //matrix reconstruction
    Xp[0]=V[0]*V[0]*Ep[0]+V[1]*V[1]*Ep[1];
    Xp[1]=V[0]*V[1]*(Ep[0]-Ep[1]);
    Xp[2]=V[0]*V[0]*Ep[1]+V[1]*V[1]*Ep[0];
  }
}


