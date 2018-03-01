#include <string.h> /* needed for memcpy() */
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
 * %
 * %  Author: stamatis.lefkimmiatis@epfl.ch
 * %
 * % =========================================================================*/


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


static double hypot2(double x, double y) {
  return sqrt(x*x+y*y);
}

// Symmetric Householder reduction to tridiagonal form.

static void tred2(double V[9], double d[3], double e[3]) {
  
//  This is derived from the Algol procedures tred2 by
//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
//  Fortran subroutine in EISPACK.
  
  int i, j, k;
  double f, g, h, hh;
  for (j = 0; j < 3; j++) {
    d[j] = V[2+3*j];
  }
  
  // Householder reduction to tridiagonal form.
  
  for (i = 2; i > 0; i--) {
    
    // Scale to avoid under/overflow.
    
    double scale = 0.0;
    double h = 0.0;
    for (k = 0; k < i; k++) {
      scale = scale + fabs(d[k]);
    }
    if (scale == 0.0) {
      e[i] = d[i-1];
      for (j = 0; j < i; j++) {
        d[j] = V[i-1+3*j];
        V[i+3*j] = 0.0;
        V[j+3*i] = 0.0;
      }
    } else {
      
      // Generate Householder vector.
      
      for (k = 0; k < i; k++) {
        d[k] /= scale;
        h += d[k] * d[k];
      }
      f = d[i-1];
      g = sqrt(h);
      if (f > 0) {
        g = -g;
      }
      e[i] = scale * g;
      h = h - f * g;
      d[i-1] = f - g;
      for (j = 0; j < i; j++) {
        e[j] = 0.0;
      }
      
      // Apply similarity transformation to remaining columns.
      
      for (j = 0; j < i; j++) {
        f = d[j];
        V[j+3*i] = f;
        g = e[j] + V[j+3*j] * f;
        for (k = j+1; k <= i-1; k++) {
          g += V[k+3*j] * d[k];
          e[k] += V[k+3*j] * f;
        }
        e[j] = g;
      }
      f = 0.0;
      for (j = 0; j < i; j++) {
        e[j] /= h;
        f += e[j] * d[j];
      }
      hh = f / (h + h);
      for (j = 0; j < i; j++) {
        e[j] -= hh * d[j];
      }
      for (j = 0; j < i; j++) {
        f = d[j];
        g = e[j];
        for (k = j; k <= i-1; k++) {
          V[k+3*j] -= (f * e[k] + g * d[k]);
        }
        d[j] = V[i-1+3*j];
        V[i+3*j] = 0.0;
      }
    }
    d[i] = h;
  }
  
  // Accumulate transformations.
  
  for (i = 0; i < 2; i++) {
    V[2+3*i] = V[4*i];
    V[4*i] = 1.0;
    h = d[i+1];
    if (h != 0.0) {
      for (k = 0; k <= i; k++) {
        d[k] = V[k+3*(i+1)] / h;
      }
      for (j = 0; j <= i; j++) {
        g = 0.0;
        for (k = 0; k <= i; k++) {
          g += V[k+3*(i+1)] * V[k+3*j];
        }
        for (k = 0; k <= i; k++) {
          V[k+3*j] -= g * d[k];
        }
      }
    }
    for (k = 0; k <= i; k++) {
      V[k+3*(i+1)] = 0.0;
    }
  }
  for (j = 0; j < 3; j++) {
    d[j] = V[2+3*j];
    V[2+3*j] = 0.0;
  }
  V[8] = 1.0;
  e[0] = 0.0;
}

// Symmetric tridiagonal QL algorithm.

static void tql2(double V[9], double d[3], double e[3]) {
  
//  This is derived from the Algol procedures tql2, by
//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
//  Fortran subroutine in EISPACK.
  
  int i, j, m, l, k;
  double g, p, r, dl1, h, f, tst1, eps;
  double c, c2, c3, el1, s, s2;
  
  for (i = 1; i < 3; i++) {
    e[i-1] = e[i];
  }
  e[2] = 0.0;
  
  f = 0.0;
  tst1 = 0.0;
  eps = pow(2.0, -52.0);
  for (l = 0; l < 3; l++) {
    
    // Find small subdiagonal element
    
    tst1 = max(tst1, fabs(d[l]) + fabs(e[l]));
    m = l;
    while (m < 3) {
      if (fabs(e[m]) <= eps*tst1) {
        break;
      }
      m++;
    }
    
    // If m == l, d[l] is an eigenvalue,
    // otherwise, iterate.
    
    if (m > l) {
      int iter = 0;
      do {
        iter = iter + 1;  // (Could check iteration count here.)
        
        // Compute implicit shift
        
        g = d[l];
        p = (d[l+1] - g) / (2.0 * e[l]);
        r = hypot2(p, 1.0);
        if (p < 0) {
          r = -r;
        }
        d[l] = e[l] / (p + r);
        d[l+1] = e[l] * (p + r);
        dl1 = d[l+1];
        h = g - d[l];
        for (i = l+2; i < 3; i++) {
          d[i] -= h;
        }
        f = f + h;
        
        // Implicit QL transformation.
        
        p = d[m];
        c = 1.0;
        c2 = c;
        c3 = c;
        el1 = e[l+1];
        s = 0.0;
        s2 = 0.0;
        for (i = m-1; i >= l; i--) {
          c3 = c2;
          c2 = c;
          s2 = s;
          g = c * e[i];
          h = c * p;
          r = hypot2(p, e[i]);
          e[i+1] = s * r;
          s = e[i] / r;
          c = p / r;
          p = c * d[i] - s * g;
          d[i+1] = h + s * (c * g + s * d[i]);
          
          // Accumulate transformation.
          
          for (k = 0; k < 3; k++) {
            h = V[k+3*(i+1)];
            V[k+3*(i+1)] = s * V[k+3*i] + c * h;
            V[k+3*i] = c * V[k+3*i] - s * h;
          }
        }
        p = -s * s2 * c3 * el1 * e[l] / dl1;
        e[l] = s * p;
        d[l] = c * p;
        
        // Check for convergence.
        
      } while (fabs(e[l]) > eps*tst1);
    }
    d[l] = d[l] + f;
    e[l] = 0.0;
  }
  
  // Sort eigenvalues and corresponding vectors.
  
  for (i = 0; i < 2; i++) {
    k = i;
    p = d[i];
    for (j = i+1; j < 3; j++) {
      if (d[j] < p) {
        k = j;
        p = d[j];
      }
    }
    if (k != i) {
      d[k] = d[i];
      d[i] = p;
      for (j = 0; j < 3; j++) {
        p = V[j+3*i];
        V[j+3*i] = V[j+3*k];
        V[j+3*k] = p;
      }
    }
  }
}



void eigensym3x3(double A[6], double V[9], double d[3]) {
  int i;
  double e[3];
  
  /* V[i][j]=V[i+n*j]
   * V[0][0]=A[0] V[0][1]=A[1] V[0][2]=A[2]
   * V[1][0]=A[1] V[1][1]=A[3] V[1][2]=A[4]
   * V[2][0]=A[2] V[2][1]=A[4] V[2][2]=A[5]
   */
  for (i = 0; i < 3; i++)
    V[i] = A[i];
  for (i = 4; i < 6; i++)
    V[i] = A[i-1];
  for (i = 7; i < 9; i++)
    V[i] = A[i-3];
  V[3]=A[1];
  V[6]=A[2];
  
  tred2(V, d, e);
  tql2(V, d, e);
}


void eigen3x3SymRec(double X[6], double V[9], double E[3]){
  X[0]=V[0]*V[0]*E[0]+V[3]*V[3]*E[1]+V[6]*V[6]*E[2];
  X[1]=V[0]*V[1]*E[0]+V[3]*V[4]*E[1]+V[6]*V[7]*E[2];
  X[2]=V[0]*V[2]*E[0]+V[3]*V[5]*E[1]+V[6]*V[8]*E[2];
  X[3]=V[1]*V[1]*E[0]+V[4]*V[4]*E[1]+V[7]*V[7]*E[2];
  X[4]=V[1]*V[2]*E[0]+V[4]*V[5]*E[1]+V[7]*V[8]*E[2];
  X[5]=V[2]*V[2]*E[0]+V[5]*V[5]*E[1]+V[8]*V[8]*E[2];
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
  normF=sqrt(X[0]*X[0]+2*X[1]*X[1]+2*X[2]*X[2]+X[3]*X[3]+2*X[4]*X[4]+X[5]*X[5]);
  
  if (normF <= rho){
    //matrix reconstruction
    memcpy(Xp, X, 6*sizeof(double));}
  else{
    //matrix projection
    for (j=0; j < 6; j++)
      Xp[j]=(X[j]/normF)*rho;
  }
}


void proxS2(double *Xp, double *X, double rho){
  int j;
  double normF=sqrt(X[0]*X[0]+X[3]*X[3]+X[5]*X[5]+2*(X[1]*X[1]+X[2]*X[2]+X[4]*X[4]));
  if (normF==0)
    for (j=0; j < 6; j++)
      Xp[j]=0;
  else{
    double k=max(normF-rho, 0)/normF;
    for (j=0; j < 6; j++)
      Xp[j]=X[j]*k;}
}


void projectSinf(double *Xp, double *X, double rho, double p){
  
  double norm_linf;
  double E[3]; // Eigenvalues
  double Ep[3];// Projected Eigenvalues
  double V[9]; // Eigenvectors
  
  eigensym3x3(X,V,E);
  
  norm_linf=pnorm(E,3,p);//Infinity norm.
  
  if (norm_linf <= rho){
    //matrix reconstruction
    memcpy(Xp, X, 6*sizeof(double));}
  else{
    //matrix projection
    //Projection of the eigenvalues
    Ep[0]=sign(E[0])*min(fabs(E[0]), rho);
    Ep[1]=sign(E[1])*min(fabs(E[1]), rho);
    Ep[2]=sign(E[2])*min(fabs(E[2]), rho);
    //matrix reconstruction
    eigen3x3SymRec(Xp, V, Ep);
  }
}


void projectSp(double *Xp, double *X, double rho, double p, double *c, double c0){
  
  double norm_lp;
  double E[3]; // Eigenvalues
  double Ep[3];// Projected Eigenvalues
  double V[9]; // Eigenvector
  
  eigensym3x3(X,V,E);
  
  norm_lp=pnorm(E,3,p);
  
  if (norm_lp <= rho){
    //matrix reconstruction
    memcpy(Xp, X, 6*sizeof(double));}
  else{
    //matrix projection
    //Projection of the eigenvalues -- projectlp(Ep,E,2,rho,p,c,c0);
    int steps[2];
    double q=p/(p-1);
    double norm_lq=pnorm(E,3,q);//lq norm.
    double rho_opt=rootfind(Ep, c, steps, E, 3, p, c0, rho, 0, norm_lq, 1e-8);
    epp(Ep, c, steps, E, 3, rho_opt, p, c0);
    
    //matrix reconstruction
    eigen3x3SymRec(Xp, V, Ep);
  }
}

void proxSp(double *Xp, double *X, double rho, double p, double c0){
  
  double E[3]; // Eigenvalues
  double Ep[3];// Projected Eigenvalues
  double V[9]; // Eigenvector
  
  eigensym3x3(X,V,E);
  
  int steps[2];
  double c[1];
  epp(Ep, c, steps, E, 3, rho, p, c0);
  
  //matrix reconstruction
  eigen3x3SymRec(Xp, V, Ep);
}
