#ifndef _SR_MATH_EX
#define _SR_MATH_EX

#define LOGZERO	(-1.0E20)   		/* ~log(0) */
#define LSMALL		(-0.5E20)   		/* log values < LSMALL are set to LOGZERO */
#define DSMALL	(1E-3)   			/* A small value */

int Cholesky(double **A, double **L, int size);

double InvertCov(double** c, double **invc, int n);

double LogCovDet(double** c, int n);
double Determinant(double **a,int n);

int minvert(double **a, double **ainv, int len);

int ludcmp_dbl(double **a,
	       int   n,
	       int   *indx,
	       double *d);

int lubksb_dbl(
		double **a, 
       int   n, 
       int   *indx, 
       double b[],
       double x[]
	      );
#endif //_SR_MATH_EX
