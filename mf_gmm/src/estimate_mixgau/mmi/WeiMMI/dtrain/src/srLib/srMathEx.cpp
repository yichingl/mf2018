#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <memory.h>

#include "srShellEx.h"
#include "srMathEx.h"


/* Choleski for symmetric matrix: 
/* Place lower triangular choleski factor of A in L.*/
/* Return 0 if matrix singular or not +definite */
int Cholesky(double **A, double **L, int size)
{
	int i,j,k;
	double sum;

	for (i = 0; i < size; i++) {
		for (j = 0; j <= i; j++) {
			sum = A[i][j];
			for (k = 0; k < j; k++)
				sum -= (L[i][k]*L[j][k]);
			if ((i == j) && (sum <= 0.0)) 
				return 0;
			else if (i == j)
				sum = sqrt(sum);
			else if (L[j][j] == 0.0)
				return 0;
			else
				sum /= L[j][j];
			L[i][j] = sum;
		}
	}
	for (i = 0; i < size; i++) 
		for (j = i + 1; j < size; j++) { 
			L[i][j] = 0.0;
		}
	return 1;
}


/* MSolve: solve Ly=e^i and L^t x = y, where e^i is a unit vector */
// nr - size of square matrix L
// i   - 0 ~ nr
void MSolve(double **L, int i, double *x, double *y, int nr)
{
	int j, k;
	double sum;
	
	// solve Ly=e^i
	for (j = 0; j < i; j++) y[j] = 0.0;  /* forward sub */
	y[i] = 1.0 / L[i][i];
	for (j = i+1; j < nr; j++){
		sum = 0.0;
		for (k = i; k < j; k++)
			sum -= L[j][k]*y[k];
		y[j] = sum / L[j][j];
	}
	
	//solve L^t x = y,
	x[nr - 1] = y[nr - 1] / L[nr - 1][nr - 1];         /* backward sub */
	for (j = nr - 2; j > 0; j--) {
		sum = y[j];
		for (k = j + 1; k <= nr; k++)
			sum -= L[k][j] * x[k];
		x[j] = sum / L[j][j];
	}
}


double Determinant(double **a,int n)
{
	int i,j,j1,j2;
	double det = 0;
	double **m = NULL;

	if (n < 1) { /* Error */

	} else if (n == 1) { /* Shouldn't get used */
		det = a[0][0];
	} else if (n == 2) {
		det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
	} else {
		det = 0;
		for (j1=0;j1<n;j1++) {
			m = (double **)malloc((n-1)*sizeof(double *));
			for (i=0;i<n-1;i++)
				m[i] = (double *)malloc((n-1)*sizeof(double));
			for (i=1;i<n;i++) {
				j2 = 0;
				for (j=0;j<n;j++) {
					if (j == j1)
						continue;
					m[i-1][j2] = a[i][j];
					j2++;
				}
			}
			det += pow(-1.0,1.0+j1+1.0) * a[0][j1] * Determinant(m,n-1);
			for (i=0;i<n-1;i++)
				free(m[i]);
			free(m);
		}
	}
	if (det <  DSMALL)
		det = DSMALL;
	return(det);
}

/* 1 - invc = c^-1, return log(Det(c)) */
// c must be positive definite
double InvertCovStpd(double** c, double **invc, int n)
{
	double ldet = .0;
	int i,j;

	double **l = (double **) alloc2d(n, n, sizeof(double));
	double *x = (double *) calloc(n, sizeof(double));
	double *y = (double *) calloc(n, sizeof(double));
	if (Cholesky(c, l, n)){
		for (j = 0; j < n; j++){
			MSolve(l , j, x, y, n);
			for (i = 0; i < n; i++)
				invc[i][j] = x[i];
			ldet += log(l[j][j]);
		}
	} else {
		printf("CovDet: [%f ...] not invertible", c[1][1]);
	}
	free2d((char **)l);
	free(x);
	free(y);
	return 2.0*ldet;
}

double InvertCov(double** c, double **invc, int n)
{
	double ldet = 0.0;
	int i,j;

	minvert(c, invc, n);
	ldet = Determinant(c, n);
	
	return log(ldet);
}

// return log(det(A)) = 2 * log(det(L)) where A = LL'
double LogCovDet(double** c, int n)
{
	int j;
	double ldet = 0.0;
	double **l = (double **) alloc2d(n, n, sizeof(double));
	if (Cholesky(c, l, n)) {
		for (j = 0; j < n; j++)
			ldet += log(l[j][j]);
	} else {
		printf("CovDet: [%f ...] not invertible", c[1][1]);
	}
	free2d((char **)l);
	return 2.0*ldet;
}

int minvert(double **a, double **ainv, int len)
{
	int i, j;
	int *indx;
	double d;
	double *col;
	double **adcmp;

//	int  ludcmp_dbl();
//	int  lubksb_dbl();

	indx = (int *)calloc(len, sizeof(int));
	col = (double *)calloc(len, sizeof(double));
	adcmp = (double **)alloc2d(len, len, sizeof(double));

	for (i = 0; i < len; i++) {
		for (j = 0; j < len; j++) {
			adcmp[i][j] = a[i][j];
		}
	}

	ludcmp_dbl(adcmp, len, indx, &d);
	for (j = 0; j < len; j++) {
		for (i = 0; i < len; i++)
			col[i] = 0;
		col[j] = 1;

		lubksb_dbl(adcmp, len, indx, col, col);
		for (i = 0; i < len; i++) {
			ainv[i][j] = col[i];
		}
	}

	free(indx);
	free(col);
	free2d((char **)adcmp);

	return 0;
}

int ludcmp_dbl(double **a,
	       int   n,
	       int   *indx,
	       double *d)
{
	int    i,imax,j,k;
	double  big,dum,sum,t1;
	double  *vv;  	     /* vv stores the implicit scaling of each row */

	if((vv=(double *) calloc(n,sizeof(double))) == NULL)
	{
		printf("Unable to allocate space for vv\n");fflush(stdout);
		exit(0);
	}

		      *d=1.0;           /* No row interchanges yet */
		      for (i = 0; i < n; i++) { /*Loop over the rows to get the implicit scaling */
			      big = 0.0;    /*information */
			      for (j = 0; j < n; j++) {
				      if((t1 = fabs(a[i][j])) > big)
					      big=t1;
			      }
			      if (big == 0.0) {
				      /* No nonzero largest element. */
				      printf("ERROR: Singular matrix in routine ludcmp\n");fflush(stdout);
			      }
			      vv[i] = 1.0/big; /* Save the scaling */
		      }
		      for (j = 0; j < n; j++) {
			      for (i = 0; i < j; i++) {
				      sum = a[i][j];
				      for (k = 0; k < i; k++)
					      sum -= a[i][k]*a[k][j];
				      a[i][j] = sum;
			      }
			      big = 0.0;
			      for (i = j; i < n;i++) {
				      sum = a[i][j];
				      for (k = 0; k < j; k++)
					      sum -= a[i][k] * a[k][j];
				      a[i][j] = sum;
				      if ( (dum = vv[i] * fabs(sum)) >= big) {
					      big = dum;
					      imax = i;
				      }
			      }
			      if (j != imax) {
				      for (k = 0; k < n; k++) {
					      dum = a[imax][k];
					      a[imax][k] = a[j][k];
					      a[j][k] = dum;
				      }
				      *d = -(*d);
				      vv[imax] = vv[j];
			      }
			      indx[j] = imax;
			      if (a[j][j] == 0.0)
				      a[j][j]=1.0e-20;

			      if (j != n-1) {
				      dum = 1.0 / a[j][j];
				      for (i = j+1; i < n; i++)
					      a[i][j] *= dum;
			      }
		      }
    
		      free(vv);

		      return 0;
}


int lubksb_dbl(
		double **a, 
       int   n, 
       int   *indx, 
       double b[],
       double x[]
	      )
{
	int   i,ii=0,ip,j,done=0;
	double sum;

	for (i = 0; i < n; i++) {
		ip = indx[i];
		sum = b[ip];
		b[ip] = b[i];
		if (done) {
			for (j = ii; j < i; j++)
				sum -= a[i][j] * b[j];
		}
		else if (sum) {
			ii=i;
			done=1;
		}
		b[i]=sum;
	}
	for (i=0;i<n;i++) {
		x[i] = b[i];
	}
    
	for (i = n-1; i >= 0; i--) {
		sum = x[i];
	
		for (j = n-1; j > i; j--)
			sum -= a[i][j]*x[j];
	
		x[i]=sum/a[i][i];
	}

	return 0;
}

