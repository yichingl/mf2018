#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <gsl/gsl_cblas.h>
#include "header.h"

float **matmult(float **a, int rowsa, int colsa,
                float **b, int rowsb, int colsb,
                float **c)
{
    assert(colsa==rowsb);
    if (c == NULL) c = (float**) alloc2d(colsb, rowsa, sizeof(float));
    cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                rowsa, colsb, colsa,
                1, a[0], colsa, b[0], colsb, 0, c[0], colsb);
    return c;
}

float **matmulttxpA(float **a, int rowsa, int colsa,
                    float **b, int rowsb, int colsb,
                    float **c)
{
    assert(rowsa==rowsb);
    if (c==NULL) c = (float **) alloc2d(colsa, colsb, sizeof(float));

    cblas_sgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                colsa, colsb, rowsa,
                1, a[0], colsa, b[0], colsb, 0, c[0], colsb);
    return c;
}

float **matmulttxpB(float **a, int rowsa, int colsa,
                    float **b, int rowsb, int colsb,
                    float **c)
{
    assert(colsa==colsb);
    if (c==NULL) c = (float **) alloc2d(rowsa, rowsb, sizeof(float));

    cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                rowsa, rowsb, colsa,
                1, a[0], colsa, b[0], colsb, 0, c[0], rowsb);
    return c;
}

/* Matrix inversion from Numerical Recipes in C */

int ludcmp(float **a, int   n, int   *indx, float *d)
{
    int    i,imax,j,k;
    float  big,dum,sum,t1;
    float  *vv;  	     /* vv stores the implicit scaling of each row */

    if((vv=(float *) calloc(n,sizeof(float))) == NULL)
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


int lubksb( float **a, int   n, int   *indx, float b[], float x[])
{
    int   i,ii=0,ip,j,done=0;
    float sum;

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


int
matrixinvert(float **ainv, float **a, int len)
{
    int i, j;
    int *indx;
    float d;
    float *col;
    float **adcmp;

    int  ludcmp();
    int  lubksb();

    indx = calloc(len, sizeof(int));
    col = calloc(len, sizeof(float));
    adcmp = (float **)alloc2d(len, len, sizeof(float));

    for (i = 0; i < len; i++) {
	for (j = 0; j < len; j++) {
	    adcmp[i][j] = a[i][j];
	}
    }

    ludcmp(adcmp, len, indx, &d);
    for (j = 0; j < len; j++) {
	for (i = 0; i < len; i++)
	    col[i] = 0;
	col[j] = 1;

	lubksb(adcmp, len, indx, col, col);
	for (i = 0; i < len; i++) {
	    ainv[i][j] = col[i];
	}
    }

    free(indx);
    free(col);
    free2d((void **)adcmp);

    return 0;
}


/*------------------------------------------------------------------------
Eigen decomposition routine!

Using numerical recipes jacobi routine to compute eigen values/vectors 
--------------------------------------------------------------------------*/

void jacobi(float **a, int n, float d[], float **v, int *nrot);

/* Stub routine to interface to NR code which has 1-to-n indexing */
int eigendecompose(float **a, int n, float *eigval, float **eigvec)
{
    float **localA, **localeigvec, *localeigval, maxeig;
    int i,j,nrot,maxid;

    localA = (float **)alloc2d(n+1,n+1,sizeof(float));   
    localeigvec = (float **)alloc2d(n+1,n+1,sizeof(float));   
    localeigval = (float *)calloc(n+1,sizeof(float));   

    for (i=0;i<n;i++) 
	for (j=0;j<n;j++) 
	    localA[i+1][j+1] = a[i][j];

    jacobi(localA,n,localeigval,localeigvec,&nrot);

    for (i=0;i<n;i++) {
        // Sort by eigenvalue: find max eig. Using stupidsort.
        maxeig = localeigval[1]; maxid = 1;
        for (j=2;j<=n;j++) 
            if (localeigval[j] > maxeig) {maxeig = localeigval[j]; maxid = j;}
        localeigval[maxid] = -1e+30;  // Kill this eigval
	for (j=0;j<n;j++) 
	    eigvec[i][j] = localeigvec[j+1][maxid];
	eigval[i] = maxeig;
    }

    free(localeigval);
    free2d((void**)localA);
    free2d((void**)localeigvec);

    return nrot;
}


#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	a[k][l]=h+s*(g-h*tau);

/*------------------------------------------------------------------------
Computes all eigenvalues and eigenvectors of a real symmetric matrix 
a[1..n][1..n]. On output, elements of a above the diagonal are destroyed. 
d[1..n] returns the eigenvalues of a.
v[1..n][1..n] is a matrix whose columns contain, on output, the normalized 
              eigenvectors of a. 
nrot returns the number of Jacobi rotations that were required.
--------------------------------------------------------------------------*/

void nrerror(char *error_text)
{
	void exit();

	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

float *vector(int nl,int nh)
{
	float *v;

	v=(float *)malloc((unsigned) (nh-nl+1)*sizeof(float));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl;
}

void free_vector(float *v, int nl, int nh)
{
	free((char*) (v+nl));
}


void jacobi(float **a, int n, float d[], float **v, int *nrot)
{
	int j,iq,ip,i;
	float tresh,theta,tau,t,sm,s,h,g,c,*b,*z,*vector();
	void nrerror(),free_vector();

	b=vector(1,n);
	z=vector(1,n);
	for (ip=1;ip<=n;ip++) {
		for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
		v[ip][ip]=1.0;
	}
	for (ip=1;ip<=n;ip++) {
		b[ip]=d[ip]=a[ip][ip];
		z[ip]=0.0;
	}
	*nrot=0;
	for (i=1;i<=50;i++) {
		sm=0.0;
		for (ip=1;ip<=n-1;ip++) {
			for (iq=ip+1;iq<=n;iq++)
				sm += fabs(a[ip][iq]);
		}
		if (sm == 0.0) {
			free_vector(z,1,n);
			free_vector(b,1,n);
			return;
		}
		if (i < 4)
			tresh=0.2*sm/(n*n);
		else
			tresh=0.0;
		for (ip=1;ip<=n-1;ip++) {
			for (iq=ip+1;iq<=n;iq++) {
				g=100.0*fabs(a[ip][iq]);
				if (i > 4 && (float)(fabs(d[ip])+g) == (float)fabs(d[ip])
					&& (float)(fabs(d[iq])+g) == (float)fabs(d[iq]))
					a[ip][iq]=0.0;
				else if (fabs(a[ip][iq]) > tresh) {
					h=d[iq]-d[ip];
					if ((float)(fabs(h)+g) == (float)fabs(h))
						t=(a[ip][iq])/h;
					else {
						theta=0.5*h/(a[ip][iq]);
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq]=0.0;
					for (j=1;j<=ip-1;j++) {
						ROTATE(a,j,ip,j,iq)
					}
					for (j=ip+1;j<=iq-1;j++) {
						ROTATE(a,ip,j,j,iq)
					}
					for (j=iq+1;j<=n;j++) {
						ROTATE(a,ip,j,iq,j)
					}
					for (j=1;j<=n;j++) {
						ROTATE(v,j,ip,j,iq)
					}
					++(*nrot);
				}
			}
		}
		for (ip=1;ip<=n;ip++) {
			b[ip] += z[ip];
			d[ip]=b[ip];
			z[ip]=0.0;
		}
	}
	nrerror("Too many iterations in routine JACOBI");
}

#undef ROTATE

/*----------------------------------------------------------------------
 Use Eigen decomposition to compute whitening matrix
 Assumes input is a covariance matrix; otherwise bad things happen
------------------------------------------------------------------------*/
float **whiteningmatrix(float **cov, int dim, int projecteddim)
{
    int i, j;
    float *eigval = (float *)calloc(dim, sizeof(float));
    float **eigvec = (float **) alloc2d(dim, dim, sizeof(float));
    float **whiteningmat = (float **) alloc2d(projecteddim, dim, sizeof(float));

    eigendecompose(cov, dim, eigval, eigvec);

    // The whitening matrix is eigvec' * sqrt(eigval)
    for (i=0; i<projecteddim; i++) {
        for (j=0; j<dim; j++) 
            whiteningmat[i][j] = sqrt(fabs(eigval[i])) * eigvec[i][j];
    }
    free2d((void**) eigvec); free(eigval);

    return(whiteningmat);
}


/* A variant which cuts off dimensions with very low eigen values */
void computewhiteningmatrix(float **covar, int dim, float threshold, 
                            int *outdim, float ***whitener)
{
    int i, j, odim;
    float *eigval = (float *)calloc(dim, sizeof(float));
    float **eigvec = (float **) alloc2d(dim, dim, sizeof(float));
    float **whiteningmat;

    eigendecompose(covar, dim, eigval, eigvec);
    odim = 1;
    for (i=1; i<dim; i++) {
        if (eigval[i] <= eigval[0] * threshold) break;
        odim++;
    }
    *outdim = odim;

    whiteningmat = (float **) alloc2d(odim, dim, sizeof(float));
    *whitener = whiteningmat;

    // The whitening matrix is eigvec' * sqrt(eigval)
    for (i=0; i<odim; i++) {
        for (j=0; j<dim; j++) 
            whiteningmat[i][j] = sqrt(fabs(eigval[i])) * eigvec[i][j];
    }
    free2d((void**) eigvec); free(eigval);

}

