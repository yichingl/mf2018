#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "header.h"

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


int main(int argc, char **argv)
{
    int i, j;
    int ndim = 10;
    float *eig = (float *) calloc(ndim, sizeof(float));
    float **eigvec = (float **)alloc2d(ndim,ndim,sizeof(float));
    float **a = (float **)alloc2d(ndim,ndim,sizeof(float));
    float **w;

    for (i=0; i<ndim; i++) for (j=0; j<=i; j++) a[i][j] = a[j][i] = i-j;
    for (i=0; i<ndim; i++) a[i][i] += 10;

    printf("A: \n");
    for (i=0; i<ndim; i++) {
        for (j=0; j<ndim; j++) printf("%f ", a[i][j]);
        printf("\n");
    }

    eigendecompose(a, ndim, eig, eigvec);
    printf("EIGENVALUES: \n");
    for(i=0;i<ndim;i++) printf("%e ",eig[i]); printf("\n");
    printf("EIGENVECTORS \n");
    for (i=0; i<ndim; i++) {
        for (j=0; j<ndim; j++) printf("%7.5f ", eigvec[j][i]);//print row x col
        printf("\n");
    }
    w = whiteningmatrix(a, ndim, ndim);
    printf("WHITENING MATRIX \n");
    for (i=0; i<ndim; i++) {
        for (j=0; j<ndim; j++) printf("%7.5f ", w[i][j]);
        printf("\n");
    }
}
