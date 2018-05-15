#include <ayesha.h>

int
invert(float **ainv,
       float **a,
       int len)
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

int ludcmp(float **a,
             int   n,
             int   *indx,
             float *d)
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


int lubksb(
       float **a, 
       int   n, 
       int   *indx, 
       float b[],
       float x[]
      )
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
