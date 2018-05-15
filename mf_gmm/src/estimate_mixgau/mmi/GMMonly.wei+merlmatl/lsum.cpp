#include "mex.h"
#include <math.h>
#include <iostream.h>

#ifdef __WIN32
#define finite(a) isfinite(a)
#endif

void mexFunction( int nlhs, mxArray **plhs, int nrhs, const mxArray **prhs)
{
	// Check
	if( nrhs != 1)
		mexErrMsgTxt( "Wrong number of inputs!");

	// Remember these
	int m = mxGetM( prhs[0]), n = mxGetN( prhs[0]);
	double *x = mxGetPr( prhs[0]);

	// Allocate output
	plhs[0] = mxCreateDoubleMatrix( m, 1, mxREAL);
	double *y = mxGetPr( plhs[0]);
	
	// Do your thing ...
	for( int r = 0 ; r != m ; r++)
		y[r] = *x++;

	for( int c = 1 ; c != n ; c++){
		for( int r = 0 ; r != m ; r++){
			if( finite( *x)){
				if( finite( y[r])){
					double z = *x + log( 1. + exp( y[r] - *x));
					if( !finite( z))
						z = y[r] + log( 1. + exp( *x - y[r]));
					if( finite( z))
						y[r] = z;
				}else
					y[r] = *x;
			}
			*x++;
		}
	}
}
