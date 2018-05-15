/*------------------------------------------------------------------------------
	Function defining the gaussian density function given mean and variance
------------------------------------------------------------------------------*/
#include <math.h>

float gauss(float *x, float *mean, float *hafinvvar, float corprod, int Ndim)
{
    int	        j;
    float	density,t;

    density = corprod;
    for (j=0;j<Ndim;++j) {
        t = *x++ - *mean++;
	density -= t*t*(*hafinvvar++);
    }
    return(density);
}
