/*------------------------------------------------------------------------------
	Function defining the gaussian density function given mean and variance
------------------------------------------------------------------------------*/
#include <math.h>

float gauss(x,mean,hafinvvar,corprod,Ndim)
float	*x,*mean,*hafinvvar,corprod;
int	Ndim;
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
