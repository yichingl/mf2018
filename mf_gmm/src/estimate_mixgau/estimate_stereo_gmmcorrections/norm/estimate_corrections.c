#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "header.h"

void estimate_gmm_stereo_corrections(float **x,   /* observation vectors */
                                     float **y,   /* Stereo counterparts */
                                     int N,       /* No. of vectors */
                                     int Ndim,    /* Dimensionality of data */
                                     int K,       /* No. of Gaussians */
                                     float **mean,/* Means of distribution */
                                     float **var, /* Variances */
                                     float *c,    /* Mixture weights */
                                     int   R,     /* "Lookahead frames" */
                                     int   L,     /* "Lookbehind frames" */
                                     float ***meany, /* Means of corrections */
                                     float *count  /* Counts to normalize */
                                    )
{
    float    **hafinvvar, *Tau, *corprod, Const;

    int     i,j,k;
    float logprob;

   /*
    * The Constant term that occurs in the computation of the gaussians.
    * As it turns out, we never use it. This one is a dummy :-)
    */
    Const = pow(2*M_PI,(Ndim/2));

   /*
    * some initializations
    */

   /*
    * Allocate spaces for the local parameter arrays
    */
    if ((hafinvvar = (float **) alloc2d(K,Ndim,sizeof(float))) == NULL)
        QUIT(("Unable to allocate for Newvar array\n"));

    if ((Tau = (float *)malloc(K*sizeof(float))) == NULL)
        QUIT(("Unable to allocate %d floats for Tau array \n",K));
    
    if ((corprod = (float *)malloc(K*sizeof(float))) == NULL)
        QUIT(("Unable to allocate %d floats for corprod array \n",K));


    /* 
     * Compute the ratio of the prior probabilities of the modulus of the
     * variance. We do this operation in the log domain to keep the
     * dynamic range in check
     */

    for (k=0;k<K;++k) {
        corprod[k] = log((double)c[k]);
        for (j=0;j<Ndim;++j){
             corprod[k] -= 0.5*log((double)var[k][j]);
             if (var[k][j] > 0) hafinvvar[k][j] = 1.0/(2.0*var[k][j]);
             else hafinvvar[k][j] = 1e+20;
        }
    }

   /*
    * Estimate means and variances and priors while computing overall
    * Likelihood
    * Note: variance estimated as Sum((x-mean)*(x-mean))/count
    * rather than Sum(x*x)/count - mean*mean
    * because the latter is an unstable formula and tends to give -ve
    * variances due to numerical errors
    */
    logprob = 0;
    for (i=L;i<N-R;++i) { // Just to simplify, start at L and and at end-R
        int ctx;
        logprob += Expectation(Tau,x[i],c,mean,hafinvvar,corprod,K,Ndim);
        for (k=0; k<K ; ++k) count[k] += Tau[k];
        for (ctx = -L; ctx <= R; ctx++) {
            float correction[8192];
            for (j=0;j<Ndim;++j) correction[j] = y[i+ctx][j] - x[i][j];
            for (k=0; k<K ; ++k) {
                if (Tau[k] > 0){
                    for (j=0;j<Ndim;++j)
                        meany[ctx+L][k][j] += Tau[k] * correction[j];
                }
            }
        }
    }
    /* Free local arrays */
    free(Tau);
    free(corprod);
    return;
}



/*----------------------------------------------------------------------------
  Routine to compute the a-posteriori probabilities of each of the modes in
  the Mixture distribution defined by c,mean and var for the the vector
  defined by x[]. Returns the posteriors in Tau[]
------------------------------------------------------------------------------*/
float    Expectation(float *Tau,   /* Mode probablity of observation vector x */
                     float *x,     /* Observation vector x */ 
                     float *c,     /* Mixing proportions of modes */
                     float **mean, /* Means of the modes */
                     float **hafinvvar,  /* 1.0/(2.0*variance) of the modes */
                     float *corprod, /* 1/sqrt(2PI^Ndim*|variance|) */
                     int K,        /* Number of modes */
                     int Ndim)     /* Dimensionality of data */
{
    float    Probability,
             max,sum,invsum;

    int      i;

    Tau[0] = max = gauss(x,mean[0],hafinvvar[0],corprod[0],Ndim);
    for (i=1;i<K;++i) {
        /*
         * Note: prior c[i] has already been factored into corprod
         */
        Tau[i] = gauss(x,mean[i],hafinvvar[i],corprod[i],Ndim);
        if (Tau[i] > max) max = Tau[i];
    }

    /* convert Tau from log to linear */
    sum = 0.0;
    for (i=0;i<K;++i) {
        if (Tau[i] > max - 22) { /* only consider probs within 10^10 of max */
            Tau[i] = exp(Tau[i] - max + 10);
            sum += Tau[i];
        }
        else Tau[i] = 0;
    }
    invsum = 1.0/sum;
    for (i=0;i<K;++i)
        Tau[i] *= invsum;

    Probability = max + log(sum) - 10;
    return(Probability);
}


/*------------------------------------------------------------------------------
	Function defining the gaussian density function given mean and variance
------------------------------------------------------------------------------*/

float gauss(float *x,
            float *mean,
            float *hafinvvar,
            float corprod,
            int Ndim)
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
