#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "header.h"

void    estimate_multi_modals(float **x,   /* The observation vectors */
                              float **y,   /* Stereo counterpart of x */
                              int N,       /* No. of vectors */
                              int Ndim,    /* Dimensionality of data */
                              int K,       /* No. of Gaussians to estimate */
                              float **mean,/* Means of distribution */
                              float **var, /* Variances of distribution */
                              float *c,    /* Mixture weights */
                              int   R,     /* "Lookahead prediction frames" */
                              int   L,     /* "Lookbehind frames that predict */
                              float ***covyx, /* Full covar of x */
                              float ****covyx, /* Covar of y & x at all lags */
                              float ***meany, /* Means of y at each lag */
                              float **count  /* Count terms to normalize */
                              )
{
    float    **hafinvvar,
             *Tau,
             *Newc,
             *corprod,
             Const,
             SumNewc;

    const float varthresh = 1e-06;

    int     i,j,k,l;

    float   Expectation();
            
    FILE    *temp,*fopen();

   /*
    * The Constant term that occurs in the computation of the gaussians.
    * As it turns out, we never use it. This one is a dummy :-)
    */
    Const = pow(2*PI,(Ndim/2));

   /*
    * some initializations
    */
    Threshold = 1.0e-08;
    Improvement = 100;

   /*
    * Allocate spaces for the local parameter arrays
    */
    if ((hafinvvar = (float **) alloc2d(K,Ndim,sizeof(float))) == NULL)
        QUIT(("Unable to allocate for Newvar array\n"));

    if ((Newc = (float *)malloc(K*sizeof(float))) == NULL)
        QUIT(("Unable to allocate %d floats for Newc array \n",K));

    if ((Tau = (float *)malloc(K*sizeof(float))) == NULL)
        QUIT(("Unable to allocate %d floats for Tau array \n",K));
    
    if ((corprod = (float *)malloc(K*sizeof(float))) == NULL)
        QUIT(("Unable to allocate %d floats for corprod array \n",K));


    /* 
     * Compute the ratio of the prior probabilities of the modulus of the
     * variance. We do this operation in the log domain to keep the
     * dynamic range in check
     */

    for (k=0;k<K;++k)
    {
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
    Prevlogprob = 0;
    for (i=0;i<N;++i)
    {
        Prevlogprob += Expectation(Tau,x[i],c,mean,hafinvvar,corprod,K,Ndim);
        for (k=0; k<K ; ++k)
        {
            if (Tau[k] > 0){
                Newc[k] += Tau[k];
                for (j=0;j<Ndim;++j)
                    Newmean[k][j] += Tau[k] * x[i][j];
            }
        }
    }
    for (k=0; k<K ; ++k)
        for (j=0;j<Ndim;++j)
	    Newmean[k][j] /= Newc[k];
    for (i=0;i<N;++i)
    {
        Expectation(Tau,x[i],c,mean,hafinvvar,corprod,K,Ndim);
        for (k=0; k<K ; ++k)
        {
            if (Tau[k] > 0) {
                for (j=0;j<Ndim;++j)
		    Newvar[k][j] += Tau[k] * (x[i][j] - Newmean[k][j]) 
                                       * (x[i][j] - Newmean[k][j]);
            }
        }
    }

   /*
    * Free local arrays
    */
    free(Tau);
    free(corprod);
    free2d((void**)Newvar);
    free2d((void**)Newmean);
    free(Newc);
    return;
}



/*----------------------------------------------------------------------------
  Routine to compute the a-posteriori probabilities of each of the modes in
  the Mixture distribution defined by c,mean and var for the the vector
  defined by x[]. Returns the posteriors in Tau[]
------------------------------------------------------------------------------*/
float    Expectation(Tau,x,c,mean,hafinvvar,corprod,K,Ndim)
    float   Tau[],           /* Mode probablity of observation vector x */
            x[],             /* Observation vector x */
            c[],             /* Mixing proportions of modes */
            **mean,          /* Means of the modes */
            **hafinvvar,     /* 1.0/(2.0*variance) of the modes */
            corprod[];       /* Precomputed 1/sqrt(2PI**Ndim * mod(variance)) */

    int     K,               /* Number of modes */
            Ndim;            /* Dimensionality of observation/distribution */
{
    float    Probability,
             max,sum,invsum;

    int      i;

    float   gauss();


    Tau[0] = max = gauss(x,mean[0],hafinvvar[0],corprod[0],Ndim);
    for (i=1;i<K;++i)
    {
        /*
         * Note: prior c[i] has already been factored into corprod
         */
        Tau[i] = gauss(x,mean[i],hafinvvar[i],corprod[i],Ndim);
        if (Tau[i] > max) max = Tau[i];
    }

    /* convert Tau from log to linear */
    sum = 0.0;
    for (i=0;i<K;++i)
    {
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
