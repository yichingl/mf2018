#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "header.h"


/*----------------------------------------------------------------------------
  Routine to compute the a-posteriori probabilities of each of the modes in
  the Mixture distribution defined by c,mean and var for the the vector
  defined by x[]. Returns the posteriors in Tau[]
------------------------------------------------------------------------------*/
float    Expectation(float *Tau, float *x, float *c, float **mean,
                     float **hafinvvar, float *corprod, int Ngau, int Ndim)
{
    float    Probability,
             max,sum,invsum;

    int      i;

    Tau[0] = max = gauss(x,mean[0],hafinvvar[0],corprod[0],Ndim);
    for (i=1;i<Ngau;++i)
    {
        /*
         * Note: prior c[i] has already been factored into corprod
         */
        Tau[i] = gauss(x,mean[i],hafinvvar[i],corprod[i],Ndim);
        if (Tau[i] > max) max = Tau[i];
    }

    /* convert Tau from log to linear */
    sum = 0.0;
    for (i=0;i<Ngau;++i)
    {
        if (Tau[i] > max - 22) { /* only consider probs within 10^10 of max */
            Tau[i] = exp(Tau[i] - max + 10);
            sum += Tau[i];
        }
        else Tau[i] = 0;
    }
    invsum = 1.0/sum;
    for (i=0;i<Ngau;++i)
        Tau[i] *= invsum;

    Probability = max + log(sum) - 10;
    return(Probability);
}

void    mgau_em(float **x, int Nvec, int Ndim, int Ngau,
                float **mean, float **var, float *c, int numiters)
{
    float    **Newvar, **hafinvvar, **Newmean, *Newc, *Tau, *corprod;
    float    Const, SumNewc, Prevlogprob, LogProb, Threshold, Improvement;
    const float varthresh = 1e-06;
    int     i,j,k,iter=0;

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
    * Initialize everything
    */
    initialize_mgau(x, Nvec, Ndim, mean, var, c, Ngau);


   /*
    * Allocate spaces for the local parameter arrays
    */
    if ((Newmean = (float **)alloc2d(Ngau,Ndim,sizeof(float))) == NULL)
        QUIT(("Unable to allocate Newmean array\n"));

    if ((Newvar = (float **) alloc2d(Ngau,Ndim,sizeof(float))) == NULL)
        QUIT(("Unable to allocate for Newvar array\n"));

    if ((hafinvvar = (float **) alloc2d(Ngau,Ndim,sizeof(float))) == NULL)
        QUIT(("Unable to allocate for Newvar array\n"));

    if ((Newc = (float *)malloc(Ngau*sizeof(float))) == NULL)
        QUIT(("Unable to allocate %d floats for Newc array \n",Ngau));

    if ((Tau = (float *)malloc(Ngau*sizeof(float))) == NULL)
        QUIT(("Unable to allocate %d floats for Tau array \n",Ngau));
    
    if ((corprod = (float *)malloc(Ngau*sizeof(float))) == NULL)
        QUIT(("Unable to allocate %d floats for corprod array \n",Ngau));

    /*
     * Initialize all New values to 0
     * Note the array position computation for Newmean, as it is a 1-D array
     */
    for (k=0; k<Ngau ; ++k)
    {
        Newc[k] = 0;
        for (j=0;j<Ndim;++j)
	{
            Newmean[k][j] = 0;
	    Newvar[k][j] = 0;
	}
    }

    /* 
     * Compute the ratio of the prior probabilities of the modulus of the
     * variance. We do this operation in the log domain to keep the
     * dynamic range in check
     */

    for (k=0;k<Ngau;++k)
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
    for (i=0;i<Nvec;++i)
    {
        Prevlogprob += Expectation(Tau,x[i],c,mean,hafinvvar,corprod,Ngau,Ndim);
        for (k=0; k<Ngau ; ++k)
        {
            if (Tau[k] > 0){
                Newc[k] += Tau[k];
                for (j=0;j<Ndim;++j)
                    Newmean[k][j] += Tau[k] * x[i][j];
            }
        }
    }
    for (k=0; k<Ngau ; ++k)
        for (j=0;j<Ndim;++j)
	    Newmean[k][j] /= Newc[k];
    for (i=0;i<Nvec;++i)
    {
        Expectation(Tau,x[i],c,mean,hafinvvar,corprod,Ngau,Ndim);
        for (k=0; k<Ngau ; ++k)
        {
            if (Tau[k] > 0) {
                for (j=0;j<Ndim;++j)
		    Newvar[k][j] += Tau[k] * (x[i][j] - Newmean[k][j]) 
                                       * (x[i][j] - Newmean[k][j]);
            }
        }
    }
    printf("EM : Initial log probablity = %e \n",Prevlogprob);

    while ((Improvement > Threshold)&&(iter < numiters))
    {
        /* 
         * We use SumNewc instead of Nvec as in the formula because the
         * Newc's may not sum to Nvec, because of accuracy errors of the 
         * computer for large Nvec.
         */
        SumNewc = 0;
        for (k=0;k<Ngau;++k)
            SumNewc += Newc[k];

        for (k=0; k<Ngau ; ++k)
        {
            for (j=0;j<Ndim;++j)
	    {
                mean[k][j] = Newmean[k][j];
		var[k][j] = Newvar[k][j]/Newc[k];
                if (var[k][j] < varthresh) var[k][j] = varthresh;
            }
            c[k] = Newc[k]/SumNewc;
        }

        /*
         * Initialize all New values to 0
         */
        for (k=0; k<Ngau ; ++k)
        {
            Newc[k] = 0;
            for (j=0;j<Ndim;++j)
	    {
                Newmean[k][j] = 0;
		Newvar[k][j] = 0;
            }
        }

        /* 
         * Compute the ratio of the prior probabilities of the modulus of the
         * variance. We do this operation in the log domain to keep the
         * dynamic range in check
         */

        for (k=0;k<Ngau;++k)
        {
            corprod[k] = log((double)c[k]);
            for (j=0;j<Ndim;++j) {
                corprod[k] -= 0.5*log((double)var[k][j]);
                if (var[k][j] > 0) hafinvvar[k][j] = 1.0/(2.0*var[k][j]);
                else hafinvvar[k][j] = 1e+20;
            }
        }

        LogProb = 0;
        for (i=0;i<Nvec;++i)
        {
            LogProb += Expectation(Tau,x[i],c,mean,hafinvvar,corprod,Ngau,Ndim);
            for (k=0; k<Ngau ; ++k)
            {
                if (Tau[k] > 0) {
                    Newc[k] += Tau[k];
                    for (j=0;j<Ndim;++j) {
                        Newmean[k][j] += Tau[k] * x[i][j];
                        Newvar[k][j] += Tau[k]*(x[i][j] - mean[k][j])*
                                               (x[i][j] - mean[k][j]);
                    }
                }
            }
        }
        for (k=0; k<Ngau ; ++k)
            for (j=0;j<Ndim;++j)
                Newmean[k][j] /= Newc[k];

/*
        for (i=0;i<Nvec;++i)
        {
            Expectation(Tau,x[i],c,mean,hafinvvar,corprod,Ngau,Ndim);
            for (k=0; k<Ngau ; ++k)
            {
                for (j=0;j<Ndim;++j)
		    Newvar[k][j] += Tau[k] * (x[i][j] - Newmean[k][j])
		                           * (x[i][j] - Newmean[k][j]);
            }
        }
*/

        Improvement = (LogProb - Prevlogprob)/LogProb;
        if (LogProb < 0) Improvement = -Improvement;
        ++iter;
        printf("EM : Log Prob = %e, improvement = %e after %d iterations\n",
                     LogProb,Improvement,iter);  fflush(stdout);
        Prevlogprob = LogProb;
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

