#include "EMheader.h"

/* First, a function to read distributions (that were stored using 
   store distributions */
int  read_distribution(char *file, int *Nmodes, int *Ndim,
                       REAL **c, REAL ***mean, REAL ***var)
{
    int     i, j, Nm, Nd;
    REAL    *lc, **lmean, **lvar;
    FILE    *fopen(),*inputfile;

    if ((inputfile = fopen(file,"r")) == NULL) return(-1);

    fscanf(inputfile,"%d %d",&Nm, &Nd);
    *Nmodes = Nm; *Ndim = Nd;

    *c = lc = (REAL *)calloc(Nm, sizeof(REAL));
    *mean = lmean = (REAL **) alloc2d(Nm, Nd, sizeof(REAL));
    *var = lvar = (REAL **) alloc2d(Nm, Nd, sizeof(REAL));

    for (i=0;i<Nm;++i) {
        fscanf(inputfile,"%f",&lc[i]);
        for (j=0;j<Nd;++j) fscanf(inputfile,"%f", &lmean[i][j]);
        for (j=0;j<Nd;++j) fscanf(inputfile,"%f", &lvar[i][j]);
    }
    fclose(inputfile);
    return(0);
}



/* First, a function to store distributions */
int  store_distribution(char *file, int Nmodes, int Ndim,
                        REAL *c, REAL **mean, REAL **var)
{
    int     i,j;
    FILE    *fopen(),*outputfile;

    if ((outputfile = fopen(file,"w")) == NULL) return(-1);

    fprintf(outputfile,"%d %d\n\n",Nmodes, Ndim);
    
    for (i=0;i<Nmodes;++i) {
        fprintf(outputfile,"%e\n",c[i]);
        for (j=0;j<Ndim;++j)
            fprintf(outputfile,"%e ",mean[i][j]);
        fprintf(outputfile,"\n");
        for (j=0;j<Ndim;++j)
            fprintf(outputfile,"%e ",var[i][j]);
        fprintf(outputfile,"\n\n");
    }
    fclose(outputfile);
    return(0);
}


/*------------------------------------------------------------------------------
	Function defining the gaussian density function given mean and variance
------------------------------------------------------------------------------*/
#include <math.h>

REAL gauss(REAL *x, REAL *mean, REAL *halfinvvar, REAL corprod, int Ndim)
{
    int	        j;
    REAL	density,t;

    density = corprod;
    for (j=0;j<Ndim;++j) {
        t = *x++ - *mean++;
	density -= t*t*(*halfinvvar++);
    }
    return(density);
}

/*----------------------------------------------------------------------------
  Routine to compute the a-posteriori probabilities of each of the modes in
  the Mixture distribution defined by c,mean and var for the the vector
  defined by x[]. Returns the posteriors in Tau[]
------------------------------------------------------------------------------*/
REAL     Expectation(REAL *Tau,         /* A posteriori Gaussian probs. */
                     REAL *x,           /* Observation vector */
                     REAL *c,           /* Mixture weights */
                     REAL **mean,       /* Means of Gaussians */
                     REAL **halfinvvar, /* 1.0/(2.0*variance) of Gaussians */
                     REAL *corprod,     /* log(mixw/sqrt(2PI^Ndim*det(cov))) */
                      int K,            /* No. of Gaussians in mixture */
                      int Ndim          /* Data dimension */
                    )
{
    REAL    Probability,
             max,sum,invsum;

    int      i;

    Tau[0] = max = gauss(x,mean[0],halfinvvar[0],corprod[0],Ndim);
    for (i=1;i<K;++i)
    {
        /*
         * Note: prior c[i] has already been factored into corprod
         */
        Tau[i] = gauss(x,mean[i],halfinvvar[i],corprod[i],Ndim);
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


/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*
    Program to compute a Multimodal gaussian distribution for a given
    data set. Estimates DIAGNOAL covariances ONLY (i.e all modes are assumed
    to have diagonal covariance matrices).
    Does NOT compute optimal number of modes; has to be given the number of
    modes required.
    The routine runs a variable number of iterations of EM each time it 
    is called. The number of iterations is passed as a variable

    The parameters passed are 
	x[][] - the data set (REAL)
        N - the number of vectors in the data set (int)
	Ndim - the dimensionality of the dataset (int)
	K - the number of modes required in the computed distribution (int)
	mean[][] - Some initial values for the means of the gaussians(REAL)
	var[][] - Initial values for the Variances of the modes(REAL)
	c[] - Initial values for the a-priori probabilities(REAL)
        tempfile - File to store partially converged distributions to
        numiters - Number of EM iterations to run

      Before passing to this routing mean, var and c have to be initialized.
      The mean, var and c values can be initialized by doing a quick VQ on
      the data and using the means,variances and population ratios of each
      of the bins for initial values.

    The returned values are
	c[] - the a-priori probabilities of each of the modes, i.e the
	      mixing proportions (REAL)
        mean[][] - the means of the modes (REAL)
	var[][] - the variances of the modes (REAL)

    Note - This program requires that the varibles x[][], mean[][] and var[][]
           have been allocated using the alloc2d() routine!
*vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "header.h"

REAL expectation_maximization(REAL **x,       /* The observation vectors */
                               int N,         /* Number of obs. vectors */
                               int Ndim,      /* Dimensions in obs. vectors */
                               int K,         /* Number of Gaussians needed */
                               REAL **mean,    /* The Gaussian means */
                               REAL **var,     /* The Gaussian covarinances */
                               REAL *c,        /* Mixture weights */
                               char *tempfile, /* Temporary file name */
                                int numiters   /* No. EM iterations needed */
                             )
{
    REAL    **halfinvvar,
            *Sharedcov,
            **Newmean,
            *Newc,
            *Tau,
            *corprod,
            Const,
            SumNewc,
            Prevlogprob,LogProb,Threshold,Improvement;

    int     i,j,k,iter=0;

   /*
    * The Constant term that occurs in the computation of the gaussians.
    * As it turns out, we never use it. This one is a dummy :-)
    */
    Const = pow(2*PI,(Ndim/2));

   /*
    * some initializations
    */
    Threshold = 1.0e-06;
    Improvement = 100;

   /*
    * Allocate spaces for the local parameter arrays
    */
    if ((Newmean = (REAL **)alloc2d(K,Ndim,sizeof(REAL))) == NULL)
        QUIT(("Unable to allocate Newmean array\n"));

    if ((Sharedcov = (REAL *) calloc(Ndim,sizeof(REAL))) == NULL)
        QUIT(("Unable to allocate for Newvar array\n"));

    if ((halfinvvar = (REAL **) alloc2d(K,Ndim,sizeof(REAL))) == NULL)
        QUIT(("Unable to allocate for Newvar array\n"));

    if ((Newc = (REAL *)malloc(K*sizeof(REAL))) == NULL)
        QUIT(("Unable to allocate %d REALs for Newc array \n",K));

    if ((Tau = (REAL *)malloc(K*sizeof(REAL))) == NULL)
        QUIT(("Unable to allocate %d REALs for Tau array \n",K));
    
    if ((corprod = (REAL *)malloc(K*sizeof(REAL))) == NULL)
        QUIT(("Unable to allocate %d REALs for corprod array \n",K));

    /*
     * Initialize all New values to 0
     * Note the array position computation for Newmean, as it is a 1-D array
     */
    for (k=0; k<K ; ++k)
    {
        Newc[k] = 0;
        for (j=0;j<Ndim;++j)
	{
            Newmean[k][j] = 0;
	}
    }
    for (j=0;j<Ndim;++j) Sharedcov[j] = 0;

    /* 
     * Compute the ratio of the prior probabilities of the modulus of the
     * variance. We do this operation in the log domain to keep the
     * dynamic range in check
     */

    for (k=0;k<K;++k)
    {
        corprod[k] = log((float)c[k]);
        for (j=0;j<Ndim;++j){
             corprod[k] -= 0.5*log((float)var[k][j]);
             if (var[k][j] > 0) halfinvvar[k][j] = 1.0/(2.0*var[k][j]);
             else halfinvvar[k][j] = 1e+20;
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
        Prevlogprob += Expectation(Tau,x[i],c,mean,halfinvvar,corprod,K,Ndim);
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
        Expectation(Tau,x[i],c,mean,halfinvvar,corprod,K,Ndim);
        for (k=0; k<K ; ++k)
        {
            if (Tau[k] > 0) {
                for (j=0;j<Ndim;++j)
		    Sharedcov[j] += Tau[k] * (x[i][j] - Newmean[k][j]) 
                                       * (x[i][j] - Newmean[k][j]);
            }
        }
    }
    printf("EM : Initial log probablity = %f \n",Prevlogprob);

    while ((Improvement > Threshold)&&(iter < numiters))
    {
        /* 
         * We use SumNewc instead of N as in the formula because the
         * Newc's may not sum to N, because of accuracy errors of the 
         * computer for large N.
         */
        SumNewc = 0;
        for (k=0;k<K;++k)
            SumNewc += Newc[k];

        for (k=0; k<K ; ++k)
        {
            for (j=0;j<Ndim;++j)
	    {
                mean[k][j] = Newmean[k][j];
		var[k][j] = Sharedcov[j]/SumNewc;
            }
            c[k] = Newc[k]/SumNewc;
        }

        if (store_distribution(tempfile,K,Ndim,c,mean,var) == -1) {
            fprintf(stderr,"Error storing distribution in %s\n",tempfile);
            fflush(stderr);
        }

        /*
         * Initialize all New values to 0
         */
        for (k=0; k<K ; ++k)
        {
            Newc[k] = 0;
            for (j=0;j<Ndim;++j)
	    {
                Newmean[k][j] = 0;
            }
        }
        for (j=0;j<Ndim;++j) Sharedcov[j] = 0;

        /* 
         * Compute the ratio of the prior probabilities of the modulus of the
         * variance. We do this operation in the log domain to keep the
         * dynamic range in check
         */

        for (k=0;k<K;++k)
        {
            corprod[k] = log((float)c[k]);
            for (j=0;j<Ndim;++j) {
                corprod[k] -= 0.5*log((float)var[k][j]);
                if (var[k][j] > 0) halfinvvar[k][j] = 1.0/(2.0*var[k][j]);
                else halfinvvar[k][j] = 1e+20;
            }
        }

        LogProb = 0;
        for (i=0;i<N;++i)
        {
            LogProb += Expectation(Tau,x[i],c,mean,halfinvvar,corprod,K,Ndim);
            for (k=0; k<K ; ++k)
            {
                if (Tau[k] > 0) {
                    Newc[k] += Tau[k];
                    for (j=0;j<Ndim;++j) {
                        Newmean[k][j] += Tau[k] * x[i][j];
                        Sharedcov[j] += Tau[k]*(x[i][j] - mean[k][j])*
                                               (x[i][j] - mean[k][j]);
                    }
                }
            }
        }
        for (k=0; k<K ; ++k)
            for (j=0;j<Ndim;++j)
                Newmean[k][j] /= Newc[k];

        Improvement = (LogProb - Prevlogprob)/LogProb;
        if (LogProb < 0) Improvement = -Improvement;
        ++iter;
        printf("EM : Log Prob = %f, improvement = %f after %d iterations\n",
                     LogProb,Improvement,iter);  fflush(stdout);
        Prevlogprob = LogProb;
    }
   /*
    * Free local arrays
    */
    free(Tau);
    free(corprod);
    free(Sharedcov);
    free2d((void**)Newmean);
    free(Newc);
    return LogProb;
}


void estimate_mixture_gaussian(REAL **x,       /* The observation vectors */
                                int N,         /* Number of obs. vectors */
                                int Ndim,      /* Dimensions in obs. vectors */
                                int K,         /* Number of Gaussians needed */
                               REAL **mean,    /* The Gaussian means */
                               REAL **var,     /* The Gaussian covarinances */
                               REAL *c,        /* Mixture weights */
                               char *tempfile, /* Temporary file name */
                                int inputNgau  /* If we're continuing from a
                                                  previous distribution. If
                                                  we should begin from scratch,
                                                  set this to 0 */
                              )
{
    int d, i, j, Ngauss, NumGauToSplit;
    float lkd;
    char tdistfile[1024];

    if (inputNgau == 0) {
        int t;

        // First set mean[0] and variance[0] to global mean and variance;
        for (t=0; t<N; t++) {
            for (d=0; d<Ndim; d++) {
                mean[0][d] += x[t][d];
            }
        }
        for (d=0; d<Ndim; d++) mean[0][d] /= (REAL) N;
        for (t=0; t<N; t++) {
            for (d=0; d<Ndim; d++) {
                var[0][d] += (x[t][d] - mean[0][d]) * (x[t][d] - mean[0][d]);
            }
        }
        for (d=0; d<Ndim; d++) var[0][d] /= (REAL) N;
        c[0] = 1.0;
        Ngauss = 1;
    }
    else 
        Ngauss = inputNgau;


    while (Ngauss < K) {
        // Find out how many Gaussians to split
        NumGauToSplit = 2*Ngauss < K ? Ngauss : K - Ngauss;

        // Split the first NumGauToSplit Gaussians
        for (i=0, j=Ngauss; i<NumGauToSplit; i++, j++) {
            c[j] = c[i] * 0.5; c[i] *= 0.5;
            for (d = 0; d < Ndim; d++) {
                mean[j][d] = mean[i][d] + 0.01*(REAL)sqrt((float)var[i][d]);
                mean[i][d] = mean[i][d] - 0.01*(REAL)sqrt((float)var[i][d]);
                var[j][d] = var[i][d];
            }
        }
        Ngauss += NumGauToSplit;
        sprintf(tdistfile,"%s.%d",tempfile,Ngauss);
        lkd = expectation_maximization(x,N,Ndim,Ngauss,mean,var,c,
                                       tdistfile,NUMITERS);
        if (store_distribution(tdistfile,Ngauss,Ndim,c,mean,var) == -1) {
            fprintf(stderr,"Error storing distribution in %s\n",tdistfile);
            fflush(stderr);
        }
        else {
            fprintf(stdout,"%d Gaussians trained,likelihood = %f\n",Ngauss,lkd);
            fprintf(stdout,"Stored distribution in %s\n",tdistfile);
            fflush(stdout);
        }
    }
    return;
}
