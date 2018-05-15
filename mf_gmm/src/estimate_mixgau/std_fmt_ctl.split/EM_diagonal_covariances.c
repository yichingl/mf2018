#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define REAL float
#define NUMITERS 30

#define QUIT(X) {printf("%s(%d):",__FILE__,__LINE__);\
                 printf X;fflush(stdout);exit(0);}

int read_cep(char *file, float***cep, int *numframes, int cepsize);


char          **alloc2d(int dim1,	/* "x" dimension */ 
			int dim2,	/* "y" dimension */
			int size	/* number of bytes each entry takes */
			)
{
	int             i;		/* loop control variable */
	unsigned        nelem;		/* total number of elements */

	char	       *p,		/* pointer to matrix memory */
		      **pp;		/* pointer to matrix mem table */

	/*
	 * Compute total number of elements needed for the two-dimensional
	 * matrix
	 */
	nelem = (unsigned) dim1 * dim2;

	/*
	 * Allocate the memory needed for the matrix
	 */
	p = calloc(nelem, (unsigned) size);

	/*
	 * If the allocation were not successful, return a NULL pointer
	 */
	if (p == NULL) return (NULL);

	/*
	 * Now allocate a table for referencing the matrix memory
	 */
	pp = (char **) calloc((unsigned) dim1, (unsigned) sizeof(char *));

	/*
	 * If the allocation were not successful, return a NULL pointer
	 * and free the previously allocated memory
	 */
	if (pp == NULL) {
		free(p);
		return (NULL);
	}

	/*
	 * Fill the table with locations to where each row begins
	 */
	for (i = 0; i < dim1; i++)
		pp[i] = p + (i * dim2 * size);

	return (pp);
}

void free2d(void **p)
{
   if (p!=NULL) {
      free(p[0]);
      free(p);
   }
   return;
}


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

    Coded by Bhiksha Raj

    Incoroporated excellent changes by PJM - Begin working in the log domain
    to eliminate corprod scaling problems. Dec 95
*vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv*/

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
    REAL    **Newvar,**halfinvvar,
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
    Const = pow(2*M_PI,(Ndim/2));

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

    if ((Newvar = (REAL **) alloc2d(K,Ndim,sizeof(REAL))) == NULL)
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
	    Newvar[k][j] = 0;
	}
    }

    /* 
     * Compute the ratio of the prior probabilities of the modulus of the
     * variance. We do this operation in the log domain to keep the
     * dynamic range in check
     */

    for (k=0;k<K;++k)
    {
        corprod[k] = log((double)c[k]);
        for (j=0;j<Ndim;++j){
             if (var[k][j] < 1e-10) var[k][j] = 1e-10; // Floor
             corprod[k] -= 0.5*log((double)var[k][j]);
             halfinvvar[k][j] = 1.0/(2.0*var[k][j]);
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
		    Newvar[k][j] += Tau[k] * (x[i][j] - Newmean[k][j]) 
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

        for (k=0; k<K ; ++k)
        {
            for (j=0;j<Ndim;++j)
	    {
                mean[k][j] = Newmean[k][j];
		var[k][j] = Newvar[k][j]/Newc[k];
                if (var[k][j] < 1e-10) var[k][j] = 1e-10; // Floor
            }
        }
        SumNewc = 0;
        for (k=0;k<K;++k) {
            Newc[k] += 1e-06/(float)K;  // Smooth slightly
            SumNewc += Newc[k];
        }
        for (k=0;k<K;++k) c[k] = Newc[k]/SumNewc;

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
		Newvar[k][j] = 0;
            }
        }

        /* 
         * Compute the ratio of the prior probabilities of the modulus of the
         * variance. We do this operation in the log domain to keep the
         * dynamic range in check
         */

        for (k=0;k<K;++k)
        {
            corprod[k] = log((double)c[k]);
            for (j=0;j<Ndim;++j) {
                corprod[k] -= 0.5*log((double)var[k][j]);
                halfinvvar[k][j] = 1.0/(2.0*var[k][j]);
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
                        Newvar[k][j] += Tau[k]*(x[i][j] - mean[k][j])*
                                               (x[i][j] - mean[k][j]);
                    }
                }
            }
        }
        for (k=0; k<K ; ++k)
            for (j=0;j<Ndim;++j)
                Newmean[k][j] /= Newc[k];

/*
        for (i=0;i<N;++i)
        {
            Expectation(Tau,x[i],c,mean,halfinvvar,corprod,K,Ndim);
            for (k=0; k<K ; ++k)
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
        printf("EM : Log Prob = %f, improvement = %f after %d iterations\n",
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
    return LogProb;
}

int findheaviest(REAL *c, int j)
{
    int i, m=0;
    for (i=1; i<j; i++) if (c[i] > c[m]) m = i;
    return m;
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
    int d, i, j, k, Ngauss, NumGauToSplit;
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

        // Split NumGauToSplit heaviest Gaussians
        for (i=0, j=Ngauss; i<NumGauToSplit; i++, j++) {
            k = findheaviest(c,j);
            c[j] = c[k] * 0.5; c[k] *= 0.5;
            for (d = 0; d < Ndim; d++) {
                float jitter = 0.03*drand48()*(REAL)sqrt((double)var[k][d]);
                mean[j][d] = mean[k][d] + jitter;
                mean[k][d] = mean[k][d] - jitter;
                var[j][d] = var[k][d];
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



/******************************************************************************
*    Main routine for computing the distribution
******************************************************************************/


#define DESCRIPTION "\nThis program computes a mixture gaussian distribution for a set of speech\ncepstral files. The dimensionality of the cepstral files is assumed to be 13\n"

#define USAGE "USAGE:\n%s -e <filename extension (def none)> \\\n\t-c <control file> \\\n\t-i <indir (def none)> \\\n\t-o <output file> \\\n\t-m <Total no of modes in distribution>\\\n\t-d <dimensionality of the data (default 13)> \\\n\t-n <For CMN> \\\n\t-l <To include delta features> \\\n\t-r <fraction of lowest-0th-term vectors to discard (def 0)>\n"

int computefeat(float ***buff, int CMS, int DELTA, float FILTERLOWEGY,
                int Ndim, int *length)
{
    int i, j;
    float **lbuff = *buff;

    if (FILTERLOWEGY > 0) { /* Remove lowest egy vectors within 5% of min egy */
        float min = 1e+30;
        float max = -1e+30;
        float c0thresh;
        int k;
        float **ltmp = *buff;

        for (i=0; i<*length; i++) {
            if (ltmp[i][0] < min) min = ltmp[i][0];
            if (ltmp[i][0] > max) max = ltmp[i][0];
        }
        c0thresh = min + (max - min)*FILTERLOWEGY;
        for (i=0,k=0; i<*length; i++) {
            if (ltmp[i][0] < c0thresh) continue;
            for (j=0; j<Ndim; j++) ltmp[k][j] = ltmp[i][j];
            k++;
        }
        *length = k;
    }

    if (CMS) {
        float *mean;
        mean = (float *)calloc(Ndim,sizeof(float));
        for (i=0; i<*length; i++)
            for (j=0; j<Ndim; j++) mean[j] += lbuff[i][j];
        for (j=0; j<Ndim; j++) mean[j] /= (float)(*length);
        for (i=0; i<*length; i++)
            for (j=0; j<Ndim; j++) lbuff[i][j] -= mean[j];
        free(mean);
    }
    if (DELTA) {
        float **tmp;
        *length -= 4;
        *buff = tmp = (float **) alloc2d(*length,2*Ndim,sizeof(float));
        for (i=0; i<*length; i++) {
            for (j=0; j<Ndim; j++) tmp[i][j] = lbuff[i+2][j];
            for (j=0; j<Ndim; j++)
                tmp[i][j+Ndim] = lbuff[i+4][j] - lbuff[i][j];
        }
        free2d((void**)lbuff);
    }

    return 0;
}

int main(int argc, char **argv)
{
    float   **vector,**buff;
    float   **mean,**variance,*c;
    float   FILTERLOWEGY=0;
    int     numspch,numvecs,Ndim,Nmodes,maxlength=0,CMS=0, DELTA=0;
    int     i,j,length;
    char    ctlfile[512],*infileext=NULL,outfile[512];
    char    basefile[512],filename[512], sw, *indir=NULL;

    FILE    *fopen(),*listfile;
    float   scale();
    char    **alloc2d();

    void    estimate_multi_modals(),vector_quantize();
    int    read_backup_distribution(),store_distribution();


    /*
     * by default I assume the programs is to be run from data
     * only with no initial set of Gaussians. i.e. go to vq and them 
     * em.
     */
    /* SET DEFAULTS */
    Ndim = 13;

	/* Print out command, just to be nice */
	for (i=0; i<argc; i++) printf("%s ",argv[i]); printf("\n"); fflush(stdout);
    if (argc == 1)
    {
        printf(DESCRIPTION);
        QUIT ((USAGE,argv[0]));
    }
    for (i=1;i<argc;++i)
    {
        if (argv[i][0]!='-')
            QUIT ((USAGE,argv[0]));
        sw = argv[i][1];
        switch(sw)
        {
            case 'e': case 'E':
            {
                infileext = argv[++i];
                break;
            }
            case 'i':
            {
                indir = argv[++i];
                break;
            }
            case 'c': case 'C':
            {
                strcpy(ctlfile,argv[++i]);
                break;
            }
            case 'o': case 'O':
            {
                strcpy(outfile,argv[++i]);
                break;
            }
            case 'm': case 'M':
            {
                Nmodes = atoi(argv[++i]);
                break;
            }
            case 'd': case 'D':
            {
                Ndim = atoi(argv[++i]);
                break;
            }
            case 'n':
            {
                CMS = 1;
                break;
            }
            case 'l':
            {
                DELTA = 1;
                break;
            }
            case 'r': 
            {
                sscanf(argv[++i],"%f",&FILTERLOWEGY);
                break;
            }
            default:
            {
                printf(DESCRIPTION);
                QUIT ((USAGE,argv[0]));
            }

        }
    }


    listfile = fopen(ctlfile,"r");
    if (listfile == NULL) 
        QUIT(("Unable to open control file %s\n",ctlfile));

    numvecs = 0; maxlength = 0;
    while (fscanf(listfile,"%s\n",basefile) != EOF)
    {
        /*
         * If an extension has been specified for the files, append it
         */
        if (infileext) {
            if (indir) sprintf(filename,"%s/%s.%s",indir,basefile,infileext);
            else sprintf(filename,"%s.%s",basefile,infileext);
        }
        else {
            if (indir) sprintf(filename,"%s/%s",indir,basefile);
            else sprintf(filename,"%s",basefile);
        }

        if (read_cep(filename,&buff,&length,Ndim) == -1)
        {
        /*
         * What happens on a datafile read error? We can either just go on to
         * the next cepfile, or stop the program.  Set switch in header.h.
         */
#ifdef IGNORE_READ_ERR
            printf("READ FAILURE; skipping cepfile\n");
            fflush(stdout);
#else
            QUIT(("Unable to read %s\n",filename));
#endif
        }
        else
        {
            if (!DELTA)
                maxlength += length;
            else
                maxlength += length - 4;
	    free2d((void**)buff);
        }
    }
    printf("%d vectors in all\n",maxlength);
    rewind(listfile);

    /*
     * For our particular case it is needed that the Maximum dimensionality
     * is the same as the actual dimension of the data. Else the read_cep
     * routine fails!
     */
    if (DELTA) 
        vector = (float **) alloc2d(maxlength,2*Ndim,sizeof(float));
    else
        vector = (float **) alloc2d(maxlength,Ndim,sizeof(float));
    if(!vector) QUIT(("Unable to Allocate Space for Vector Array\n"));

    while ((fscanf(listfile,"%s\n",basefile) != EOF) && (numvecs <= maxlength))
    {
        /*
         * If an extension has been specified for the files, append it
         */
        if (infileext) {
            if (indir) sprintf(filename,"%s/%s.%s",indir,basefile,infileext);
            else sprintf(filename,"%s.%s",basefile,infileext);
        }
        else {
            if (indir) sprintf(filename,"%s/%s",indir,basefile);
            else sprintf(filename,"%s",basefile);
        }
        if (read_cep(filename,&buff,&length,Ndim) == -1)
        {
        /*
         * What happens on a datafile read error? We can either just go on to
         * the next cepfile, or stop the program.  Set switch in header.h.
         */
#ifdef IGNORE_READ_ERR
            printf("READ FAILURE; skipping cepfile\n");
            fflush(stdout);
#else
            QUIT(("Unable to read %s\n",filename));
#endif
        }
        else
        {
#ifdef DEBUG
	    printf("Success in reading %s feature file\n", filename);
#endif
            computefeat(&buff,CMS,DELTA,FILTERLOWEGY,Ndim,&length);
            for (i=0;i<length;i++)
	    {
		for (j=0;j<Ndim;++j)
		    vector[numvecs][j] = buff[i][j];
                if (DELTA)
		    for (j=0;j<Ndim;++j)
		        vector[numvecs][j+Ndim] = buff[i][j+Ndim];
                ++numvecs;
	    }
	    free2d((void**)buff);
        }
        if (numvecs > maxlength)
            QUIT(("**** Too many frames? Check your dimensionality!! ****\n"));
    }
    fclose(listfile);

    if (numvecs == 0)
        QUIT(("Only 0 vectors to compute a DISTRIBUTION!\nI am quitting!\n"));

    numspch = numvecs;
    if (DELTA) Ndim *= 2;

    c = (float *) malloc(Nmodes*sizeof(float));
    mean = (float **) alloc2d(Nmodes,Ndim,sizeof(float));
    if((variance = (float **) alloc2d(Nmodes,Ndim,sizeof(float))) == NULL)
        QUIT(("Unable to allocate space for variances\n"));

    estimate_mixture_gaussian(vector,numspch,Ndim,Nmodes,
                              mean,variance,c,outfile,0);
    if (store_distribution(outfile,Nmodes,Ndim,c,mean,variance) != 0) {
        printf("Unable to open %s to store distribution\n",outfile);
        exit (0);
    }

    free2d((void**)vector);
    free (c);
    free2d((void**)mean);
    free2d((void**)variance);

    return 0;
}
