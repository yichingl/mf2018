/*---------------------------------------------------------------------------
* Program that actually trains a gaussian mixture distribution using EM
* given a control file and the parameters (no. of mixture components needed)
* The mfc file is read, and a combination of mfc and delta-mfc features is
* obtained and the joint feature is used for the distribution estimation.
*----------------------------------------------------------------------------*/

#include "trainclass.h"


#define USAGE "USAGE:\n%s -c <control file> \\\n\t-m <for CMN> \\\n\t-l <for delta feats> \\\n\t-g <input model file> \\\n\t-o <output model file> \\\n\t-a <interpolation alpha(def 0.5)>\n"

int CMN, DELTA;

void  getfeatures(char*,float***,int*,int*,int);
void  initialize(float**,int,int,int,float**,float***,float***);
void  estimate_mixture_gau(float**,int,int,int,float**,float**,float*,float**,float**,float*,int);
void  store_distribution(char*,int,int,float*,float**,float**);
void  alloc_params(int K, int Ndim, float ***nmean, float ***nvar, float **nmixw);

int main(int argc, char **argv)
{
    float **feat,*mixw,**mean,**vars;
    float *nmixw,**nmean,**nvars, alpha;
    int  nframes,featdim,nmodes,cepdim,i,j;
    char *modelfile,*ctlfile,*outdir,*indir,sw,filename[1024],filehead[1024];
    FILE *fp;

    CMN = 0; DELTA = 0; alpha = 0.5;
    if (argc == 1) {
        printf(USAGE,argv[0]); exit (0);
    }
    for (i=1;i<argc;++i) {
        if (argv[i][0]!='-')
            QUIT ((USAGE,argv[0]));
        sw = argv[i][1];
        switch(sw)
        {
            case 'c': case 'C': {
                ctlfile = argv[++i];
                break;
            }
            case 'm': {
                CMN = 1;
                break;
            }
            case 'l': {
                DELTA = 1;
                break;
            }
            case 'i': case 'I': {
                indir = argv[++i];
                break;
            }
            case 'o': case 'O': {
                outdir = argv[++i];
                break;
            }
            case 'g': {
                modelfile = argv[++i];
                break;
            }
            case 'a': {
                sscanf(argv[++i],"%f",&alpha);
                break;
            }
            default: {
                printf(USAGE,argv[0]); exit(0);
            }

        }
    }

    read_classes(modelfile,&nmodes,&featdim,&mean,&vars,&mixw);
    alloc_params(nmodes,featdim,&nmean,&nvars,&nmixw);

    if ((fp = fopen(ctlfile,"r")) == NULL)
        QUIT(("Unable to open %s for reading\n",ctlfile))

    if (DELTA) {
        cepdim = featdim/2;
        if (2*cepdim != featdim) 
            QUIT(("Model dim. of %s odd; does not support deltas\n",modelfile))
    }
    getfeatures(ctlfile,&feat,&nframes,&featdim,cepdim); 


    estimate_mixture_gau(feat,nframes,featdim,nmodes,mean,vars,mixw,nmean,nvars,nmixw,1);
    for (i=0;i<nmodes;i++) {
        for (j=0; j<featdim; j++)
            nmean[i][j] = (1 - alpha)*mean[i][j] + alpha*nmean[i][j];
        nmixw[i] = (1 - alpha)*mixw[i][j] + alpha*nmixw[i][j];
    }
    store_distribution(outfile,nmodes,featdim,mixw,mean,vars);

    // Gotta free stuff here. Too lazy for now.
    return 0;
}



/*==========================================================================*/
/* Routine to read mfc file and get feature mfc+deltamfc */

void  getfeatures(char *ctlfile, char *infile, float ***feat, int *nfr, int *featdim, int cepdim, char ***tags)
{
    float **cep, **features=NULL, **outfeat, *cepmean;
    int i,j,k,beg,end,frmno,nfrm,outdim,totframes=0,begseg;
    char fname[1024], prevfile[1024], fnamebase[256], *dummy="/dev/null";
    FILE *fp;

    outdim = DELTA ? 2*cepdim : cepdim;
    if ((fp=fopen(ctlfile,"r"))==NULL) QUIT(("Unable to open %s\n",ctlfile))
    while(fscanf(fp,"%s %d %d",fname,&beg,&end)!=EOF)
        totframes += end - beg + 1;
    rewind(fp);

    if (CMN) cepmean = (float*) calloc(outdim, sizeof(float));

    outfeat = (float **)alloc2d(totframes,outdim,sizeof(float));
    *featdim = outdim; *feat = outfeat;

    strcpy(prevfile,dummy); frmno = 0;
    while(fscanf(fp,"%s %d %d",fname,&beg,&end)!=EOF){
        if (beg < 1) beg = 1;
        if (strcmp(fname,prevfile) != 0){
            printf("Reading from %s\n",fname); fflush(stdout);
            if (features != NULL) free2d((void **)features);
            if (read_cep(fname,&cep,&nfrm,cepdim) == IO_ERR)
                QUIT(("Unable to read %s\n",fname))
            features = (float **)alloc2d(nfrm,2*cepdim,sizeof(float));
            for (i=0;i<nfrm;i++){
                for(j=0;j<cepdim;j++) features[i][j] = cep[i][j];
            }
            if (DELTA) {
                for (i=2;i<nfrm-2;i++){
                    for(k=0,j=cepdim;k<cepdim;k++,j++) {
                        features[i][j] = 0.67*cep[i+1][k]+0.33*cep[i+2][k] -
                                 0.67*cep[i-1][k]-0.33*cep[i-2][k];
                    }
                }
                for(j=0;j<cepdim;j++) {
                    features[0][j] = cep[0][j];
                    features[1][j] = cep[1][j];
                    features[nfrm-1][j] = cep[nfrm-1][j];
                    features[nfrm-2][j] = cep[nfrm-2][j];
                    features[0][j+cepdim]=features[2][j+cepdim];
                    features[1][j+cepdim]=features[2][j+cepdim];
                    features[nfrm-1][j+cepdim]=features[nfrm-3][j+cepdim];
                    features[nfrm-2][j+cepdim]=features[nfrm-3][j+cepdim];
                }
            }
            free2d((void **)cep);
        }
        end = end >= nfrm ? DELTA ? nfrm - 2 : nfrm : end;
        if (end < beg)
            QUIT(("Beg %d, End %d and nfr %d error %s\n",beg,end,nfrm,fname))
        for (j=0;j<outdim;j++) cepmean[j] = 0; begseg = frmno;
        for (i=beg;i<end;i++){
            for (j=0;j<outdim;j++) outfeat[frmno][j] = features[i][j];
            if (CMN) {
                for (j=0;j<outdim;j++) cepmean[j] += features[i][j];
            }
            frmno++;
        }
        if (CMN) {
            for (j=0;j<outdim;j++) cepmean[j] /= (float)(end-beg+1);
            for (i=begseg;i<frmno;i++){
               for (j=0;j<cepdim;j++) outfeat[i][j] -= cepmean[j];
            }
        }
        strcpy(prevfile,fname);
    }
    if (features != NULL) free2d((void **)features);
    if (CMN) free(cepmean);
    printf("%d of %d total frames obtained from cepfiles\n",frmno,totframes);
    *nfr = frmno;
    return;
}

/*==========================================================================*/
// Routine to estimate mixture gaussians from an initial point
float gauss(float*,float*,float*,float,int);

void  alloc_params(int K, int Ndim, float ***nmean, float ***nvar, float **nmixw)
{
    float    **Newvar, **Newmean, *Newc;

    if ((Newmean = (float **)alloc2d(K,Ndim,sizeof(float))) == NULL)
        QUIT(("Unable to allocate Newmean array\n"));

    if ((Newvar = (float **) alloc2d(K,Ndim,sizeof(float))) == NULL)
        QUIT(("Unable to allocate for Newvar array\n"));

    if ((Newc = (float *) calloc(K,sizeof(float))) == NULL)
        QUIT(("Unable to allocate for Newc array\n"));

    *nmean = Newmean;
    *nvar = Newvar;
    *nmixw = Newc;

    return;
}

void    estimate_mixture_gau(
    float    **x,                /* The observation vectors */
    int      N,                  /* Number of observation vectors */
    int      Ndim,               /* Dimensionality of observations */
    int      K,                  /* Number of modes in the distribution */
    float    **omean,             /* The means of the various modes */
    float    **ovar,              /* The variances of all the modes */
    float    *oc,                 /* A-priori probabilities of each of the */
    float    **mean,             /* The means of the various modes */
    float    **var,              /* The variances of all the modes */
    float    *c,                 /* A-priori probabilities of each of the
                                    modes, or mixing proportion */
    int      numiters            /* Number of iterations of EM to run */
   )
{
    float    **Newvar,
             **invvar,
             **Newmean,
             *Newc,
             *Tau,
             *corprod,
             Const,
             SumNewc,
             tdiff,Prevlogprob,LogProb,Threshold,Improvement;

    int     i,j,k,iter=0;

    float   Expectation();
    char    **alloc2d();
    void    free2d();
            
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
    if ((Newmean = (float **)alloc2d(K,Ndim,sizeof(float))) == NULL)
        QUIT(("Unable to allocate Newmean array\n"));

    if ((Newvar = (float **) alloc2d(K,Ndim,sizeof(float))) == NULL)
        QUIT(("Unable to allocate for Newvar array\n"));

    if ((invvar = (float **) alloc2d(K,Ndim,sizeof(float))) == NULL)
        QUIT(("Unable to allocate for invvar array\n"));

    if ((Newc = (float *)malloc(K*sizeof(float))) == NULL)
        QUIT(("Unable to allocate %d floats for Newc array \n",K));

    if ((Tau = (float *)malloc(K*sizeof(float))) == NULL)
        QUIT(("Unable to allocate %d floats for Tau array \n",K));
    
    if ((corprod = (float *)malloc(K*sizeof(float))) == NULL)
        QUIT(("Unable to allocate %d floats for corprod array \n",K));

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
    for (k=0; k<K ; ++k)
    {
        c[k] = oc[k];
        for (j=0;j<Ndim;++j)
	{
            mean[k][j] = omean[k][j];
	    var[k][j] = ovar[k][j];
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
             if (var[k][j] < 1.0e-20) var[k][j] = 1.0e-20;
             invvar[k][j] = 1.0/var[k][j];
             corprod[k] -= 0.5*log((double)var[k][j]);
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
        Prevlogprob += Expectation(Tau,x[i],c,mean,invvar,corprod,K,Ndim);
        for (k=0; k<K ; ++k)
        {
            if (Tau[k] > 1.0e-07){
                Newc[k] += Tau[k];
                for (j=0;j<Ndim;++j)
                    Newmean[k][j] += Tau[k] * x[i][j];
            }
        }
    }
    for (k=0; k<K ; ++k){
        if (Newc[k] < 1.0e-20) Newc[k] = 1.0e-20;
        for (j=0;j<Ndim;++j)
	    Newmean[k][j] /= Newc[k];
    }
    for (i=0;i<N;++i)
    {
        Expectation(Tau,x[i],c,mean,invvar,corprod,K,Ndim);
        for (k=0; k<K ; ++k)
        {
            if (Tau[k] > 1.0e-07){
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
        SumNewc = 0;
        for (k=0;k<K;++k)
            SumNewc += Newc[k];

        for (k=0; k<K ; ++k)
        {
            for (j=0;j<Ndim;++j)
	    {
                mean[k][j] = Newmean[k][j];
//		var[k][j] = Newvar[k][j]/Newc[k]; // Not updating variances
            }
            c[k] = Newc[k]/SumNewc;
        }

        /* * Initialize all New values to 0 */
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

/* Not updating variances
        for (k=0;k<K;++k)
        {
            corprod[k] = log((double)c[k]);
            for (j=0;j<Ndim;++j){
                if (var[k][j] < 1.0e-20) var[k][j] = 1.0e-20;
                invvar[k][j] = 1.0/var[k][j];
                corprod[k] -= 0.5*log((double)var[k][j]);
            }
        }
*/

        LogProb = 0;
        for (i=0;i<N;++i)
        {
            LogProb += Expectation(Tau,x[i],c,mean,invvar,corprod,K,Ndim);
            for (k=0; k<K ; ++k)
            {
                if (Tau[k] > 1.0e-07){
                    Newc[k] += Tau[k];
                    for (j=0;j<Ndim;++j){
                        Newmean[k][j] += Tau[k] * x[i][j];
                        tdiff = x[i][j] - mean[k][j];
		        Newvar[k][j] += Tau[k] * tdiff * tdiff;
                    }
                }
            }
        }
        for (k=0; k<K ; ++k){
            if (Newc[k] < 1.0e-20) Newc[k] = 1.0e-20;
            for (j=0;j<Ndim;++j)
                Newmean[k][j] /= Newc[k];
        }

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
    free2d((void **)Newvar);
    free2d((void **)invvar);
    free2d((void **)Newmean);
    free(Newc);
    return;
}


/*----------------------------------------------------------------------------
  Routine to compute the a-posteriori probabilities of each of the modes in
  the Mixture distribution defined by c,mean and var for the the vector
  defined by x[]. Returns the posteriors in Tau[]
------------------------------------------------------------------------------*/
float    Expectation(Tau,x,c,mean,invvar,corprod,K,Ndim)
    float   Tau[],           /* Mode probablity of observation vector x */
            x[],             /* Observation vector x */
            c[],             /* Mixing proportions of modes */
            **mean,          /* Means of the modes */
            **invvar,        /* Variances of the modes */
            corprod[];       /* Precomputed 1/sqrt(2PI**Ndim * mod(variance)) */

    int     K,               /* Number of modes */
            Ndim;            /* Dimensionality of observation/distribution */
{
    float    Probability,
             max,sum;

    int      i;

    Tau[0] = max = gauss(x,mean[0],invvar[0],corprod[0],Ndim);
    for (i=1;i<K;++i) {
        /* Note: prior c[i] has already been factored into corprod */
        Tau[i] = gauss(x,mean[i],invvar[i],corprod[i],Ndim);
        if (Tau[i] > max) max = Tau[i];
    }

    /* convert Tau from log to linear */
    sum = 0.0;
    for (i=0;i<K;++i) {
         if (Tau[i] > max-15) Tau[i] = exp(Tau[i] - max + 10);
         else Tau[i] = 0;
         sum += Tau[i];
    }
    for (i=0;i<K;++i) Tau[i] /= sum;

    Probability = max + log(sum) - 10;
    return(Probability);
}


float gauss(float *x,float *mean,float *invvar,float corprod,int Ndim)
{
    auto int	j;
    auto float	power=0,density;

    for (j=0;j<Ndim;++j)
	power += (x[j]-mean[j])*(x[j]-mean[j])*invvar[j];
    power /= -2.0;
    density = corprod + power;
    return(density);
}


/*========================================================================*/
/* Routine to store the distributions computed */

void  store_distribution(char *file, int Nmodes, int Ndim, 
                         float *c, float **mean, float **var)
{
    int     i,j;
    FILE    *fopen(),*outputfile;

    if ((outputfile = fopen(file,"w")) == NULL) {
        printf("Unable to open %s for writing!\n",file);
        printf("Will store models in /tmp/DISTRIB instead\n");
        if ((outputfile = fopen("/tmp/DISTRIB","w")) == NULL) 
            QUIT(("Unable to open /tmp/DISTRIB either! No output!\n"))
    }
    
    fprintf(outputfile,"%d %d\n",Nmodes,Ndim);
    for (i=0;i<Nmodes;++i)
    {
        fprintf(outputfile,"%e\n",c[i]);
        for (j=0;j<Ndim;++j)
            fprintf(outputfile,"%e ",mean[i][j]);
        fprintf(outputfile,"\n");
        for (j=0;j<Ndim;++j)
            fprintf(outputfile,"%e ",var[i][j]);
        fprintf(outputfile,"\n");
    }
    fclose(outputfile);
    return;
}

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
