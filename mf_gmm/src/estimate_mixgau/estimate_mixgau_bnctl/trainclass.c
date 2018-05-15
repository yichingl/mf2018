/*---------------------------------------------------------------------------
* Program that actually trains a gaussian mixture distribution using EM
* given a control file and the parameters (no. of mixture components needed)
* The mfc file is read, and a combination of mfc and delta-mfc features is
* obtained and the joint feature is used for the distribution estimation.
*----------------------------------------------------------------------------*/
#include <ayesha.h>

void  getfeatures(char*,float***,int*,int*);
void  initialize(float**,int,int,int,float**,float***,float***);
void  estimate_mixture_gau(float**,int,int,int,float**,float**,float*,int);
void  store_distribution(char*,int,int,float*,float**,float**);
void  read_distribution(char *,int*,int*,float**,float***,float***);

int main(int argc, char **argv)
{
    float **feat,*mixw,**mean,**vars;
    int  nframes,featdim,nmodes,modeldim;
    char *ctlfile,*outfile, *initialdistfile;

    if (argc != 4 && argc != 5)
        QUIT(("USAGE:\n %s <ctlfile> <nmodes> <outfile> [initialdist]\nIf you use initialdist, the actual value of nmodes is read from it\n", argv[0]))
    ctlfile = argv[1]; sscanf(argv[2],"%d",&nmodes);outfile=argv[3];
    if (argc == 5)
    	initialdistfile = argv[4];

    getfeatures(ctlfile,&feat,&nframes,&featdim); 
    if (argc == 4)
        initialize(feat,nframes,featdim,nmodes,&mixw,&mean,&vars);
    else{
        read_distribution(initialdistfile,&nmodes,&modeldim,&mixw,&mean,&vars);
	if (modeldim != featdim)
            QUIT(("Mismatch between data dim %d and model dim %d\n",
	           featdim, modeldim))
    }
    estimate_mixture_gau(feat,nframes,featdim,nmodes,mean,vars,mixw,10);
    store_distribution(outfile,nmodes,featdim,mixw,mean,vars);
    return 0;
}



/*==========================================================================*/
/* Routine to read mfc file and get feature mfc+deltamfc */

#define CEPDIM 13
void  getfeatures(char *ctlfile, float ***feat, int *nfr, int *featdim)
{
    float **cep, **features=NULL, **outfeat, cepmean[100];
    int i,j,k,beg,end,frmno,nfrm,cepdim,outdim,totframes=0,begseg;
    char fname[256],prevfile[256],*dummy="/dev/null";
    FILE *fp;

    cepdim = (int)CEPDIM; outdim = cepdim;
    if ((fp=fopen(ctlfile,"r"))==NULL) QUIT(("Unable to open %s\n",ctlfile))
    while(fscanf(fp,"%s %d %d",fname,&beg,&end)!=EOF)
        totframes += end - beg + 1;
    rewind(fp);

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
            free2d((void **)cep);
        }
        if (end >= nfrm) end = nfrm;
        if (end <= beg) continue; // Ignore segment
//            QUIT(("Beg %d, End %d and nfr %d error %s\n",beg,end,nfrm,fname))
        for (j=0;j<outdim;j++) cepmean[j] = 0; begseg = frmno;
        for (i=beg;i<end;i++){
            for (j=0;j<outdim;j++) outfeat[frmno][j] = features[i][j];
//            for (j=0;j<outdim;j++) cepmean[j] += features[i][j];
            frmno++;
        }
//        for (j=0;j<outdim;j++) cepmean[j] /= (float)(end-beg+1);
//        for (i=begseg;i<frmno;i++){
//           for (j=0;j<cepdim;j++) outfeat[i][j] -= cepmean[j];
//        }
        strcpy(prevfile,fname);
    }
    if (features != NULL) free2d((void **)features);
    printf("%d of %d total frames obtained from cepfiles\n",frmno,totframes);
    *nfr = frmno;
    return;
}

/*=========================================================================*/
/* The  Routine that intializes the mixture gaussian */

void  vector_quantize(float**,int,float**,int,int,int*);    
void initialize(
    float **vector,       // Vector sequence for mixgau
    int   numvecs,        // Number of vectors
    int   Ndim,           // Dimension of vectors
    int   Nmodes,         // Number of modes in Gaussian Mixture
    float **omixw,        // Mixture weights
    float ***omean,       // Means of modes
    float ***ovars        // Variances of modes
)
{
    float   **mean,**variance,*c;
    float   atemp;
    int     i,j,k,*nbin,*bin;

    /*  allocate the mean and variance and c arrays.  */ 
    c = (float *) malloc(Nmodes*sizeof(float));
    mean = (float **) alloc2d(Nmodes,Ndim,sizeof(float));
    if((variance = (float **) alloc2d(Nmodes,Ndim,sizeof(float))) == NULL)
        QUIT(("Unable to allocate space for variances\n"));

    nbin = (int *)calloc(Nmodes,sizeof(int)); /* no of vectors in a mode */
    if (nbin==NULL)
       QUIT(("Unable to allocate %d space for nbin array\n",Nmodes));

    /* The vector_quantize routine performs VQ with a mahalonobis metric
     * and returns the codes as the means and the binning information for
     * each of the vectors. We use these to find variances and mixture weights
     */

    if ((bin = (int *)calloc(numvecs,sizeof(int))) == NULL)
       QUIT(("Unable to allocate %d ints for bin array\n",numvecs));

    vector_quantize(mean,Nmodes,vector,numvecs,Ndim,bin);

    for (i=0;i<Nmodes;++i) c[i] = 1.0/(float)Nmodes;
    for (k=0;k<Nmodes;++k) {
        nbin[k] = 0;
        for (i=0;i<Ndim;++i) variance[k][i] = 0;
    }
    for (i=0;i<numvecs;++i) {
        for (j=0;j<Ndim;++j) {
            atemp = (vector[i][j] - mean[bin[i]][j]);
            variance[bin[i]][j] += atemp*atemp;
        }
        ++nbin[bin[i]];
    }
    for (k=0;k<Nmodes;++k)
        for(j=0;j<Ndim;++j) variance[k][j] /= nbin[k];

    free(bin);   
    free(nbin); 
	
    *omixw = c; *omean = mean; *ovars = variance;
    return;
}

/*==========================================================================*/
/* THE VQ ROUTINE */

void  vector_quantize(
    float **codes,           // The VQ codewords 
    int Ncodes,              // Number of VQ codewords
    float **vector,          // The vectors that we compute the codewords off
    int Nvecs,               // Number of vectors
    int Ndim,                // Dimensionality of the vectors
    int *bin                 // Array indicting codeword each vector belongs to
   )    
{
    float    distortion,improvement,prevdist,threshold;
    int    iter=0;
    extern void    Initial_codebook(),Update_codebook();
    extern float    Distortion_and_cluster();

    improvement = 100;
    threshold = 1e-06;


    Initial_codebook(Ncodes,codes,Nvecs,vector,Ndim,bin);
    prevdist = Distortion_and_cluster(bin,vector,Nvecs,codes,Ncodes,Ndim);

    printf("Initial Distortion = %f\n",prevdist);
    while (((improvement > threshold)||(improvement < 0))&&(iter<5))
    {
        Update_codebook(bin,vector,Nvecs,codes,Ncodes,Ndim);
        distortion = 
            Distortion_and_cluster(bin,vector,Nvecs,codes,Ncodes,Ndim);
        improvement = (prevdist-distortion)/distortion;
        prevdist = distortion;
        printf("VQ distortion = %f, improvement = %f after %d iterations\n",distortion,improvement,iter); fflush(stdout);
        ++iter;
    }
    printf("VQ exited after %d iterations with total distortion %f\n",iter,prevdist); fflush(stdout);

    return;
}


/*------------------------------------------------------------------------------
    Function for generating initial codebook
------------------------------------------------------------------------------*/
void Initial_codebook(Ncodes,codes,Nvecs,vector,Ndim,bin)
    float    **codes,**vector;
    int      Ncodes,Nvecs,Ndim,*bin;
    {
    auto int    i,j,nbin;
    auto float  mindist,codedist;
    float       prune_distance();

    for (i=0;i<Ncodes;++i){
        for (j=0;j<Ndim;++j)
            codes[i][j] = vector[i][j];
        bin[i] = 1;
    }

    for (i=Ncodes;i<Nvecs;++i){
        nbin = 0; mindist = 1.0e+32;
        mindist = prune_distance(vector[i],codes[0],Ndim,mindist);
        for (j=1;j<Ncodes;++j){
            codedist = prune_distance(vector[i],codes[j],Ndim,mindist);
            if (codedist < mindist){
                mindist = codedist;
                nbin = j;
            }
        }
        for (j=0;j<Ndim;++j){
            codes[nbin][j] = codes[nbin][j]*bin[nbin]+vector[i][j];
            codes[nbin][j] /= (bin[nbin]+1);
        }
        ++bin[nbin];
    }
    return;
}


/*------------------------------------------------------------------------------
    Function to cluster the training vectors around the codes
    and calculate the sum distortion
------------------------------------------------------------------------------*/

float  Distortion_and_cluster(bins,vector,Nvecs,codes,Ncodes,Ndim)
    float    **vector,**codes;
    int    Nvecs,Ncodes,*bins,Ndim;
    {
    auto int    i,j,nearest;
    auto float    distortion,codedist,mindist;
    float        prune_distance();

    distortion = 0;
    for (i=0;i<Nvecs;++i){
	mindist = 1.0e+32;
        nearest = 0;
        mindist = prune_distance(vector[i],codes[0],Ndim,mindist);
        for (j=1;j<Ncodes;++j){
            codedist = prune_distance(vector[i],codes[j],Ndim,mindist);
            if (codedist < mindist){
                mindist = codedist;
                nearest = j;
            }
        }
        bins[i] = nearest;
        distortion += mindist;
    }
    return(distortion);
}


/*------------------------------------------------------------------------------
    Function to update codebook according to the clustering information
    receive through the array bin
------------------------------------------------------------------------------*/

void  Update_codebook(bin,vector,Nvecs,codes,Ncodes,Ndim)
    float    **vector,**codes;
    int    Nvecs,Ncodes,Ndim,*bin;
    {
    auto int    i,j,k,l,*Noincode;
    auto int    nfull,nmpt,*Mpt,*Full;
    int         Largest();

    if ((Noincode = (int *)malloc(Ncodes*sizeof(int))) == NULL)
        QUIT(("Unable to allocate space for array Noincode\n"));

    if ((Mpt = (int *)malloc(Ncodes*sizeof(int))) == NULL)
        QUIT(("Unable to allocate space for array Mpt\n"));

    if ((Full = (int *)malloc(Ncodes*sizeof(int))) == NULL)
        QUIT(("Unable to allocate space for array Full\n"));

    for (i=0;i<Ncodes;++i){
        Noincode[i] = 0;
        for (j=0;j<Ndim;++j)
            codes[i][j] = 0;
    }
    
    for (i=0;i<Nvecs;++i){
        j = bin[i];
        for (k=0;k<Ndim;++k)
            codes[j][k] += vector[i][k];
        ++Noincode[j];
    }

    nfull = 0;
    nmpt = 0;
    for (i=0;i<Ncodes;++i){
        if (Noincode[i] > 1){
            Full[nfull] = i;
            ++nfull;
            for (j=0;j<Ndim;++j)
                codes[i][j] /= (float)Noincode[i];
        }
        else{
            Mpt[nmpt] = i;
            ++nmpt;
        }
    }
    
    if (nmpt<=nfull){
        for (i=0;i<nmpt;++i){
            j = Mpt[i];
            l = Largest(&nfull,Full,Noincode);
            for (k=0;k<Ndim;++k){
                codes[j][k] = codes[l][k] * 0.90;
                codes[l][k] *= 1.10;
            }
        }
    }
    else
        QUIT(("CATASTROPHE, CATASTROPHE\nRerun with better initialization\n"));


    free(Noincode);
    free(Mpt);
    free(Full);
    return;
}

/* Function to obtain largest bin in case of an empty bin, to split */

int Largest(nfulptr,Full,Noincode)
    int    *Full,*Noincode,*nfulptr;
    {
    auto int i,j,ibig,ifull;

    ifull = 0; ibig = Full[ifull];
    for (i=1;i<*nfulptr;++i){
        j = Full[i];
        if (Noincode[j]>Noincode[ibig]){
            ibig = j;
            ifull = i;
        }
    }

    for (i=ifull;i<*nfulptr-1;++i)
        Full[i] = Full[i+1];
    --(*nfulptr);
    
    return(ibig);
}


/*------------------------------------------------------------------------------
    This function defines the distance metric being used and
    can be varied as per requirement. Metric used here euclidean.
    This particular function also checks if the distance exeeds a
    previously given minimum and returns if it does. Useful to
    speed up clustering..
------------------------------------------------------------------------------*/

float  prune_distance(x,y,Ndim,mindist)
    float    x[],y[],mindist;
    int    Ndim;
    {
    auto int    i;
    auto float    dist;
    dist = (x[0]-y[0])*(x[0]-y[0]);
    for (i=1;i<Ndim;++i){
        dist += (x[i]-y[i])*(x[i]-y[i]);
        if (dist > mindist) break;
    }
    return(dist);
}

/*==========================================================================*/
// Routine to estimate mixture gaussians from an initial point
float gauss(float*,float*,float*,float,int);

void    estimate_mixture_gau(
    float    **x,                /* The observation vectors */
    int      N,                  /* Number of observation vectors */
    int      Ndim,               /* Dimensionality of observations */
    int      K,                  /* Number of modes in the distribution */
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
		var[k][j] = Newvar[k][j]/Newc[k];
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

        for (k=0;k<K;++k)
        {
            corprod[k] = log((double)c[k]);
            for (j=0;j<Ndim;++j){
                if (var[k][j] < 1.0e-20) var[k][j] = 1.0e-20;
                invvar[k][j] = 1.0/var[k][j];
                corprod[k] -= 0.5*log((double)var[k][j]);
            }
        }

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

/*========================================================================*/
/* Routine to read initial distributions, if any */

void  read_distribution(char *file, int *nmodes, int *ndim, 
                         float **c, float ***mean, float ***var)
{
    int     i,j;
    float *lc, **lmean, **lvar;
    FILE    *fopen(),*inputfile;

    if ((inputfile = fopen(file,"r")) == NULL) 
        QUIT(("Unable to open %s for reading\n",file))
    
    fscanf(inputfile,"%d %d",nmodes,ndim);
    *c = lc = (float *)calloc(*nmodes,sizeof(float));
    *mean = lmean = (float**) alloc2d(*nmodes,*ndim,sizeof(float));
    *var = lvar = (float**) alloc2d(*nmodes,*ndim,sizeof(float));

    for (i=0;i<*nmodes;++i)
    {
        fscanf(inputfile,"%f",&lc[i]);
        for (j=0;j<*ndim;++j)
            fscanf(inputfile,"%f",&lmean[i][j]);
        for (j=0;j<*ndim;++j)
            fscanf(inputfile,"%f",&lvar[i][j]);
    }
    fclose(inputfile);
    return;
}
