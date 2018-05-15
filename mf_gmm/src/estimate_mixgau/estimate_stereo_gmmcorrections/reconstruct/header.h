/*#define maxdim 13 */

#define QUIT(X) {printf("%s(%d):",__FILE__,__LINE__);\
                 printf X;fflush(stdout);exit(0);}

void estimate_gmm_stereo_corrections(float **x,   /* observation vectors */
                                     float ***y,  /* Reconstruction */
                                     int N,       /* No. of vectors */
                                     int Ndim,    /* Dimensionality of data */
                                     int K,       /* No. of Gaussians */
                                     float **mean,/* Means of distribution */
                                     float **var, /* Variances */
                                     float *c,    /* Mixture weights */
                                     int   R,     /* "Lookahead frames" */
                                     int   L,     /* "Lookbehind frames" */
                                     float ***meany /* Means of corrections */
                                    );

float    Expectation(float *Tau,   /* Mode probablity of observation vector x */
                     float *x,     /* Observation vector x */ 
                     float *c,     /* Mixing proportions of modes */
                     float **mean, /* Means of the modes */
                     float **hafinvvar,  /* 1.0/(2.0*variance) of the modes */
                     float *corprod, /* 1/sqrt(2PI^Ndim*|variance|) */
                     int K,        /* Number of modes */
                     int Ndim);    /* Dimensionality of data */

float gauss(float *x,
            float *mean,
            float *hafinvvar,
            float corprod,
            int Ndim);


int    read_input_distribution(char *filen,
                               float ***mean,
                               float ***variance,
                               float **pi,
                               int *Nmodes,
                               int *Ndim,
                               int *L, int *R,
                               float ****meany);

int write_file(int ncopies, int nfr, int ndim, float ***v, char *fname);


char **alloc2d(int dim1,int dim2, int size);
void free2d(void **p);

char ***alloc3d(int dim1,int dim2,int dim3,int size);
void free3d(void ***X);

int read_cep(char *file, float***cep, int *numframes, int cepsize);
