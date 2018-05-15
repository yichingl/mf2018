#include <ayesha.h>

void computeposteriors(float *x,
                       float ***mean,       // Class, gaussian, dim
                       float ***halfinvvar, // Class, gaussian, dim
                       float  **mixw,       // Class, gaussian
                       float  **corprod,    // Class, gaussian
                       float  *logclassprior,
                       int    nclass,
                       int    ngauss,
                       int    ndim,
                       int    classid,      // The class for the current vector
                       float  *withinclassposterior,
                       float  *betweenclassposterior);


void accumulateforclass(float **x,
                        int   ndata,
                        float ***mean,       // Class, gaussian, dim
                        float ***halfinvvar, // Class, gaussian, dim
                        float **mixw,        // Class, gaussian
                        float **corprod,     // Class, gaussian
                        float *logclassprior,
                        int   nclass,
                        int   ngauss,
                        int   ndim,
                        int   classid,       // The class for the data
                        float **firstmoment,
                        float **secondmoment,
                        float *thetaclass,
                        float *thetaglobal,
                        float *gamma);


float estimateD(float **mean,
                float **var,
                float **firstmoment,
                float **secondmoment,
                float *gamma,
                int   ngauss,
                int   ndim);

void computeclassparams(float **mean,
                        float **var,
                        float *mixw,
                        float **newmean,
                        float **newvar,
                        float *newmixw,
                        int   ngauss,
                        int   ndim,
                        float **firstmoment,
                        float **secondmoment,
                        float *gamma,
                        float *thetaclass,
                        float *thetaglobal,
                        float D);


void reestimateallmodels(float ****data,
                         int   nclass,
                         int   *nseqs,
                         int   **seqlen,
                         int   ndim,
                         float ***mean,
                         float ***var,
                         float **mixw,
                         int   *nmodes,
                         float ****onewmean,
                         float ****onewvar,
                         float ***onewmixw,
                         float *classprior);


void readallclassdata(char *ctllist, int dim,
                      float *****odata, // nclass, nseqs, len, dim
                      int *onclass,
                      int **onseqs,     // nclass
                      int ***olenseq);  // nclass, nseqs

void readalldistributions(char *ctllist, int *onclass,
                          int **onmodes, int *ondim,
                          float ***oc, float ****omean, float ****ovar);

void storealldistributions(char *ctllist, int nclass,
                           int *nmodes, int ndim,
                           float **c, float ***mean, float ***var);

void readpriors(char *priorsfile, int nclass, float **oclasspriors);
