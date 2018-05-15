#include "discriminative.h"

float chklog(x) {float y; y = x > 1e-20 ? log(x) : log(1e-20); return y;}

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
                         float *classprior)
{
    int i, j, k;
    float ***newmean, ***newvar, **newmixw;
    float ***halfinvvar, **corprod, *logclassprior, D;
    float **firstmoment, **secondmoment, *thetaclass, *thetaglobal, *gamma;

    // Allocate for new means vars and mixws
    newmean = (float***)alloc3d(nclass,nmodes[0],ndim,sizeof(float));
    *onewmean = newmean;
    newvar = (float***)alloc3d(nclass,nmodes[0],ndim,sizeof(float));
    *onewvar = newvar;
    newmixw = (float**)alloc2d(nclass,nmodes[0],sizeof(float));
    *onewmixw = newmixw;
    
    // First compute halfinvvar, corprod, and logpriors  for all models
    halfinvvar = (float ***)alloc3d(nclass, nmodes[0], ndim, sizeof(float));
    corprod = (float **)alloc2d(nclass, nmodes[0], sizeof(float));
    logclassprior = (float *)calloc(nclass, sizeof(float));
    for (i=0; i<nclass; i++) {
        logclassprior[i] = chklog(classprior[i]);
        for (j=0; j<nmodes[i]; j++) {
            corprod[i][j] = chklog(mixw[i][j]);
            for (k=0; k<ndim; k++) {
                halfinvvar[i][j][k] = 0.5/var[i][j][k];
                corprod[i][j] -= 0.5*chklog(var[i][j][k]);
            }
        }
    }
    for (i=0; i<nclass; i++) {
        printf("Reestimating %dth class.. ",i); fflush(stdout);
        // First empty buffers
        firstmoment = (float **)alloc2d(nmodes[i],ndim,sizeof(float));
        secondmoment = (float **)alloc2d(nmodes[i],ndim,sizeof(float));
        thetaclass = (float *)calloc(nmodes[i],sizeof(float));
        thetaglobal = (float *)calloc(nmodes[i],sizeof(float));
        gamma = (float *)calloc(nmodes[i],sizeof(float));

        for (j=0; j<nseqs[i]; j++) {
            accumulateforclass(data[i][j], seqlen[i][j],
                               mean, halfinvvar, mixw, corprod,
                               logclassprior, nclass, nmodes[i], ndim,
                               i, firstmoment, secondmoment,
                               thetaclass, thetaglobal, gamma);
        }
        D = estimateD(mean[i],var[i],firstmoment,secondmoment,gamma,
                nmodes[i],ndim);
        // Reestimate parameters
        computeclassparams(mean[i],var[i],mixw[i],
                           newmean[i],newvar[i],newmixw[i],
                           nmodes[i],ndim,firstmoment,secondmoment,
                           gamma,thetaclass,thetaglobal,D);

        free2d((void**)firstmoment); free2d((void**)secondmoment);
        free(thetaclass); free(thetaglobal); free(gamma);
        printf("done!\n"); fflush(stdout);
    }
    free3d((void***)halfinvvar);free2d((void**)corprod);free(logclassprior);

    return;
}
