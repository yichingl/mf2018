#include "discriminative.h"

void accumulateforclass(float **x,
                        int   ndata,
                        float ***mean,       // Class, gaussian, dim
                        float ***halfinvvar, // Class, gaussian, dim
                        float **mixw,       // Class, gaussian
                        float **corprod,    // Class, gaussian
                        float *logclassprior,
                        int   nclass,
                        int   ngauss,
                        int   ndim,
                        int   classid,      // The class for the data
                        float **firstmoment,
                        float **secondmoment,
                        float *thetaclass,
                        float *thetaglobal,
                        float *gamma)
{
    int t, i, j;
    float  *withinclassposterior, *betweenclassposterior, gammai;

    withinclassposterior = (float *)calloc(ngauss, sizeof(float));
    betweenclassposterior = (float *)calloc(ngauss, sizeof(float));

    for (t=0; t<ndata; t++) {
        computeposteriors(x[t],mean,halfinvvar,mixw,corprod,
                          logclassprior,nclass,ngauss,ndim,classid,
                          withinclassposterior,betweenclassposterior);

        for (i=0; i<ngauss; i++) {
            thetaclass[i] += withinclassposterior[i];
            thetaglobal[i] += betweenclassposterior[i];

            gammai = withinclassposterior[i] - betweenclassposterior[i];
            gamma[i] += gammai;
            for (j=0; j<ndim; j++) {
                firstmoment[i][j] += gammai * x[t][j];
                secondmoment[i][j] += gammai * x[t][j] * x[t][j];
            }
        }
    }

    return;
}
