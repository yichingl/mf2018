#include "discriminative.h"

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
                       float  *betweenclassposterior)
{
    int i, j, k;
    float **logscore, tmp, max, maxinclass, totallkd, totalbtwclasslkd;

    // Find log scores for all Gaussians of all classes. 
    logscore = (float **) alloc2d(nclass,ngauss,sizeof(float));
    max = -1e+30; 
    for (i=0; i<nclass; i++){
        for (j=0; j<ngauss; j++){
            logscore[i][j] = corprod[i][j] + logclassprior[i];
            for (k=0; k<ndim; k++) {
                tmp = x[k] - mean[i][j][k];
                logscore[i][j] -= tmp*tmp*halfinvvar[i][j][k];
            }
            if (max < logscore[i][j]) max = logscore[i][j];
        }
    }

    // Find within class a posteriors
    maxinclass = logscore[classid][0];
    for (j=1; j<ngauss; j++)
        if (maxinclass < logscore[classid][j])
            maxinclass = logscore[classid][j];

    totallkd = 0;
    for (j=0; j<ngauss; j++) {
        withinclassposterior[j] = exp(logscore[classid][j] - maxinclass);
        totallkd += withinclassposterior[j];
    }
    for (j=0; j<ngauss; j++) 
        withinclassposterior[j] /= totallkd;

    // Find between class a posteriors
    totalbtwclasslkd = 0;
    for (i=0; i<nclass; i++) {
        for (j=0; j<ngauss; j++) {
            logscore[i][j] = exp(logscore[i][j] - max);
            totalbtwclasslkd += logscore[i][j];
        }
    }
    for (j=0; j<ngauss; j++) 
        betweenclassposterior[j] = logscore[classid][j]/totalbtwclasslkd;

    free2d((void**)logscore);
}
