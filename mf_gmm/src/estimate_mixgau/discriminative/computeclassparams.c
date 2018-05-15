#include "discriminative.h"

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
                        float D)
{
    int i, j;
    float summixw, xsqbar, denom;
    // float minmixw, sumthetaclass, sumthetaglobal;

    for (i=0; i<ngauss; i++) {
        denom = gamma[i] + D;
        for (j=0; j<ndim; j++) {
            xsqbar = mean[i][j]*mean[i][j] + var[i][j];

            newmean[i][j] = (firstmoment[i][j] + D*mean[i][j]) / denom;
            newvar[i][j] = (secondmoment[i][j] + D*xsqbar) / denom -
                            newmean[i][j]*newmean[i][j];
        }
    }

    /**
    sumthetaclass = 0;
    sumthetaglobal = 0;
    for (i=0; i<ngauss; i++) {
        sumthetaclass += thetaclass[i];
        sumthetaglobal += thetaglobal[i];
    }

    minmixw = 1e+30;
    for (i=0; i<ngauss; i++) {
        newmixw[i] = thetaclass[i]/sumthetaclass-thetaglobal[i]/sumthetaglobal;
        if (minmixw > newmixw[i]) newmixw[i] = minmixw;
    }

    // Adjust to maintain positive
    if (minmixw < 0) {
        for (i=0; i<ngauss; i++) newmixw[i] += fabs(minmixw) + 1.0;
    }
    summixw = 0;
    for (i=0; i<ngauss; i++) summixw += newmixw[i];
    for (i=0; i<ngauss; i++) newmixw[i] /= summixw;
    **/

/* Just dont!!
    for (i=0, summixw=0; i<ngauss; i++){
        newmixw[i] = thetaglobal[i] > 1e-20 ? 
		mixw[i]*thetaclass[i]/thetaglobal[i]: mixw[i];
	summixw += newmixw[i];
    }
    if (summixw > 0)
        for (i=0; i<ngauss; i++) newmixw[i] /= summixw;
    else
*/
        for (i=0; i<ngauss; i++) newmixw[i] = mixw[i];

    return;
}
