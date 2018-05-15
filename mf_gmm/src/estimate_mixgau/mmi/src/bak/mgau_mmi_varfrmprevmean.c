#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "header.h"


/*--------------------------------------------------------------------------------
 * The log of a Gaussian 
 *------------------------------------------------------------------------------*/

float computeloggauss(float *x, float *m, float *iv, float cp, int dim)
{
    int i;
    float val = cp;

    for (i=0; i<dim; i++) {
        float xm = x[i] - m[i];
        val -= xm * xm * iv[i];
    }

    return val;
}


/*---------------------------------------------------------------------------
 * Computing the D term
 *--------------------------------------------------------------------------*/

float compute_D(float **nmN, float **nmD, float **nvN, float **nvD,
                float *denN, float *denD, float **mean, float **var, 
                int Ngauss, int Ndim)
{
    int gaus, i;
    float D, Di, Dflr;
    const float shrinkage = 0.9; // Never shrink var by more than this
    const float invshrinkage = shrinkage/(1-shrinkage);

    // Solve a quadratic for D for each dim. Pick the highest D.
    D = 0; // Default error flag value

    for (gaus = 0; gaus < Ngauss; gaus++) {
        for (i=0; i<Ndim; i++) {
            float thetav = nvN[gaus][i] - nvD[gaus][i];
            float gamma = denN[gaus] - denD[gaus];
  
            //if (gamma < 0) Dflr = 1-gamma;
            Dflr = denD[gaus];

            Di = invshrinkage * gamma - thetav/((1.0-shrinkage)*var[gaus][i]);
            if (Di < Dflr) Di = Dflr;
            if (Di > D) D = Di;
        }
    }
    D += 1.0; // A tiny safety buffer

    return D;
}

/*--------------------------------------------------------------------------------
 * One iteration of an MMI Estimator
 * ------------------------------------------------------------------------------*/

float iterGMM_MMI(float **x, int Nvec, int Ndim, int Nclass, int *Clslbl, 
                  int Ngauss, float ***mean, float ***var, float **mixw, 
                  float *classprior)
{
    float **Tau, *Tauc;
    float **corprod, ***halfinvvar;
    float ***newmeanN, ***newmeanD, ***newvarN, ***newvarD;
    float **denN, **denD, totallogposterior;
    int t, i, j, k;
    int cls, gaus;

    Tau = (float**) alloc2d(Nclass, Ngauss, sizeof(float));
    corprod = (float**) alloc2d(Nclass, Ngauss, sizeof(float));
    denN = (float**) alloc2d(Nclass, Ngauss, sizeof(float));
    denD = (float**) alloc2d(Nclass, Ngauss, sizeof(float));
    halfinvvar = (float***) alloc3d(Nclass, Ngauss, Ndim, sizeof(float));
    newmeanN = (float***) alloc3d(Nclass, Ngauss, Ndim, sizeof(float));
    newmeanD = (float***) alloc3d(Nclass, Ngauss, Ndim, sizeof(float));
    newvarN = (float***) alloc3d(Nclass, Ngauss, Ndim, sizeof(float));
    newvarD = (float***) alloc3d(Nclass, Ngauss, Ndim, sizeof(float));
    Tauc = (float*) calloc(Ngauss, sizeof(float));

    // Compute corpord
    for (i=0; i<Nclass; i++) {
        for (j=0; j<Ngauss; j++) {
            corprod[i][j] = log(classprior[i]) + log(mixw[i][j]);
            for (k=0; k<Ndim; k++) {
                corprod[i][j] -= 0.5*log(var[i][j][k]);
                halfinvvar[i][j][k] = 0.5/var[i][j][k];
            }
        }
    }

    totallogposterior = 0.0;
    for (t = 0; t < Nvec; t++) {
        int classlabel = Clslbl[t];
        float **m = mean[0];
        float **iv = halfinvvar[0];
        float *cp = corprod[0];

        // First compute the log gaussians for all gaussians of all classes
        float sumclasstau, sumtau;
        float maxclslblgau = 1e-33, maxgau = -1e+33;

        for (cls = 0; cls < Nclass; cls++) {
            float maxclassgau = -1e+33;
            for (gaus = 0; gaus < Ngauss; gaus++) {
                Tau[cls][gaus] = computeloggauss(x[t],*m++,*iv++,*cp++,Ndim);
                if (Tau[cls][gaus] > maxclassgau) maxclassgau = Tau[cls][gaus];
            }
            if (maxgau < maxclassgau) maxgau = maxclassgau;
            if (cls == classlabel) maxclslblgau = maxclassgau;
        }

        // Now compute class-specific posterior
        sumclasstau = 0;
        for (gaus = 0; gaus < Ngauss; gaus++) {
            Tauc[gaus] = exp(Tau[classlabel][gaus] - maxclslblgau);
            sumclasstau += Tauc[gaus];
        }
        for (gaus = 0; gaus < Ngauss; gaus++) Tauc[gaus] /= sumclasstau;

        // Compute posteriors for all gaussians of all classes
        sumtau = 0;
        for (cls = 0; cls < Nclass; cls++) {
            for (gaus = 0; gaus < Ngauss; gaus++) {
                Tau[cls][gaus] = exp(Tau[cls][gaus] - maxgau);
                sumtau += Tau[cls][gaus];
            }
        }
        for (cls = 0; cls < Nclass; cls++) {
            for (gaus = 0; gaus < Ngauss; gaus++) {
                Tau[cls][gaus] /= sumtau;
            }
        }

        totallogposterior += log(sumclasstau) + maxclslblgau - log(sumtau) - maxgau;

        // Now we're set. Just update buffers
        // First the self class
        for (gaus = 0; gaus < Ngauss; gaus++) {
            float gaucontrib = Tauc[gaus];
            for (j = 0; j < Ndim; j++) {
                float xadj = x[t][j] - mean[classlabel][gaus][j];
                newmeanN[classlabel][gaus][j] += gaucontrib * x[t][j];
                newvarN[classlabel][gaus][j] += gaucontrib * xadj * xadj;
            }
            denN[classlabel][gaus] += gaucontrib;
        }

        // The discriminative component
        for (cls = 0; cls < Nclass; cls++) {
            for (gaus = 0; gaus < Ngauss; gaus++) {
                float gaucontrib = Tau[cls][gaus];
                for (j = 0; j < Ndim; j++) {
                    float xadj = x[t][j] - mean[cls][gaus][j];
                    newmeanD[cls][gaus][j] += gaucontrib * x[t][j];
                    newvarD[cls][gaus][j] += gaucontrib * xadj * xadj;
                }
                denD[cls][gaus] += gaucontrib;
            }
        }
    }

    // Now update
    // For mixw we use the formula:
    // newmw = curmw(gammanum/sum(gammanum) - gammaden/sum(gammaden) + C)*cnst
    // Note: Valtchev 97.
    for (cls = 0; cls < Nclass; cls++) {
        float **nmN = newmeanN[cls];
        float **nmD = newmeanD[cls];
        float **nvN = newvarN[cls];
        float **nvD = newvarD[cls];
        float *dNcls = denN[cls];
        float *dDcls = denD[cls];
        float **mcls = mean[cls];
        float **vcls = var[cls];
        float *mwbuf = (float *) calloc(Ngauss, sizeof(float));
        float C, D, summw, minmixw = 1e-10;
        float sumdN = 0;
        float sumdD = 0;

        D = compute_D(nmN,nmD,nvN,nvD,dNcls,dDcls,mcls,vcls,Ngauss,Ndim);
        for (gaus = 0; gaus < Ngauss; gaus++) {
            float den, m, v;

            sumdN += dNcls[gaus]; sumdD += dDcls[gaus];
            den = dNcls[gaus] - dDcls[gaus] + D;

            for (j = 0; j < Ndim; j++) {
                m = mean[cls][gaus][j];
                v = var[cls][gaus][j];
                mean[cls][gaus][j] = (nmN[gaus][j] - nmD[gaus][j] + D*m) / den;
                var[cls][gaus][j] = (nvN[gaus][j] - nvD[gaus][j] + D*v) / den;
            }
        }

        // The mixture weights
        /*
        C = 0;
        for (gaus = 0; gaus < Ngauss; gaus++) {
            mwbuf[gaus] = dNcls[gaus]/sumdN - dDcls[gaus]/sumdD;
            if (mwbuf[gaus]+C < minmixw) C = minmixw - mwbuf[gaus];
        }

        summw = 0; 
        for (gaus = 0; gaus < Ngauss; gaus++) {
            mixw[cls][gaus] = (mwbuf[gaus] + C) * mixw[cls][gaus];
            summw += mixw[cls][gaus];
        }
        for (gaus = 0; gaus < Ngauss; gaus++) mixw[cls][gaus] /= summw;
        */
    }

    free2d((void**)Tau); free(Tauc);
    free2d((void**)corprod);
    free2d((void**)denN);
    free2d((void**)denD);
    free3d((void***)halfinvvar);
    free3d((void***)newmeanN);
    free3d((void***)newmeanD);
    free3d((void***)newvarN);
    free3d((void***)newvarD);

    return totallogposterior;
}


int mgau_mmi(float **x, int Nvec, int Ndim, int Nclass, int *Clslbl,
             int Ngauss, float ***mean, float ***var, float **mixw,
             float *classprior, int niter)
{
    int iter;
    float score, prevscore, improvement;

    prevscore = iterGMM_MMI(x, Nvec, Ndim, Nclass, Clslbl, Ngauss,
                            mean, var, mixw, classprior);
    printf("Iter 1: Log posterior = %f\n",prevscore);
    for (iter = 2; iter <= niter; iter++) {
        score = iterGMM_MMI(x, Nvec, Ndim, Nclass, Clslbl, Ngauss,
                            mean, var, mixw, classprior);
        improvement = 100.0*(score-prevscore)/fabs(prevscore);
        printf("Iter %d: Log posterior = %f, improvement = %f%%\n",
               iter, score, improvement);
        prevscore = score;
    }

    return 0;
}
