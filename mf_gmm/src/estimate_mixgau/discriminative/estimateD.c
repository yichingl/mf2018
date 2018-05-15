#include "discriminative.h"

float estimateD(float **mean,
                float **var,
                float **firstmoment,
                float **secondmoment,
                float *gamma,
                int   ngauss,
                int   ndim)
{
    int i, j;
    float bstD, A, B, C, D, sqterm;

    bstD = -1e+30;
    for (i=0; i<ngauss; i++) {
        for (j=0; j<ndim; j++) {
            A = var[i][j];
            B = gamma[i]*(var[i][j] + mean[i][j]*mean[i][j]) +
                secondmoment[i][j] - 2.0*firstmoment[i][j]*mean[i][j];
            C = gamma[i]*secondmoment[i][j] -
                firstmoment[i][j] * firstmoment[i][j];

            sqterm = B*B - 4*A*C;
            if (sqterm < 0) // Any D will do since vars are +ve for all D
                D = 1.0;
            else {
                sqterm = sqrt(sqterm);
                D = (-B + sqterm)/2*A + 1.0; // The larger solution
                if (D < 1.0) D = 1.0; // Ensure postiveness
            }
            if (D > bstD) bstD = D;
        }
    }

    return 2.0*bstD;  // Increase D for safety
}
