/*
 * PNCC feature extraction function
 *
 * Programmed by Kazumasa Yamamoto, May 9, 2012.
 *
 * Based on Chanwoo Kim's MATLAB code
 * for the IEEETran Speech, Audio, and Langauge Processing and ICASSP 2012
 * (chanwook@cs.cmu.edu)
 */

/*
 * Important: The input should be mono and be sampled at "16 kHz".
 *
 * If you want to use logarithmic nonlinearity instead of the power
 * nonlinearity, change bPowerLaw to 0 (line XX)
 *
 * PNCC_IEEETran(OutFile, InFile)
 *
 * iFiltTyep 2 Gamma
 * iFiltType 1 HTK MEL
 * iFiltType 0 Snaley MEL
 * Default
 * 0.5, 0.01, 2, 4
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "matrix.h"
#include "sigproc.h"
#include "pncc.h"

#define max(x,y) ((x)>(y)?(x):(y))
#define min(x,y) ((x)<(y)?(x):(y))

extern matrix matH; /* Filter response */
extern double dLambda3;

/*
 *
 * Array Queue Ring-buffer
 *
 */
matrix Queue_matP;
int    Queue_iHead;
int    Queue_iTail;
int    Queue_iWindow;
int    Queue_iNumElem;

void Queue_offer(vector vec, int iLen)
{
  int i;

  for (i = 0; i < iLen; i++) {
    Queue_matP[Queue_iTail][i] = vec[i];
  }
    
  Queue_iTail = (Queue_iTail + 1) % Queue_iWindow;
  Queue_iNumElem++;

  if (Queue_iNumElem > Queue_iWindow) {
    fprintf(stderr, "Queue overflow\n");
    exit(EXIT_FAILURE);
  }

  return;
}


// void Queue_poll(vector vec, int iLen)
void Queue_poll()
{
  if (Queue_iNumElem <= 0) {
    fprintf(stderr, "No elements\n");
    exit(EXIT_FAILURE);
  }

#if 0
  {
    int i;
    for (i = 0; i < iLen; i++) {
      vec[i] = Queue_matP[Queue_iHead][i];
    }
  }
#endif
    
  Queue_iHead = (Queue_iHead + 1) % Queue_iWindow;
  Queue_iNumElem--;

  return;
}


void Queue_avg(vector vecMean, int iLen)
{
  int i, j, iPos;

  for (i = 0; i < iLen; i++) {
    iPos = Queue_iHead;
    vecMean[i] = 0.0;
    for (j = 0; j < Queue_iNumElem; j++) {
      vecMean[i] += Queue_matP[iPos][i];
      iPos = (iPos + 1) % Queue_iWindow;
    }
    vecMean[i] = vecMean[i] / Queue_iNumElem;
  }
}


double max_vec(vector vec, int iLen)
{
  int i;
  double max;

  max = vec[0];
  for (i = 1; i < iLen; i++) {
    if (max < vec[i]) max = vec[i];
  }

  return max;
}


void preemphasis(vector vec, int iLen, double dCoef)
{
  int i;

  for (i = iLen - 1; i > 0; i--) {
    vec[i] -= (dCoef * vec[i - 1]);
  }
  // vec[0] *= dCoef; /* bug? */

  return;
}


matrix calc_pncc(vector vecWav, int length, int iNumFilts, int iNumBin, int *iNumFrames)
{
  int i, j, k, m;
  double max;

  double dSampRate    = SAMPRATE;
  double dFrameLen    = FRAMELEN;
  double dFramePeriod = FRAMEPERIOD;

  int iFFTSize = iNumBin * 2;

  double dPowerCoeff = POWERCOEF;

  int    iFiltType       = FILTTYPE;
  double dGammaThreshold = GAMMATHRESH;

  int iSMType = SMTYPE;
  int iM      = M;
  int iN      = N;

  double dLambda   = LAMBDA;
  double dLambda2  = LAMBDA2;
  double dDelta1   = DELTA;
  double dDelta2   = DELTA2;
  double dFactor   = FACTOR;
  double dLambda_L = LAMBDA_L;
  double dLambda_S = LAMBDA_S;

  double dOrg;

  int iFL        = floor(dFrameLen    * dSampRate);
  int iFP        = floor(dFramePeriod * dSampRate);
  int iSpeechLen = length;

  int i_FI     = 0;
  int i_FI_Out = 0;

  double dMeanPower = 1e+10;
  double dMean      = 4e+07; // 4.0638e+07  --> Mean from WSJ0-si84  (Important!!!)
  double dTmp;
    
  /*
   * Flags
   */
  int bPreem    = 1;
  int bSSF      = 1;
  int bPowerLaw = 1;
  int bDisplay  = 0;

  *iNumFrames = floor((length - iFL) / iFP) + 1; // Number of frames

  matrix matP       = new_matrix(*iNumFrames, iNumFilts);
  matrix matX       = new_matrix(*iNumFrames, iNumBin);
  vector vecQ       = new_vector(iNumFilts);
  vector vecQ_Out   = new_vector(iNumFilts);
  vector vecQMVAvg2 = new_vector(iNumFilts);
  vector vecQMVAvg3 = new_vector(iNumFilts);
  vector vecPBias   = new_vector(iNumFilts);
  vector vecW       = new_vector(iNumFilts);
  vector vecW_sm    = new_vector(iNumFilts);
  vector vecX_st_r  = new_vector(iFFTSize);
  vector vecX_st_i  = new_vector(iFFTSize);
  vector vecTmpDCT  = new_vector(iNumFilts);

  if (bSSF == 1) {
    *iNumFrames -= 2 * iM; // The number of frames decreases when using SSF.
  }

  matrix matP_Out    = new_matrix(*iNumFrames, iNumFilts);
  matrix matSpec     = new_matrix(*iNumFrames, iNumFilts);
  matrix matDCT      = new_matrix(*iNumFrames, CEPDIM);
  vector vecSumPower = new_vector(*iNumFrames);

	 
  /*
   * Initialization for Queue
   */
  Queue_iWindow  = 2 * iM + 1;
  Queue_matP     = new_matrix(Queue_iWindow, iNumFilts);
  Queue_iHead    = 0;
  Queue_iTail    = 0;
  Queue_iNumElem = 0;


  /*
   *
   * Pre-emphasis using H(z) = 1 - (PREEMPHASIS) z ^ -1
   *
   */
  if (bPreem == 1) {
    preemphasis(vecWav, length, PREEMPHASIS);
  }


  /*
   *
   * Obtaning the gammatone coefficient. 
   *
   * Based on M. Snelly's auditory toolbox. 
   * In actual C-implementation, we just use a table
   *
   */
  if (iFiltType == 2) {
    // matH[][] has been already read from file until reached here.
    for (i = 0; i < iNumFilts; i++) {
      max = max_vec(matH[i], iNumBin);
      for (j = 0; j < iNumBin; j++) {
	if (matH[i][j] < dGammaThreshold * max) {
	  matH[i][j] = 0.0; // Cut off the coefficient under the threshold
	}
      }
    }
#if 0
    /* for DEBUG */
    for (i = 0; i < iNumFilts; i++) {
      for (j = 0; j < iNumBin; j++) {
	fprintf(stderr, "matH[%d][%d] = %g\n", i, j, matH[i][j]);
	fflush(stderr);
      }
    }
#endif
  }


  /*
   *
   * Obtaining the short-time Power Spectrum X(i, j)
   *
   */
  for (m = 0, i_FI = 0; m <= iSpeechLen - iFL; m += iFP, i_FI++) {

    for (i = m, j = 0; i < m + iFL; i++, j++) {
      vecX_st_r[j] = vecWav[i];
      vecX_st_i[j] = 0.0;
    }
    for (; j < iFFTSize; j++) {
      vecX_st_r[j] = 0.0;
      vecX_st_i[j] = 0.0;
    }
    Hamming(vecX_st_r, iFL);
#if 0
    /* for DEBUG */
    if (i_FI == 0) {
      for (i = 0; i < iFL; i++) {
	fprintf(stderr, "%g ", vecX_st_r[i]);
      }
      fprintf(stderr, "\n");
      fflush(stderr);
    }
#endif
    fft(iFFTSize, vecX_st_r, vecX_st_i);
    for (i = 0; i < iNumBin; i++) {
      // Power Spectrum
      matX[i_FI][i] = vecX_st_r[i] * vecX_st_r[i] + vecX_st_i[i] * vecX_st_i[i];
    }
#if 0
    /* for DEBUG */
    for (i = 0; i < iNumBin; i++) {
      fprintf(stderr, "%d %d %g\n", i_FI, i, log(matX[i_FI][i]));
      fflush(stderr);
    }
    fprintf(stderr, "\n");
#endif
  }


  /*
   *
   * Calculating the Filtered Spectrum P(i, j)
   *
   */
  for (m = 0, i_FI = 0; m <= iSpeechLen - iFL; m += iFP, i_FI++) {
    for (j = 0; j < iNumFilts; j++) {
      /*
       *
       * Squared integration
       *
       */
      if (iFiltType == 2) {
	matP[i_FI][j] = 0.0;
	for (i = 0; i < iNumBin; i++) {
	  matP[i_FI][j] += matX[i_FI][i] * matH[j][i] * matH[j][i]; // use power spectra
	}
      } else {
	matP[i_FI][j] = 0.0;
	for (i = 0; i < iNumBin; i++) {
	  matP[i_FI][j] += matX[i_FI][i] * matH[j][i]; // use power spectra
	}
      }
    }
#if 0
    /* for DEBUG */
    if (i_FI == 20) {
      for (i = 0; i < iNumFilts; i++) {
	fprintf(stderr, "matP[%d][%d] = %g\n", i_FI, i, matP[i_FI][i]);
	fflush(stderr);
      }
    }
#endif
  }

#if 0
  /*
   * Mean subtraction in log spectral domain instead of CMN
   */
  for (i = 0; i < iNumFilts; i++) {
    dTmp = 0.0;
    for (j = 0; j < *iNumFrames; j++) {
      matP[j][i] = log(max(matP[j][i], EPS));
      dTmp += matP[j][i];
    }
    dTmp /= *iNumFrames;
    for (j = 0; j < *iNumFrames; j++) {
      matP[j][i] = exp(matP[j][i] - dTmp);
    }
  }
#endif

  /*
   *
   * Calculating the PNCC P(i, j)
   *
   */
  i_FI     = 0;
  i_FI_Out = 0;
  for (m = 0; m <= iSpeechLen - iFL; m += iFP) {

    if (bSSF == 1) {
      /*
       *
       * Ring buffer (using a Queue)
       *
       */
      if (i_FI >= 2 * iM + 1) {
	Queue_poll();
      }
      Queue_offer(matP[i_FI], iNumFilts);

      Queue_avg(vecQ, iNumFilts);

      if (i_FI == 2 * iM) {
	for (i = 0; i < iNumFilts; i++) {
	  vecPBias[i]  = vecQ[i] * 0.9;
	}
#if 0
	/* for DEBUG */
	for (i = 0; i < iNumFilts; i++) {
	  fprintf(stderr, "vecQ[%d] = %g, vecPBias[%d] = %g\n", i, vecQ[i], i, vecPBias[i]);
	}
	fflush(stderr);
	for (i = 0; i < Queue_iNumElem; i++) {
	  fprintf(stderr, "%d: ", i);
	  for (j = 0; j < iNumFilts; j++) {
	    fprintf(stderr, "%g ", Queue_matP[i][j]);
	  }
	  fprintf(stderr, "\n");
	}
	fflush(stderr);
#endif
      }

      if (i_FI >= 2 * iM) {
	/*
	 *
	 * Bias Update
	 *
	 */
	for (i = 0; i < iNumFilts; i++) {
	  if (vecQ[i] > vecPBias[i]) {
	    vecPBias[i] = dLambda  * vecPBias[i] + (1 - dLambda)  * vecQ[i];
	  } else {
	    vecPBias[i] = dLambda2 * vecPBias[i] + (1 - dLambda2) * vecQ[i];
	  }
	}

#if 0
	/* for DEBUG */
	if (i_FI == 20) {
	  for (i = 0; i < iNumFilts; i++) {
	    fprintf(stderr, "vecQ[%d] = %g, vecPBias[%d] = %g\n", i, vecQ[i], i, vecPBias[i]);
	    fflush(stderr);
	  }
	}
#endif
	for (i = 0; i < iNumFilts; i++) {
	  vecQ_Out[i] =  max(vecQ[i] - vecPBias[i], 0);

	  if (i_FI == 2 * iM) {
	    vecQMVAvg2[i] = 0.9 * vecQ_Out[i];
	    vecQMVAvg3[i] = vecQ_Out[i];
	  }

	  if (vecQ_Out[i] > vecQMVAvg2[i]) {
	    vecQMVAvg2[i] = dLambda  * vecQMVAvg2[i] + (1 - dLambda)  *  vecQ_Out[i];
	  } else {
	    vecQMVAvg2[i] = dLambda2 * vecQMVAvg2[i] + (1 - dLambda2) *  vecQ_Out[i];
	  }

	  dOrg = vecQ_Out[i];

	  vecQMVAvg3[i] *= dLambda3;
			          
#if 0
	  /* for DEBUG */
	  fprintf(stderr, "%d %d %g %g %g %g %g %g",
		  i_FI, i, vecQ[i], dFactor * vecPBias[i], vecQMVAvg2[i],
		  dDelta1 * vecQMVAvg3[i], dDelta2 * vecQMVAvg3[i], vecQ_Out[i]);
#endif
	  if (vecQ[i] < dFactor * vecPBias[i]) {
	    vecQ_Out[i] = vecQMVAvg2[i];
	  } else {
	    if (vecQ_Out[i] <= dDelta1 * vecQMVAvg3[i]) {
	      vecQ_Out[i] = dDelta2 * vecQMVAvg3[i];
	    }
	  }
	  vecQMVAvg3[i] = max(vecQMVAvg3[i], dOrg);

	  vecQ_Out[i] = max(vecQ_Out[i], vecQMVAvg2[i]);
#if 0
	  /* for DEBUG */
	  fprintf(stderr, " %g\n", vecQ_Out[i]);
#endif
	}

	/*
	 *
	 * Spectral smoothing
	 *
	 */
	for (i = 0; i < iNumFilts; i++) {
	  vecW[i] = vecQ_Out[i] / max(vecQ[i], EPS);
	}

	for (i = 0; i < iNumFilts; i++) {
	  if (iSMType == 0) { // using this case, by default
	    vecW_sm[i] = 0.0;
	    for (j = max(i - iN, 0), k = 0; j < min(i + iN + 1, iNumFilts); j++, k++) {
	      vecW_sm[i] += vecW[j];
	    }
	    vecW_sm[i] /= k;
	  } else if (iSMType == 1) {
	    vecW_sm[i] = 0.0;
	    for (j = max(i - iN, 0), k = 0; j < min(i + iN + 1, iNumFilts); j++, k++) {
	      vecW_sm[i] += log(vecW[j]);
	    }
	    vecW_sm[i] /= k;
	    vecW_sm[i] = exp(vecW_sm[i]);
	  } else if (iSMType == 2) {
	    vecW_sm[i] = 0.0;
	    for (j = max(i - iN, 0), k = 0; j < min(i + iN + 1, iNumFilts); j++, k++) {
	      vecW_sm[i] += pow(vecW[j], dPowerCoeff);
	    }
	    vecW_sm[i] /= k;
	    vecW_sm[i] = pow(vecW_sm[i], 1.0 / dPowerCoeff);
	  } else if (iSMType == 3) {
	    vecW_sm[i] = 0.0;
	    for (j = max(i - iN, 0), k = 0; j < min(i + iN + 1, iNumFilts); j++, k++) {
	      vecW_sm[i] += pow(vecW[j], 1.0 / dPowerCoeff);
	    }
	    vecW_sm[i] /= k;
	    vecW_sm[i] = pow(vecW_sm[i], dPowerCoeff);
	  }
	}

	vecSumPower[i_FI_Out] = 0.0;
	for (i = 0; i < iNumFilts; i++) {
	  matP_Out[i_FI_Out][i] =  vecW_sm[i] * matP[i_FI - iM][i];
	  vecSumPower[i_FI_Out] += matP_Out[i_FI_Out][i];
	}

	if (vecSumPower[i_FI_Out] > dMean) {
	  dMean = dLambda_S * dMean + (1 - dLambda_S) * vecSumPower[i_FI_Out];
	} else {
	  dMean = dLambda_L * dMean + (1 - dLambda_L) * vecSumPower[i_FI_Out];
	}

	for (i = 0; i < iNumFilts; i++) {
	  matP_Out[i_FI_Out][i] = matP_Out[i_FI_Out][i] / (dMean)  * dMeanPower;
	}

	i_FI_Out++;
      }
           
    } else { // if not SSF
      vecSumPower[i_FI] = 0.0;
      for (i = 0; i < iNumFilts; i++) {
	vecSumPower[i_FI] += matP[i_FI][i];
      }

      if (vecSumPower[i_FI_Out] > dMean) {
	dMean = dLambda_S * dMean + (1 - dLambda_S) * vecSumPower[i_FI_Out];
      } else {
	dMean = dLambda_L * dMean + (1 - dLambda_L) * vecSumPower[i_FI_Out];
      }

      for (i = 0; i < iNumFilts; i++) {
	matP_Out[i_FI][i] = matP_Out[i_FI][i] / (dMean)  * dMeanPower;
      }
    }

    i_FI++;
  }

  /*
   * Mean subtraction in log spectral domain instead of CMN
   */
  for (i = 0; i < iNumFilts; i++) {
    dTmp = 0.0;
    for (j = 0; j < *iNumFrames; j++) {
      matP_Out[j][i] = log(max(matP_Out[j][i], EPS));
      dTmp += matP_Out[j][i];
    }
    dTmp /= *iNumFrames;
    for (j = 0; j < *iNumFrames; j++) {
      matP_Out[j][i] = exp(matP_Out[j][i] - dTmp);
    }
  }

  /*
   *
   * Apply the nonlinearity (pow() function)
   *
   */
  if (bPowerLaw == 1) {
    for (i = 0; i < *iNumFrames; i++) {
      for (j = 0; j < iNumFilts; j++) {
	matSpec[i][j] = pow(matP_Out[i][j], dPowerCoeff);
      }
    }
  } else {
    for (i = 0; i < *iNumFrames; i++) {
      for (j = 0; j < iNumFilts; j++) {
	matSpec[i][j] = log(matP_Out[i][j] + EPS);
      }
    }
  }
#if 0
  /* for DEBUG */
  for (i = 0; i < *iNumFrames; i++) {
    for (j = 0; j < iNumFilts; j++) {
      fprintf(stderr, "matSpec[%d][%d] = %g\n", i, j, matSpec[i][j]);
      fflush(stderr);
    }
  }
#endif

  /*
   *
   * DCT
   *
   */
  for (i = 0; i < *iNumFrames; i++) { 
    for (j = 0; j < iNumFilts; j++) {
      vecTmpDCT[j] = matSpec[i][j];
    }
    dct(iNumFilts, vecTmpDCT, CEPDIM, matDCT[i]);
  }
#if 0
  /* for DEBUG */
  for (i = 0; i < *iNumFrames; i++) {
    for (j = 0; j < CEPDIM; j++) {
      fprintf(stderr, "matDCT[%d][%d] = %g\n", i, j, matDCT[i][j]);
      fflush(stderr);
    }
  }
#endif

#if 0
  /*
   *
   * CMN
   *
   */
  for (i = 0; i < CEPDIM; i++) {
    dMean = 0.0;
    for (j = 0; j < *iNumFrames; j++) {
      dMean += matDCT[j][i];
    }
    dMean /= *iNumFrames;
    for (j = 0; j < *iNumFrames; j++) {
      matDCT[j][i] -= dMean;
    }
  }
#if 0
  /* for DEBUG */
  for (i = 0; i < *iNumFrames; i++) {
    for (j = 0; j < CEPDIM; j++) {
      fprintf(stderr, "matDCT[%d][%d] = %g\n", i, j, matDCT[i][j]);
      fflush(stderr);
    }
  }
#endif
#endif

  /*
   *
   * Display
   *
   */
  if (bDisplay == 1) {
    /*
      figure

      matSpec = idct(matDCT', iNumFilts);
      imagesc(matSpec); axis xy;
    */
  }

  /*
   * free memory
   */
  free_matrix(matP);
  free_matrix(matP_Out);
  free_matrix(matSpec);
  free_matrix(Queue_matP);
  free_matrix(matX);
  free_vector(vecQ);
  free_vector(vecQ_Out);
  free_vector(vecQMVAvg2);
  free_vector(vecQMVAvg3);
  free_vector(vecPBias);
  free_vector(vecW);
  free_vector(vecW_sm);
  free_vector(vecX_st_r);
  free_vector(vecX_st_i);
  free_vector(vecTmpDCT);
  free_vector(vecSumPower);

  return matDCT;
}
