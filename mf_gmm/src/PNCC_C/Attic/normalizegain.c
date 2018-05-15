/*
 * Normalize the Gain of the filter
 *
 * Programmed by Kazumasa Yamamoto, May 08, 2012.
 *
 * Based on Chanwoo Kim's MATLAB code
 * for the IEEETran Speech, Audio, and Langauge Processing and ICASSP 2012
 * (chanwook@cs.cmu.edu)
 */

#include <math.h>
#include <matrix.h>

void NormalizeFilterGain(matrix filter, int iNumFilts, int iFFTSize, double dSampRate)
{
  int i, j;
  double dWeight;

  for (i = 0; i < iNumFilts; i++) {
    dWeight = 0.0;
    for (j = 0; j < iFFTSize / 2; j++) {
      dWeight += fabs(filter[i][j] * filter[i][j]);
    }
    dWeight = sqrt(dWeight * dSampRate / iFFTSize);
    filter[i][j] /= dWeight;
  }

  return;
}

/*
function [aad_H] = NormalizeGain(aad_H, dSampRate)


        [M, N] = size(aad_H);


    iFFTSize  = size(aad_H, 1) * 2;
        adWeight = sqrt(sum(abs(aad_H .* aad_H)) * (dSampRate / iFFTSize));




    % adWeight = adWeight / adWeight(1);


        for i = 1 : N

                aad_H(:, i) = aad_H(:, i) / adWeight(i);

        end

%       plot(linspace(0, 8000, length(aad_H)), 20 * log10(abs(aad_H)));




end
*/
