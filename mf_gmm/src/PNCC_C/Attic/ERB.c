/*
 * Test program of ERB spaced frequency
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main()
{
  double lowFreq = 200;
  double highFreq = 8000;
  int N = 40;

  double EarQ  = 9.26449; // Glasberg and Moore Parameters
  double minBW = 24.7;
  int    order = 1;

  double coeff;
  int i;

  for (i = 0; i < N; i++) {
    coeff = -(EarQ * minBW)
      + exp((i + 1) * (-log(highFreq + EarQ * minBW) + log(lowFreq + EarQ * minBW)) / N)
      * (highFreq + EarQ * minBW);

    printf("%f\n", coeff);
  }

  return 0;
}
