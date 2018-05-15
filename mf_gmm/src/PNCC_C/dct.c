#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <matrix.h>


int dct(int n, double x[], int k, double y[])
{
  int i, j;
  double dB, dC;

  static int    last_n = 0;
  static double dA = 0.0;
  static double sqrt2 = 1.41421356;
  static double **costbl = NULL;

  if (n != last_n || n == 0) {
    last_n = n;
    if (costbl != NULL) free_matrix(costbl);
    if (n == 0) return 0;  /* free the memory of the DCT matrix */

    sqrt2 = sqrt(2.0);
    dA    = sqrt(2.0 / (double)n);
    dB    = M_PI / (double)n;

    costbl = new_matrix(n, n);
    for (i = 0; i < n; i++) {
      dC = i * dB;
      for (j = 0; j < n; j++) {
	costbl[i][j] = cos(dC * (j + 0.5));
      }
    }
  }


  // DCT
  for (i = 0; i < k; i++) {
    y[i] = 0.0;
    for (j = 0; j < n; j++) {
      y[i] += x[j] * costbl[i][j];
    }
    y[i] *= dA;
  }
  y[0] /= sqrt2;

  return 0;  /* finished successfully */
}
