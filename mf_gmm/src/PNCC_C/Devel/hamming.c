#include <stdlib.h>
#include <math.h>

int chkws = 0;

int Hamming(double *x, int winsize)
{
  int i;
  static double *window = NULL;

  if (chkws != winsize) {
    chkws = winsize;
    window = (double*)realloc(window, winsize * sizeof(double));
    if (window == NULL) {
      return 1;
    }
    for (i = 0; i < winsize; i++) {
      window[i] = 0.54 - 0.46 * cos(2.0 * M_PI * i / (winsize - 1.0));
    }
  }

  for (i = 0; i < winsize; i++) {
    x[i] *= window[i];
  }

  return 0;
}
