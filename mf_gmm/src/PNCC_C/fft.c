/*
 * FFT from "Algorithm Dictionary in C" (in Japanese)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


static void make_sintbl(int n, double sintbl[])
{
  int i, n2, n4, n8;
  double c, s, dc, ds, t;

  n2 = n / 2;  n4 = n / 4;  n8 = n / 8;
  t = sin(M_PI / n);
  dc = 2 * t * t;  ds = sqrt(dc * (2 - dc));
  t = 2 * dc;  c = sintbl[n4] = 1;  s = sintbl[0] = 0;
  for (i = 1; i < n8; i++) {
    c -= dc;  dc += t * c;
    s += ds;  ds -= t * s;
    sintbl[i] = s;  sintbl[n4 - i] = c;
  }
  if (n8 != 0) sintbl[n8] = sqrt(0.5);
  for (i = 0; i < n4; i++)
    sintbl[n2 - i] = sintbl[i];
  for (i = 0; i < n2 + n4; i++)
    sintbl[i + n2] = - sintbl[i];
}


static void make_bitrev(int n, int bitrev[])
{
  int i, j, k, n2;

  n2 = n / 2;  i = j = 0;
  for ( ; ; ) {
    bitrev[i] = j;
    if (++i >= n) break;
    k = n2;
    while (k <= j) { j -= k; k /= 2; }
    j += k;
  }
}


int fft(int n, double x[], double y[])
{
  static int    last_n = 0;     /* previous n */
  static int    *bitrev = NULL; /* bit reversal table */
  static double *sintbl = NULL; /* trigonometric function table */
  int i, j, k, ik, h, d, k2, n4, inverse;
  double t, s, c, dx, dy;

  /* preparation */
  if (n < 0) {
    n = -n;  inverse = 1; /* inverse transform */
  } else     inverse = 0;
  n4 = n / 4;
  if (n != last_n || n == 0) {
    last_n = n;
    if (sintbl != NULL) free(sintbl);
    if (bitrev != NULL) free(bitrev);
    if (n == 0) return 0;             /* free the memory */
    sintbl = malloc((n + n4) * sizeof(double));
    bitrev = malloc(n * sizeof(int));
    if (sintbl == NULL || bitrev == NULL) {
      fprintf(stderr, "Memory shortage\n");  return 1;
    }
    make_sintbl(n, sintbl);
    make_bitrev(n, bitrev);
  }

  /* bit reversal */
  for (i = 0; i < n; i++) {
    j = bitrev[i];
    if (i < j) {
      t = x[i];  x[i] = x[j];  x[j] = t;
      t = y[i];  y[i] = y[j];  y[j] = t;
    }
  }

  /* transformation */
  for (k = 1; k < n; k = k2) {
    h = 0;  k2 = k + k;  d = n / k2;
    for (j = 0; j < k; j++) {
      c = sintbl[h + n4];
      if (inverse) s = - sintbl[h];
      else         s =   sintbl[h];
      for (i = j; i < n; i += k2) {
	ik = i + k;
	dx = s * y[ik] + c * x[ik];
	dy = c * y[ik] - s * x[ik];
	x[ik] = x[i] - dx;  x[i] += dx;
	y[ik] = y[i] - dy;  y[i] += dy;
      }
      h += d;
    }
  }
  if (inverse) /* divide by n in case of the inverse transformation */
    for (i = 0; i < n; i++) { x[i] /= n;  y[i] /= n; }
  return 0;  /* finished successfully */
}



/* Test code */

#if 0

#define N 64

int main(void)
{
  int i;
  static double x1[N], y1[N], x2[N], y2[N], x3[N], y3[N];

  for (i = 0; i < N; i++) {
    x1[i] = x2[i] = 6 * cos( 6 * M_PI * i / N)
                  + 4 * sin(18 * M_PI * i / N);
    y1[i] = y2[i] = 0;
  }
  if (fft(N, x2, y2)) return EXIT_FAILURE;
  for (i = 0; i < N; i++) {
    x3[i] = x2[i];  y3[i] = y2[i];
  }
  if (fft(-N, x3, y3)) return EXIT_FAILURE;
  printf("      Org data      FFT           IFFT\n");
  for (i = 0; i < N; i++)
    printf("%4d | %6.3f %6.3f |"
	   " %6.3f %6.3f | %6.3f %6.3f\n",
	   i, x1[i], y1[i], x2[i], y2[i], x3[i], y3[i]);
  return EXIT_SUCCESS;
}

#endif
