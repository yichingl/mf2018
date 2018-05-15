#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "matrix.h"
#include "pncc.h"

matrix read_filterbank(char *strFBankFileName, int iNumFilts, int iNumBin);
vector read_speech(char *strInFileName, int *length);
matrix calc_pncc(vector speech, int length, int iNumFilts, int iNumBin, int *iNumFrames);
void write_feature(char *strOutFileName, matrix pncc, int iNumDim, int iNumFrames);

matrix matH;                /* filter bank parameters */
double dLambda3  = LAMBDA3; /* Masking decay parameter */

void usage(char *pname)
{
  fprintf(stderr,
	  "Usage: %s [-m beta]\n"
	  "          -i input-pcm(sph)-file-name\n"
	  "          -o output-sphinx-format-parameter-file-name\n"
	  "          -f filter-bank-parameter-file-name\n"
	  "Options:\n"
	  "          -m beta: Masking decay rate [default: 0.85]",
	  pname);
}


int main(int argc, char *argv[])
{
  char strInFileName[STRLEN], strOutFileName[STRLEN];
  char strFBankFileName[STRLEN];
  int iNumFilts, iNumBin, iNumFrames, iNumDim, length;
  int opt, verbose = 0;

  vector speech;
  matrix pncc;

  iNumFilts = NUMFILTS;
  iNumBin   = (FFTSIZE / 2);
  iNumDim   = CEPDIM;

  strInFileName[0] = strOutFileName[0] = strFBankFileName[0] = '\0';

  while ((opt = getopt(argc, argv, "i:o:f:m:v")) != -1) {
    switch (opt) {
    case 'i':
      strcpy(strInFileName, optarg);
      break;
    case 'o':
      strcpy(strOutFileName, optarg);
      break;
    case 'f':
      strcpy(strFBankFileName, optarg);
      break;
    case 'm':
      dLambda3 = atof(optarg);
      break;
    case 'v':
      verbose = 1;
      break;
    default:
      usage(argv[0]);
      return EXIT_FAILURE;
    }
  }

  if (strInFileName[0] == '\0'
      || strOutFileName[0] == '\0'
      || strFBankFileName[0] == '\0') {
    usage(argv[0]);
    return EXIT_FAILURE;
  }

  if (verbose) {
    fprintf(stderr, "      Input SPEHRE PCM file: %s\n", strInFileName);
    fprintf(stderr, "   Output SPHINX param file: %s\n", strOutFileName);
    fprintf(stderr, "Input FilterBank param file: %s\n", strFBankFileName);
    fprintf(stderr, "         Masking decay rate: %g\n", dLambda3);
  }

  matH = read_filterbank(strFBankFileName, iNumFilts, iNumBin);

#ifdef DEBUG
  {
    int i, j;
    fprintf(stderr, "\n");
    for (i = 0; i < iNumFilts; i++) {
      for (j = 0; j < iNumBin; j++) {
	fprintf(stderr, "matH[%d][%d] = %g\n", i, j, matH[i][j]);
      }
    }
    fprintf(stderr, "\n");
  }
#endif  

  speech = read_speech(strInFileName, &length);

  pncc = calc_pncc(speech, length, iNumFilts, iNumBin, &iNumFrames);

  write_feature(strOutFileName, pncc, iNumDim, iNumFrames);

  return EXIT_SUCCESS;
}
