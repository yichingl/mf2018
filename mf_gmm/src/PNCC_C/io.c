#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <matrix.h>


matrix read_filterbank(char *fname, int iNumFilts, int iNumBin)
{
  FILE *fp;
  int i, j;

  matrix H = new_matrix(iNumFilts, iNumBin);

  if ((fp = fopen(fname, "r")) == NULL) {
    perror("read_filterbank");
    exit(EXIT_FAILURE);
  }

  for (i = 0; i < iNumFilts; i++) {
    for (j = 0; j < iNumBin; j++) {
      if (fread(&(H[i][j]), sizeof(double), 1, fp) != 1) {
	fprintf(stderr,"read_filterbank: read error %s\n", fname);
	exit(EXIT_FAILURE);
      }
    }
  }
  fclose(fp);

  return H;
}


vector read_speech(char *fname, int *iSamp)
{
  struct stat buf;
  FILE *fp;
  vector speech;
  int length, i;
  int bytestoskip = 0;
  short tmp;

  if (stat(fname, &buf) != 0) {
    perror("read_speech");
    exit(EXIT_FAILURE);
  }

  if (strcmp(fname+strlen(fname)-3,"wav")==0)  bytestoskip = 44;
  else if (strcmp(fname+strlen(fname)-3,"sph")==0)  bytestoskip = 1024;

  length = (buf.st_size - bytestoskip) / sizeof(short);
  speech = new_vector(length);

  if ((fp = fopen(fname, "r")) == NULL) {
    perror("read_speech");
    exit(EXIT_FAILURE);
  }

  fseek(fp, bytestoskip, SEEK_SET);

  for (i = 0; i < length; i++) {
    if (fread(&tmp, sizeof(short), 1, fp) != 1) {
      fprintf(stderr,"read_speech: read error %s\n", fname);
      exit(EXIT_FAILURE);
    }
    speech[i] = tmp;
  }
  fclose(fp);

  *iSamp = length;
  return speech;
}


void write_feature(char *fname, matrix feat, int iNumDim, int iNumFrame)
{
  FILE *fp;
  int i, j, itmp;
  float ftmp;

  if ((fp = fopen(fname, "w")) == NULL) {
    perror("write_feature");
    exit(EXIT_FAILURE);
  }

  itmp = iNumDim * iNumFrame;
  if (fwrite(&itmp, sizeof(int), 1, fp) != 1) {
    fprintf(stderr,"write_feat: write error %s\n", fname);
    exit(EXIT_FAILURE);
  }

  for (i = 0; i < iNumFrame; i++) {
    for (j = 0; j < iNumDim; j++) {

      ftmp = (float)feat[i][j];
#if 0
      /* for DEBUG */
      fprintf(stderr, "feat[%d][%d] = %g\n", i, j, ftmp);
      fflush(stderr);
#endif

      if (fwrite(&ftmp, sizeof(float), 1, fp) != 1) {
	fprintf(stderr,"write_feat: write error %s\n", fname);
	exit(EXIT_FAILURE);
      }
    }
  }

  fclose(fp);

  return;
}
