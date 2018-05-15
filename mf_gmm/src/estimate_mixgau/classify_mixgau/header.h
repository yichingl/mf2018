#ifndef _AYESHA_H_
#define _AYESHA_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>


void lloyd_quantize(double *codes,int Ncodes,double *vector,int Nvecs);

char **alloc2d(int dim1, int dim2, int size);
void free2d(void **p);


/* FOR BOTH THE FOLLOWING numcep -HAS- TO BE KNOWN AND PASSED
   read_cep WILL RETURN numframes. write_cep ALSO NEEDS nframes */
#define IO_ERR	-1
int read_cep (char *file, float ***cep, int *nframes, int numcep);


#define QUIT(X) {printf("%s(%d):",__FILE__,__LINE__);\
                 printf X;fflush(stdout);exit(0);}


#endif
