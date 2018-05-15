#ifndef _TRAINCLASS_H_
#define _TRAINCLASS_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

char **alloc2d(int dim1, int dim2, int size);
void free2d(void **p);

#define IO_ERR	-1
int read_cep (char *file, float ***cep, int *nframes, int numcep);
void read_classes(char *modelfile, int *ngau, int *dim, 
                  float ***mean, float ***vars, float **c);


#define QUIT(X) {printf X;fflush(stdout);exit(0);}
#define WARN(X) {printf("***WARNING****:");printf X;fflush(stdout);}
#define RETURNERR(X) {printf X;fflush(stdout);return(IO_ERR);}

#define ERROR -1
#define SUCCESS 0
#define RETURN_ERROR(X) {printf X; fflush(stdout); return ERROR;} 

#endif 
