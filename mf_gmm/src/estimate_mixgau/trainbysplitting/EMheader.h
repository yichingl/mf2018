#ifndef header_h
#define header_h

#include <stdio.h>
#include <math.h>

#ifndef REAL              // Allows dynamic definition
#define REAL   float
#endif
#ifndef NUMITERS          // Allows dynamic definition
#define NUMITERS 10
#endif

int  read_distribution(char *file, int *Nmodes, int *Ndim,
                       REAL **c, REAL ***mean, REAL ***var);
int  store_distribution(char *file, int Nmodes, int Ndim,
                        REAL *c, REAL **mean, REAL **var);
void estimate_mixture_gaussian(REAL **x,       /* The observation vectors */
                                int N,         /* Number of obs. vectors */
                                int Ndim,      /* Dimensions in obs. vectors */
                                int K,         /* Number of Gaussians needed */
                               REAL **mean,    /* The Gaussian means */
                               REAL **var,     /* The Gaussian covarinances */
                               REAL *c,        /* Mixture weights */
                               char *tempfile, /* Temporary file name */
                                int inputNgau  /* If we're continuing from a
                                                  previous distribution. If
                                                  we should begin from scratch,
                                                  set this to 0 */
                              );



#define PI	3.1415927

#define QUIT(X) {printf("%s(%d):",__FILE__,__LINE__);\
                 printf X;fflush(stdout);exit(0);}

char          **alloc2d(int dim1,	/* "x" dimension */ 
			int dim2,	/* "y" dimension */
			int size	/* number of bytes each entry takes */
			);
void free2d(void **p);
int read_cep(char *file, float***cep, int *numframes, int cepsize);
int  write_cep(char *cepfilename, float **cep, int nframes, int numcep);

#endif
