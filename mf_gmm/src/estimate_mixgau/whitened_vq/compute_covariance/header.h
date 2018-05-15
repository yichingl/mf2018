/*#define maxdim 13 */
#define maxdim 40
#define maxmodes 4096
#define PI	3.1415927

#define QUIT(X) {printf("%s(%d):",__FILE__,__LINE__);\
                 printf X;fflush(stdout);exit(0);}

char **alloc2d(int,int,int);
int read_cep(char *file, float***cep, int *numframes, int cepsize);
void free2d(void**);
double gaussrand(void);
int matrixinvert(float **ainv, float **a, int len);

void  vector_quantize(float **codes, int Ncodes, float **vector, int Nvecs,
                      int Ndim, int maxiter);
