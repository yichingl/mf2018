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

float **matmult(float **a, int rowsa, int colsa,
                float **b, int rowsb, int colsb,
                float **c);
float **matmulttxpA(float **a, int rowsa, int colsa,
                    float **b, int rowsb, int colsb,
                    float **c);
float **matmulttxpB(float **a, int rowsa, int colsa,
                    float **b, int rowsb, int colsb,
                    float **c);

void  bagofwords(float **mean, int Nmodes, float **vector,
                 int numvecs, int Ndim, int *counts, float *projs);

