/*#define maxdim 13 */
#define maxdim 40
#define maxmodes 4096
#define PI	3.1415927

#define QUIT(X) {printf("%s(%d):",__FILE__,__LINE__);\
                 printf X;fflush(stdout);exit(0);}

char **alloc2d(int,int,int);
int get_featfile_size(char *file, int basedim);
int read_cep(char *file, float***cep, int *numframes, int cepsize);
void free2d(void**);

float **matmulttxpB(float **a, int rowsa, int colsa,
                    float **b, int rowsb, int colsb,
                    float **c);
float **whiteningmatrix(float **cov, int dim, int projecteddim);


int  readdata(char *infile, int bracket, int basedim, int maxvecs,
              float ***vector, int *numvecs, char *covarfile, int featdim,
              float ***whiteningmat);

int savefile(char *cbkfile, float **mean, float **whiteningmat,
             int Ndimfeat, int Ndimfull, int cepdim, int Nmodes);

void  vector_quantize(float **codes, int Ncodes, float **vector, int Nvecs,
                      int Ndim, int maxiter);
