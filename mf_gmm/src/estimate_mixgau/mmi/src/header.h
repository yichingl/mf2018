/*#define maxdim 13 */
#define maxdim 40
#define maxmodes 4096
#define PI	3.1415927

#define QUIT(X) {printf("%s(%d):",__FILE__,__LINE__);\
                 printf X;fflush(stdout);exit(0);}

char **alloc2d(int,int,int);
void free2d(void**);
char ***alloc3d(int,int,int,int);
void free3d(void***);

int read_cep(char *file, float***cep, int *numframes, int cepsize);

float gauss(float *x, float *mean, float *hafinvvar, float corprod, int Ndim);

int mgau_mmi(float **x, int Nvec, int Ndim, int Nclass, int *Clslbl,
             int Ngauss, float ***mean, float ***var, float **mixw,
             float *classprior, int niter);

void initialize_mgau(float **x, int nvec, int Ndim,
                     float **mean, float **var, float *c, int Nmodes);

void mgau_em(float **x, int Nvec, int Ndim, int Ngau,
             float **mean, float **var, float *c, int numiters);

int  store_distribution(char *file, int Nmodes, int Ndim,
                        float *mixw, float **mean, float **var);
