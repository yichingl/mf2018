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

int  readdata(char *infile, int bracket, int basedim, int maxvecs,
              float ***vector, int *numvecs, int *Ndimfull);
int savefile(char *cbkfile, float **mean, int Ndimfeat, int cepdim, int Nmodes);

void  vector_quantize(float **codes, int Ncodes, float **vector, int Nvecs,
                      int Ndim, int maxiter);
