/*#define maxdim 13 */

#define QUIT(X) {printf("%s(%d):",__FILE__,__LINE__);\
                 printf X;fflush(stdout);exit(0);}

char **alloc2d(int,int,int);
void free2d(void**);
int read_cep(char *file, float***cep, int *numframes, int cepsize);
int extend_vec(float **vec, int cepdim, int Nvec, int width,
               float ***extendedvec, int *Nfeat);
float  Quantize(int *bins, float **vector, int Nvecs,
                float **codes, int Ncodes, int Ndim);
int  BagofWords(int *bow, float **vector, int Nvecs,
                float **codes, int Ncodes, int Ndim);

