#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "header.h"

#define DESCRIPTION "\nThis program computes a VQ codebook from a set of feature vectors.  The data are assumed to be in mfc format. The code will create extended vectors from K consecutive vectors prior to training the codebook.\n"

#define USAGE "USAGE:\n%s -i <input control file> \\\n\t-d <dimensionality of input> \\\n\t-w <bracket -- extended vec = currentframe +- w frames>\\\n\t-m <Max vectors to use in training> \\\n\t-o <output codebook file> \\\n\t-c <No of codewords in codebook>\\\n\t-r <Covariance matrix for whitening (def NULL)> \\\n\t-D <reduced dimensions from projection (def fulldim)> \\\n\t-n <Max iterations [Default 15]>\n"


int main(int argc, char **argv)
{
    float   **vector, **mean, **whiteningmat=NULL;
    int     i, numvecs, cepdim, Ndimfull, Ndimfeat, bracket; 
    int     Nmodes, maxiter, maxvecs;
    char    *infile, *cbkfile=NULL, *covarfile=NULL;
    clock_t start, end;
    double  cpu_time_used;

    maxiter = -1;
    Nmodes = 0; maxvecs = -1; Ndimfeat = -1;
    /* Print out command, just to be nice */
    for (i=0; i<argc; i++) printf("%s ",argv[i]); printf("\n"); fflush(stdout);
    if (argc == 1) {
        printf(DESCRIPTION);
        QUIT ((USAGE,argv[0]));
    }
    for (i=1;i<argc;++i) {
        char sw;
        if (argv[i][0]!='-')
            QUIT ((USAGE,argv[0]));
        sw = argv[i][1];
        switch(sw) {
            case 'i': {
                infile = argv[++i];
                break;
            }
            case 'r': {
                covarfile = argv[++i];
                break;
            }
            case 'w': case 'W': {
                bracket = atoi(argv[++i]);
                break;
            }
            case 'd': {
                cepdim = atoi(argv[++i]);
                break;
            }
            case 'D': {
                Ndimfeat = atoi(argv[++i]);
                break;
            }
            case 'm': case 'M': {
                maxvecs = atoi(argv[++i]);
                break;
            }
            case 'o': case 'O': {
                cbkfile = argv[++i];
                break;
            }
            case 'c': case 'C': {
                Nmodes = atoi(argv[++i]);
                break;
            }
            case 'n': case 'N': {
                maxiter = atoi(argv[++i]);
                break;
            }
            default: {
                printf(DESCRIPTION);
                QUIT ((USAGE,argv[0]));
            }
        }
    }
    if (!cbkfile) QUIT(("No output file. Argument -o must be specified\n"));

    Ndimfull = cepdim*(2*bracket + 1); 
    if (Ndimfeat <= 0) Ndimfeat = Ndimfull; // Full

    start = clock(); 
    readdata(infile, bracket, cepdim, maxvecs, &vector, &numvecs, 
             covarfile, Ndimfeat, &whiteningmat);
    end = clock();
    cpu_time_used = ((double)(end-start)) / CLOCKS_PER_SEC;
    printf("%f secs to read %d vectors.\n",cpu_time_used,numvecs); fflush(stdout);


    mean = (float **) alloc2d(Nmodes,Ndimfeat,sizeof(float));

    // Save whitening matrix, basically
    savefile(cbkfile, mean, whiteningmat, Ndimfeat, Ndimfull, cepdim, Nmodes);
    vector_quantize(mean,Nmodes,vector,numvecs,Ndimfeat,maxiter);

    savefile(cbkfile, mean, whiteningmat, Ndimfeat, Ndimfull, cepdim, Nmodes);

    free2d((void**)vector); free2d((void**)mean);

    return 0;
}
