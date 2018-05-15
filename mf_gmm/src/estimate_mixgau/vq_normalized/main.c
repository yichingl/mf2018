#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gsl/gsl_cblas.h>
#include "header.h"

#define DESCRIPTION "\nThis program computes a VQ codebook from a set of feature vectors.  The data are assumed to be in mfc format. The code will create extended vectors from K consecutive vectors prior to training the codebook.\n"

#define USAGE "USAGE:\n%s -i <input control file> \\\n\t-d <dimensionality of input> \\\n\t-w <bracket -- extended vec = currentframe +- w frames>\\\n\t-o <output codebook file> \\\n\t-c <No of codewords in codebook>\\\n\t-n <Max iterations [Default 100]>\n"


float **matmult(float **a, float **b, int rowsa, int colsa, int colsb)
{
    float **c = (float **) alloc2d(colsb, rowsa, sizeof(float));

    cblas_sgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
                rowsa, colsb, colsa,
                1, a[0], rowsa, b[0], colsa, 0, c[0], rowsa);
    return c;
}

void mean_normalize(float **x, int dim, int nvec)
{
    int i, j;
    float mean[1024];

    for (j=0; j<dim; j++) mean[j] = 0;

    for (i=0; i<nvec; i++) 
        for (j=0; j<dim; j++) mean[j] += x[i][j];

    for (j=0; j<dim; j++) mean[j] /= (float) nvec;

    for (i=0; i<nvec; i++) 
        for (j=0; j<dim; j++) x[i][j] -= mean[j];

    return; 
}

float **extend_vector(float **x, int ndim, int nvecs, int bracket, int *onvecs)
{
    int i, j, k, l, m, fulldim;
    float **y;

    fulldim = ndim * (2*bracket + 1);
    *onvecs = nvecs - 2*bracket;
    if (*onvecs < 0) {*onvecs = 0; return NULL;}

    y  = (float **) alloc2d(*onvecs, fulldim, sizeof(float));
    if (!y) QUIT(("Unable to allocate %d x %d extended vecs\n",*onvecs,fulldim))

    for (i=bracket, m=0; i<nvecs-bracket; i++,m++) {
        j = 0;
        for (k=i-bracket; k<=i+bracket; k++) {
            for (l=0; l<ndim; l++) {
                y[m][j++] = x[k][l];
            }
        }
    }

    return y;
}


/* Read input features, create extended feature, whiten. */
int  readdata(char *infile, int bracket, 
              float ***vector, int *numvecs, int *Ndim, float ***invcovariance)
{
    int totalframes, nfr, extnfr, fulldim, basedim, i, j, k, nfiles;
    char cepfile[1024];
    float **cepfeat, **extcepfeat, **feat, **covar, **invcovar;
    FILE *fp;

    fprintf(stdout,"Reading input.\n");
    basedim = *Ndim;
    *Ndim = fulldim = (2*bracket + 1) * basedim;

    if ((fp = fopen(infile,"r")) == NULL) QUIT(("Cannot open %s\n",infile))

    covar = (float **) alloc2d(fulldim, fulldim, sizeof(float));
    invcovar = (float **) alloc2d(fulldim, fulldim, sizeof(float));
    *invcovariance = invcovar;

    totalframes = 0; nfiles = 0;
    while (fscanf(fp,"%s",cepfile) != EOF) {
        if (read_cep(cepfile,&cepfeat,&nfr,basedim) == -1)
            QUIT(("Unable to read feature file %s\n",cepfile))
        mean_normalize(cepfeat, basedim, nfr);
        extcepfeat = extend_vector(cepfeat, basedim, nfr, bracket, &extnfr);
        totalframes += extnfr;
        for (i=0; i<extnfr; i++) {
            for (j=0; j<fulldim; j++) {
                for (k=0; k<=j; k++) {
                    covar[j][k] += extcepfeat[i][j] * extcepfeat[i][k];
                }
            }
        }
        free2d((void**)cepfeat); free2d((void**)extcepfeat);
        nfiles++; if (nfiles%10 == 0) {fprintf(stdout,".\n"); fflush(stdout);}
    }
    rewind(fp);
    *numvecs = totalframes;
    fprintf(stdout,"Computed covariances. %dx%d data.\n",totalframes,fulldim);

    for (j=0; j<fulldim; j++) {
        for (k=0; k<=j; k++) {
            covar[j][k] /= (float)totalframes;
            covar[k][j] = covar[j][k];
        }
    }

    matrixinvert(invcovar, covar, fulldim);

    if ((feat = (float **) alloc2d(totalframes,fulldim,sizeof(float))) == NULL)
        QUIT(("Unable to allocate %d x %d data vectors\n",totalframes,fulldim))
    *vector = feat;

    k = 0; nfiles = 0;
    while (fscanf(fp,"%s",cepfile) != EOF) {
        float **rotatedextfeat;

        if (read_cep(cepfile,&cepfeat,&nfr,basedim) == -1)
            QUIT(("Unable to read feature file %s\n",cepfile))
        mean_normalize(cepfeat, basedim, nfr);
        extcepfeat = extend_vector(cepfeat, basedim, nfr, bracket, &extnfr);
        rotatedextfeat = matmult(invcovar,extcepfeat,fulldim,fulldim,extnfr);
        for (i=0; i<extnfr; i++,k++) {
            for (j=0; j<fulldim; j++) {
                feat[k][j] = rotatedextfeat[i][j];
            }
        }
        free2d((void**)cepfeat); free2d((void**)extcepfeat);
        free2d((void**)rotatedextfeat);
        nfiles++; if (nfiles%10 == 0) {fprintf(stdout,".\n"); fflush(stdout);}
    }
    fprintf(stdout,"Read and whitened all %d frames of input.\n",k);
    return 0;
}


int savefile(char *cbkfile, float **mean, float **invcov, int Nmodes, int Ndim)
{
    int i,j;
    FILE *fp;

    if ((fp = fopen(cbkfile,"w")) == NULL)
        QUIT(("Error opening %s\n",cbkfile))

    fprintf(fp,"%d\n\n",Ndim);
    for (i=0; i<Ndim; i++) {
        for (j=0; j<Ndim; j++) fprintf(fp,"%e ",invcov[i][j]);
        fprintf(fp,"\n");
    }
    for (i=0; i<Nmodes; i++) {
        for (j=0; j<Ndim; j++) fprintf(fp,"%e ",mean[i][j]);
        fprintf(fp,"\n");
    }
    fclose(fp);

    return 0;
}

int main(int argc, char **argv)
{
    float   **vector, **mean, **invcovar;
    int     i, numvecs, Ndim, bracket, Nmodes, maxiter;
    char    *infile, *cbkfile=NULL;

    void    vector_quantize();


    maxiter = 100;
    Nmodes = 0;
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
            case 'w': case 'W': {
                bracket = atoi(argv[++i]);
                break;
            }
            case 'd': case 'D': {
                Ndim = atoi(argv[++i]);
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

    readdata(infile, bracket, &vector, &numvecs, &Ndim, &invcovar);

    mean = (float **) alloc2d(Nmodes,Ndim,sizeof(float));

    vector_quantize(mean,Nmodes,vector,numvecs,Ndim,maxiter);

    savefile(cbkfile, mean, invcovar, Nmodes, Ndim);

    free2d((void**)vector); free2d((void**)invcovar); free2d((void**)mean);

    return 0;
}
