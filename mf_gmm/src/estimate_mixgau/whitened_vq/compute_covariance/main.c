#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_cblas.h>
#include "header.h"

#define DESCRIPTION "\nThis program computes a VQ codebook from a set of feature vectors.  The data are assumed to be in mfc format. The code will create extended vectors from K consecutive vectors prior to training the codebook.\n"

#define USAGE "USAGE:\n%s -i <input control file> \\\n\t-d <dimensionality of input> \\\n\t-w <bracket -- extended vec = currentframe +- w frames>\\\n\t-o <output codebook file>\n"


float **matmult(float **a, float **b, int rowsa, int colsa, int colsb)
{
    float **c = (float **) alloc2d(colsb, rowsa, sizeof(float));

    cblas_sgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
                rowsa, colsb, colsa,
                1, a[0], rowsa, b[0], colsa, 0, c[0], rowsa);
    return c;
}

float **matmulttxp(float **a, float **b, int rowsa, int colsa, int rowsb)
{
    float **c = (float **) alloc2d(rowsa, rowsa, sizeof(float));

    cblas_sgemm(CblasColMajor,CblasNoTrans,CblasTrans,
                rowsa, rowsa, colsa,
                1, a[0], rowsa, b[0], rowsb, 0, c[0], rowsa);
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

void variance_normalize(float **x, int dim, int nvec)
{
    int i, j;
    float std[1024];

    for (j=0; j<dim; j++) std[j] = 0;

    for (i=0; i<nvec; i++) 
        for (j=0; j<dim; j++) std[j] += x[i][j]*x[i][j];

    for (j=0; j<dim; j++) std[j] = sqrt(std[j]/(float) nvec);

    for (i=0; i<nvec; i++) 
        for (j=0; j<dim; j++) x[i][j] /= std[j];

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
        for (j=0,l=0; l<ndim; l++) {
            y[m][j++] = x[i][l];
        }
        if (k>0) {
            for (k=1;k<=bracket;k++) {
                for (l=0; l<ndim; l++) {
                    y[m][j++] = x[i-k][l];
                }
                for (l=0; l<ndim; l++) {
                    y[m][j++] = x[i+k][l];
                }
            }
        }
    }

    return y;
}


/* Read input features, create extended feature, whiten. */
float  **getcovariance(char *infile, int bracket, int *Ndim)
{
    int totalframes, nfr, extnfr, fulldim, basedim, i, j, k, nfiles;
    char cepfile[1024];
    float **cepfeat, **extcepfeat, **feat, **covar;
    FILE *fp;

    fprintf(stdout,"Reading input.\n");
    basedim = *Ndim;
    *Ndim = fulldim = (2*bracket + 1) * basedim;

    if ((fp = fopen(infile,"r")) == NULL) QUIT(("Cannot open %s\n",infile))

    covar = (float **) alloc2d(fulldim, fulldim, sizeof(float));

    totalframes = 0; nfiles = 0;
    while (fscanf(fp,"%s",cepfile) != EOF) {
        float **lcovar;
        if (read_cep(cepfile,&cepfeat,&nfr,basedim) == -1)
            QUIT(("Unable to read feature file %s\n",cepfile))
        if (nfr < 3*(bracket+1)) { // Bad file; skip
            fprintf(stderr,"%s has only %d frames; skipping.\n",cepfile,nfr);
            free2d((void**)cepfeat);
            continue;
        }
        mean_normalize(cepfeat, basedim, nfr);
        variance_normalize(cepfeat, basedim, nfr);
        extcepfeat = extend_vector(cepfeat, basedim, nfr, bracket, &extnfr);
        totalframes += extnfr;
        lcovar = matmulttxp(extcepfeat,extcepfeat,fulldim,extnfr,fulldim);
        for (j=0; j<fulldim; j++) {
            for (k=0; k<=j; k++) {
                covar[j][k] += lcovar[j][k];
            }
        }
        free2d((void**)cepfeat); free2d((void**)extcepfeat); 
        free2d((void**)lcovar);
        nfiles++; if (nfiles%100 == 0) {fprintf(stdout,"."); fflush(stdout);}
    }
    fprintf(stdout,"\n");
    rewind(fp);
    fprintf(stdout,"Computed covariances. %dx%d data.\n",totalframes,fulldim);

    for (j=0; j<fulldim; j++) {
        for (k=0; k<=j; k++) {
            covar[j][k] /= (float)totalframes;
            covar[k][j] = covar[j][k];
        }
    }

    fprintf(stdout,"Computed covariance from %d frames of input.\n",totalframes);
    return covar;
}

int savefile(char *cbkfile, float **cov, int Ndim)
{
    int i,j;
    FILE *fp;

    if ((fp = fopen(cbkfile,"w")) == NULL)
        QUIT(("Error opening %s\n",cbkfile))

    fprintf(fp,"%d\n\n",Ndim);
    for (i=0; i<Ndim; i++) {
        for (j=0; j<Ndim; j++) fprintf(fp,"%e ",cov[i][j]);
        fprintf(fp,"\n");
    }
    fclose(fp);

    return 0;
}



int main(int argc, char **argv)
{
    float   **covar;
    int     i, numvecs, Ndim, bracket, Nmodes, maxiter;
    char    *infile, *cbkfile=NULL;

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
            default: {
                printf(DESCRIPTION);
                QUIT ((USAGE,argv[0]));
            }
        }
    }
    if (!cbkfile) QUIT(("No output file. Argument -o must be specified\n"));
    covar = getcovariance(infile, bracket, &Ndim);

    savefile(cbkfile, covar, Ndim);

    free2d((void**)covar); 

    return 0;
}
