#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "header.h"

#define DESCRIPTION "\nThis program computes a VQ codebook from a set of feature vectors.  The data are assumed to be in mfc format. The code will create extended vectors from K consecutive vectors prior to training the codebook.\n"

#define USAGE "USAGE:\n%s -i <input control file> \\\n\t-d <dimensionality of input> \\\n\t-w <bracket -- extended vec = currentframe +- w frames>\\\n\t-v <Covariance file if any> \\\n\t-o <output codebook file> \\\n\t-c <No of codewords in codebook>\\\n\t-n <Max iterations [Default 100]>\n"


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


float **readcovarfile(char *covarfile, int cdim)
{
    float **covar, **tcovar;
    int readdim, i, j;
    FILE *fp = fopen(covarfile,"r");

    fprintf(stdout,"Reading covariances from %s\n",covarfile);
    if (!fp) QUIT(("Unable to open %s for reading.\n",covarfile))
    fscanf(fp,"%d",&readdim);
    tcovar = (float**) alloc2d(readdim,readdim,sizeof(float));
    covar = (float**) alloc2d(cdim,cdim,sizeof(float));

    for (i=0; i<readdim; i++) 
        for (j=0; j<readdim; j++) fscanf(fp,"%f",&tcovar[i][j]);

    for (i=0; i<cdim; i++) 
        for (j=0; j<cdim; j++) covar[i][j] = tcovar[i][j];

    free2d((void**) tcovar);

    return covar;
}

/* Read input features, create extended feature, whiten. */
int  readdata(char *infile, int bracket, char *covarfile, float eigenthreshold,
              int basedim, float ***vector, int *numvecs, 
              int *Ndimfeat, int *Ndimfull, float ***whiteningmatrix)
{
    int totalframes, nfr, extnfr, fulldim, featdim, i, j, k, nfiles;
    char cepfile[1024];
    float **cepfeat, **extcepfeat, **feat, **covar, **whitener;
    FILE *fp;

    fprintf(stdout,"Reading input.\n");
    *Ndimfull = fulldim = (2*bracket + 1) * basedim;

    if ((fp = fopen(infile,"r")) == NULL) QUIT(("Cannot open %s\n",infile))

    covar = (float **) alloc2d(fulldim, fulldim, sizeof(float));

    totalframes = 0; nfiles = 0;
    while (fscanf(fp,"%s",cepfile) != EOF) {
        if (read_cep(cepfile,&cepfeat,&nfr,basedim) == -1)
            QUIT(("Unable to read feature file %s\n",cepfile))
        if (nfr < 3*(bracket+1)) { // Bad file; skip
            fprintf(stderr,"%s has only %d frames; skip.\n",cepfile,nfr);
            free2d((void**)cepfeat);
            continue;
        }
        totalframes += nfr  - 2*bracket;
        if (covarfile == NULL) { // Compute the covariance
            float **lcovar;
            mean_normalize(cepfeat, basedim, nfr);
            extcepfeat = extend_vector(cepfeat, basedim, nfr, bracket, &extnfr);
            lcovar = matmulttxpA(extcepfeat,fulldim,extnfr,
                                 extcepfeat,fulldim,extnfr,NULL);
            for (j=0; j<fulldim; j++) {
                for (k=0; k<=j; k++) {
                    covar[j][k] += lcovar[j][k];
                }
            }
            free2d((void**)extcepfeat); 
            free2d((void**)lcovar);
            nfiles++; if (nfiles%10 == 0) {fprintf(stdout,"."); fflush(stdout);}
        
        }
        free2d((void**)cepfeat); 
    }

    if (covarfile == NULL){
        fprintf(stdout,"Computed covars. %dx%d data.\n",totalframes,fulldim);
        fprintf(stdout,"\n");
        for (j=0; j<fulldim; j++) {
            for (k=0; k<=j; k++) {
                covar[j][k] /= (float)totalframes;
                covar[k][j] = covar[j][k];
            }
        }
    }
    else { // Read covar file
        covar = readcovarfile(covarfile,fulldim);
    }
    rewind(fp);
    *numvecs = totalframes;

    computewhiteningmatrix(covar, fulldim, eigenthreshold, &featdim, &whitener);
    *Ndimfeat = featdim;
    *whiteningmatrix = whitener;

    if ((feat = (float **) alloc2d(totalframes,featdim,sizeof(float))) == NULL)
        QUIT(("Unable to allocate %d x %d data vectors\n",totalframes,featdim))
    *vector = feat;
    printf("Feature dimensionality %d\n",featdim);

    k = 0; nfiles = 0;
    while (fscanf(fp,"%s",cepfile) != EOF) {
        float **rotatedextfeat;

        if (read_cep(cepfile,&cepfeat,&nfr,basedim) == -1)
            QUIT(("Unable to read feature file %s\n",cepfile))
        if (nfr < 3*(bracket+1)) { // Bad file; skip
            //fprintf(stderr,"%s has only %d frames; skipping.\n",cepfile,nfr);
            free2d((void**)cepfeat);
            continue;
        }
        mean_normalize(cepfeat, basedim, nfr);
        extcepfeat = extend_vector(cepfeat, basedim, nfr, bracket, &extnfr);
        rotatedextfeat = matmulttxpB(extcepfeat,extnfr,fulldim,
                                     whitener,featdim,fulldim,NULL);
        for (i=0; i<extnfr; i++,k++) {
            for (j=0; j<featdim; j++) {
                feat[k][j] = rotatedextfeat[i][j];
            }
        }
        free2d((void**)cepfeat); free2d((void**)extcepfeat);
        free2d((void**)rotatedextfeat);
        nfiles++; if (nfiles%10 == 0) {fprintf(stdout,"."); fflush(stdout);}
    }
    fprintf(stdout,"\n");
    fprintf(stdout,"Read and whitened all %d frames of input.\n",k);
    return 0;
}


int savefile(char *cbkfile, float **mean, float **whiteningmat, 
             int Ndimfeat, int Ndimfull, int cepdim, int Nmodes)
{
    int i,j;
    FILE *fp;

    if ((fp = fopen(cbkfile,"w")) == NULL)
        QUIT(("Error opening %s\n",cbkfile))

    fprintf(fp,"%d %d %d %d\n\n",Nmodes, Ndimfeat, Ndimfull, cepdim);
    for (i=0; i<Ndimfeat; i++) {
        for (j=0; j<Ndimfull; j++) fprintf(fp,"%e ",whiteningmat[i][j]);
        fprintf(fp,"\n");
    }
    for (i=0; i<Nmodes; i++) {
        for (j=0; j<Ndimfeat; j++) fprintf(fp,"%e ",mean[i][j]);
        fprintf(fp,"\n");
    }
    fclose(fp);

    return 0;
}

int main(int argc, char **argv)
{
    float   **vector, **mean, **whiteningmat;
    int     i, numvecs, cepdim, Ndimfeat, Ndimfull;
    int     bracket, Nmodes, maxiter;
    char    *infile, *cbkfile=NULL, *covarfile=NULL;

    void    vector_quantize();

    float   eigenthreshold = 4e-04;


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
                cepdim = atoi(argv[++i]);
                break;
            }
            case 'v': case 'V': {
                covarfile = argv[++i];
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

    readdata(infile, bracket, covarfile, eigenthreshold, cepdim,
             &vector, &numvecs, &Ndimfeat, &Ndimfull, &whiteningmat);

    mean = (float **) alloc2d(Nmodes,Ndimfeat,sizeof(float));
    //* Save the initial models
    savefile(cbkfile, mean, whiteningmat, Ndimfeat, Ndimfull, cepdim, Nmodes);

    vector_quantize(mean,Nmodes,vector,numvecs,Ndimfeat,maxiter);

    savefile(cbkfile, mean, whiteningmat, Ndimfeat, Ndimfull, cepdim, Nmodes);

    free2d((void**)vector); free2d((void**)whiteningmat); free2d((void**)mean);

    return 0;
}
