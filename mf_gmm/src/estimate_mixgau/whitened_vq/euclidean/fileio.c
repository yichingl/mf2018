#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "header.h"

float **readcovarfile(char *covarfile, int cdim)
{
    float **covar, **tcovar;
    int readdim, i, j;
    FILE *fp = fopen(covarfile,"r");

    fprintf(stdout,"Reading covariances from %s\n",covarfile);
    if (!fp) QUIT(("Unable to open %s for reading.\n",covarfile))
    fscanf(fp,"%d",&readdim);
    if (readdim < cdim)
        QUIT(("Covariance matrix has fewer dims %d than required.\n",readdim))

    tcovar = (float**) alloc2d(readdim,readdim,sizeof(float));
    covar = (float**) alloc2d(cdim,cdim,sizeof(float));

    for (i=0; i<readdim; i++)
        for (j=0; j<readdim; j++) fscanf(fp,"%f",&tcovar[i][j]);

    for (i=0; i<cdim; i++)
        for (j=0; j<cdim; j++) covar[i][j] = tcovar[i][j];

    free2d((void**) tcovar);

    return covar;
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

/* This assumes its already zero mean */
void variance_normalize(float **x, int dim, int nvec)
{
    int i, j;
    float stdev[1024];

    for (j=0; j<dim; j++) stdev[j] = 0;

    for (i=0; i<nvec; i++) 
        for (j=0; j<dim; j++) stdev[j] += x[i][j]*x[i][j];

    for (j=0; j<dim; j++) stdev[j] = sqrt(stdev[j] / (float) nvec);

    for (i=0; i<nvec; i++) 
        for (j=0; j<dim; j++) x[i][j] /= stdev[j];

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

#if 1  // Conventional arrangment may be slow
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
#else
    for (i=bracket, m=0; i<nvecs-bracket; i++,m++) {
        if (k <= 0) {
            for (j=0,l=0; l<ndim; l++) {
                y[m][j++] = x[i][l];
            }
        }
        else {
            // Weird ordering for speed
            for (j=0,l=0; l<ndim; l++) {
                y[m][j++] = x[i][l];
                for (k=1;k<=bracket;k++) {
                    y[m][j++] = x[i-k][l];
                    y[m][j++] = x[i+k][l];
                }
            }
        }
    }
#endif

    return y;
}

/* Read input features, create extended feature, whiten. */
int  readdata(char *infile, int bracket, int basedim, int maxvecs,
              float ***vector, int *numvecs, char *covarfile, int featdim,
              float ***whiteningmat)
{
    int totalframes, nfr, extnfr, fulldim, i, j, k, skip, skipresidue, nfiles;
    char cepfile[1024];
    float **cepfeat, **extcepfeat, **feat, **covar, **whitener;
    FILE *fp;

    if (maxvecs <=0) maxvecs = 10000000; // Maximum 10 million vectors

    fprintf(stdout,"Reading input.\n");
    fulldim = (2*bracket + 1) * basedim;

    if ((fp = fopen(infile,"r")) == NULL) QUIT(("Cannot open %s\n",infile))

    totalframes = 0;
    while (fscanf(fp,"%s",cepfile) != EOF) {
        if ((nfr =  get_featfile_size(cepfile, basedim)) == -1)
            QUIT(("Unable to read feature file %s\n",cepfile))
        if (nfr < 3*(bracket+1)) { // Bad file; skip
            fprintf(stderr,"%s has only %d frames; skip.\n",cepfile,nfr);
            continue;
        }
        totalframes += nfr - 2*bracket;
    }
    rewind(fp);

    if (totalframes < maxvecs) maxvecs = totalframes;

    *numvecs = maxvecs;
    if ((feat = (float **) alloc2d(maxvecs,fulldim,sizeof(float))) == NULL)
        QUIT(("Unable to allocate %d x %d data vectors\n",totalframes,fulldim))
    *vector = feat;
    printf("Feature dimensionality %d\n",fulldim);

    skip = totalframes/maxvecs;
    printf("Picking %d from %d vectors; skip = %d\n",maxvecs,totalframes,skip);

    k = 0; skipresidue = 0; nfiles = 0;
    while (fscanf(fp,"%s",cepfile) != EOF) {
        nfr =  get_featfile_size(cepfile, basedim);
        if (nfr < 3*(bracket+1)) continue; // File too short
        extnfr = nfr-2*bracket;
        if (skipresidue > extnfr) {
            skipresidue -= extnfr; continue; // Not going to read from here
        }
        read_cep(cepfile,&cepfeat,&nfr,basedim);
        mean_normalize(cepfeat, basedim, nfr);
        variance_normalize(cepfeat, basedim, nfr);
        extcepfeat = extend_vector(cepfeat, basedim, nfr, bracket, &extnfr);

        for (i=skipresidue; i<extnfr && k<maxvecs; i+=skip,k++) {
            for (j=0; j<fulldim; j++) {
                feat[k][j] = extcepfeat[i][j];
            }
            skipresidue = skip - (extnfr - i); // vecs to skip in next file
        }
        free2d((void**)cepfeat); free2d((void**)extcepfeat);
        nfiles++; if (nfiles%100 == 0) {printf("\r%d",k); fflush(stdout);}
    }
    fprintf(stdout,"\rRead %d frames of input\n",k); fflush(stdout);

    whitener = NULL;
    if (covarfile) {
        covar = readcovarfile(covarfile,fulldim);
        whitener = whiteningmatrix(covar, fulldim, featdim);
        *whiteningmat = whitener;
        *vector = matmulttxpB(feat,maxvecs,fulldim,
                              whitener,featdim,fulldim,NULL);
        free2d((void**)feat);
        fprintf(stdout,"\nWhitened %d frames of input\n",k); fflush(stdout);
    }

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
    if (whiteningmat != NULL) {
        for (i=0; i<Ndimfeat; i++) {
            for (j=0; j<Ndimfull; j++) fprintf(fp,"%e ",whiteningmat[i][j]);
            fprintf(fp,"\n");
        }
    }
    else { // Identity
        for (i=0; i<Ndimfeat; i++) {
            for (j=0; j<i; j++) fprintf(fp,"0 ");
            fprintf(fp,"1.0 ");
            for (j=i+1; j<Ndimfeat; j++) fprintf(fp,"0e ");
            fprintf(fp,"\n");
        }
    }
    for (i=0; i<Nmodes; i++) {
        for (j=0; j<Ndimfeat; j++) fprintf(fp,"%e ",mean[i][j]);
        fprintf(fp,"\n");
    }
    fclose(fp);

    return 0;
}

