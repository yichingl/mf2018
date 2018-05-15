#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "header.h"

#define DESCRIPTION "\nThis program computes a VQ codebook from a set of feature vectors.  The data are assumed to be in mfc format. The code will create extended vectors from K consecutive vectors prior to training the codebook.\n"

#define USAGE "USAGE:\n%s -i <input control file> \\\n\t-c <codebook file> \\\n\t-o <output file basename>\n"

int readmodel(char *cbkfile, float ***mean, float ***whiteningmat, 
              int *Nmodes, int *Ndimfeat, int *Ndimfull, int *cepdim)
{
    int i, j, nmodes, ndimfull, ndimfeat, basedim;
    float **codes, **whitener;
    FILE *fp;

    if ((fp = fopen(cbkfile,"r")) == NULL)
        QUIT(("Error opening %s\n",cbkfile))

    fscanf(fp,"%d",&nmodes); *Nmodes = nmodes;
    fscanf(fp,"%d",&ndimfeat); *Ndimfeat = ndimfeat;
    fscanf(fp,"%d",&ndimfull); *Ndimfull = ndimfull;
    fscanf(fp,"%d",&basedim); *cepdim = basedim;

    *mean = codes = (float **) alloc2d(nmodes, ndimfeat, sizeof(float));
    *whiteningmat = whitener = (float **) alloc2d(ndimfeat, ndimfull, sizeof(float));
    
    for (i=0; i<ndimfeat; i++) {
        for (j=0; j<ndimfull; j++) fscanf(fp,"%f",&whitener[i][j]);
    }
    for (i=0; i<nmodes; i++) {
        for (j=0; j<ndimfeat; j++) fscanf(fp,"%f",&codes[i][j]);
    }
    fclose(fp);

    return 0;
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
int  readdata(char *infile, int bracket, int basedim, 
              float **whitener, int featdim, int fulldim,
              float ***vector, int *numvecs)
{
    int nfr, extnfr;
    float **cepfeat, **extcepfeat, **rotatedextfeat;

    fprintf(stdout,"Reading %s.\n", infile);
    assert(fulldim == (2*bracket + 1)*basedim);

    if (read_cep(infile,&cepfeat,&nfr,basedim) == -1) {
        fprintf(stderr,"Unable to read feature file %s\n",infile);
        return -1;
    }
    if (nfr < 3*(bracket+1)) { // Bad file; skip
        fprintf(stderr,"%s has only %d frames; skipping.\n",infile,nfr);
        free2d((void**)cepfeat);
        return -1;
    }
    mean_normalize(cepfeat, basedim, nfr);
    extcepfeat = extend_vector(cepfeat, basedim, nfr, bracket, &extnfr);
    rotatedextfeat = matmulttxpB(extcepfeat,extnfr,fulldim,
                                 whitener,featdim,fulldim,NULL);
    *vector = rotatedextfeat;
    *numvecs = extnfr;
    free2d((void**)cepfeat); free2d((void**)extcepfeat);

    return 0;
}


int main(int argc, char **argv)
{
    float   **vector, **mean, **whiteningmat, *projs;
    int     i, numvecs, Ndimfeat, Ndimfull, bracket, Nmodes, cepdim, *counts;
    char    *cbkfile=NULL,*ctlfile=NULL,*ofile=NULL;
    char    featfile[1024], outfile1[1024], outfile2[1024]; 

    void    vector_quantize();
    FILE    *fp, *ofp1, *ofp2;

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
                ctlfile = argv[++i];
                break;
            }
            case 'c': {
                cbkfile = argv[++i];
                break;
            }
            case 'o': case 'O': {
                ofile = argv[++i];
                break;
            }
            default: {
                printf(DESCRIPTION);
                QUIT ((USAGE,argv[0]));
            }
        }
    }

    readmodel(cbkfile,&mean,&whiteningmat,&Nmodes,&Ndimfeat,&Ndimfull,&cepdim);
    bracket = ((Ndimfull/cepdim) - 1)/2;
    counts = (int *) calloc(Nmodes, sizeof(int));
    projs = (float *) calloc(Nmodes, sizeof(float));

    sprintf(outfile1,"%s.counts",ofile);
    sprintf(outfile2,"%s.projs",ofile);
    if ((fp = fopen(ctlfile,"r"))==NULL)
        QUIT(("Unable to open %s for reading.\n",ctlfile))
    if ((ofp1 = fopen(outfile1,"w"))==NULL)
        QUIT(("Unable to open %s for writing.\n",outfile1))
    if ((ofp2 = fopen(outfile2,"w"))==NULL)
        QUIT(("Unable to open %s for writing.\n",outfile2))

    while(fscanf(fp,"%s",featfile) != EOF) {
        if (readdata(featfile,bracket,cepdim,whiteningmat,Ndimfeat,Ndimfull,
                     &vector,&numvecs) == -1)
            QUIT(("Unable to read %s\n",featfile))
        bagofwords(mean,Nmodes,vector,numvecs,Ndimfeat,counts,projs);
        for (i=0; i<Nmodes; i++) fprintf(ofp1,"%d ",counts[i]);
        for (i=0; i<Nmodes; i++) fprintf(ofp2,"%f ",projs[i]);
        fprintf(ofp1,"\n"); fprintf(ofp2,"\n");
        free2d((void**)vector);
    }

    fclose(fp); fclose(ofp1); fclose(ofp2);

    free2d((void**)whiteningmat); free2d((void**)mean);

    return 0;
}
