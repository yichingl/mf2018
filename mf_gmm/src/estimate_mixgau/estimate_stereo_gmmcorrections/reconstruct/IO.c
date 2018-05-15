#include <stdio.h>
#include <stdlib.h>
#include "header.h"

int    read_input_distribution(char *filen,
                               float ***mean,
                               float ***variance,
                               float **pi,
                               int *Nmodes,
                               int *Ndim,
                               int *L, int *R,
                               float ****meany)
{
    int    i,j,k,imodes,idim, lc, rc;
    float  **varbuff,**meanbuff,*cbuff, ***lmeany;

    FILE   *fopen(),*dist;

    if ((dist = fopen(filen,"r")) == NULL) {
        printf("Unable to open distribution file\n");
        fflush(stdout);
        return(-1);
    }

    fscanf(dist,"%d %d",&imodes,&idim);

    *Nmodes = imodes;
    *Ndim = idim;

    if ((varbuff = (float **) alloc2d(imodes,idim,sizeof(float))) == NULL) {
        printf("Unable to allocate space for variances\n");
        fflush(stdout);
        return(-1);
    }

    if ((meanbuff = (float **) alloc2d(imodes,idim,sizeof(float))) == NULL) {
        printf("Unable to allocate space for means\n");
        fflush(stdout);
        return(-1);
    }

    if ((cbuff = (float *) malloc(imodes*sizeof(float))) == NULL) {
        printf("Unable to allocate space for priors\n");
        fflush(stdout);
        return(-1);
    }

    for (k=0;k<imodes;++k) {
        fscanf(dist,"%f",&cbuff[k]);
        for (j=0;j<idim;++j)
            fscanf(dist,"%f ",&meanbuff[k][j]);
        for (j=0;j<idim;++j)
            fscanf(dist,"%f ",&varbuff[k][j]);
    }

    fscanf(dist,"%d %d",&lc,&rc);
    *L = lc;  *R = rc;
    lmeany = (float***) alloc3d(lc+rc+1,imodes,idim,sizeof(float));
    *meany = lmeany;
    for (k=-lc; k<=rc; k++) {
        for (i=0; i<imodes; ++i) {
            for (j=0; j<idim; ++j)
                fscanf(dist,"%e",&lmeany[k+lc][i][j]);
        }
    }
    fclose(dist);

    printf("Successfully read distribution\n");fflush(stdout);
    *variance = varbuff; *mean = meanbuff; *pi = cbuff;

    return(0);
}


int write_file(int ncopies, int nfr, int ndim, float ***v, char *fname)
{
    int i, j, k;
    FILE *fp;

    if ((fp = fopen(fname,"w")) == NULL) {
        QUIT(("Unable to open %s for writing\n",fname))
    }
    for (i=0; i<ncopies; i++){
        for (j=0; j<nfr; j++) {
            for (k=0; k<ndim; k++) {
                fprintf(fp,"%e ",v[i][j][k]);
            }
            fprintf(fp,"\n");
        } 
    }

    return 0;
}
