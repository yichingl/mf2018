#include <stdio.h>
#include <stdlib.h>
#include "header.h"

int    read_input_distribution(char *filen,
                               float ***mean,
                               float ***variance,
                               float **pi,
                               int *Nmodes,
                               int *Ndim)
{
    int    j,k,imodes,idim;
    float  **varbuff,**meanbuff,*cbuff;

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
    fclose(dist);

    printf("Successfully read distribution\n");fflush(stdout);
    *variance = varbuff; *mean = meanbuff; *pi = cbuff;

    return(0);
}



int  store_buff(char *file, int Nmodes, int Ndim,
                float *c, float **mean, float **var,
                int leftctx, int rightctx, float ***meany, float *count)
{
    int     i,j,k;

    FILE    *fopen(),*outputfile;

    if ((outputfile = fopen(file,"w")) == NULL) return(-1);
    
    fprintf(outputfile,"%d %d\n",Nmodes,Ndim);
    for (i=0;i<Nmodes;++i) {
        fprintf(outputfile,"%e\n",c[i]);
        for (j=0;j<Ndim;++j)
            fprintf(outputfile,"%e ",mean[i][j]);
        fprintf(outputfile,"\n");
        for (j=0;j<Ndim;++j)
            fprintf(outputfile,"%e ",var[i][j]);
        fprintf(outputfile,"\n");
    }

    for (i=0;i<Nmodes;++i) fprintf(outputfile,"%e ",count[i]);
    fprintf(outputfile,"\n");

    fprintf(outputfile,"\n%d %d\n",leftctx,rightctx);
    for (k=-leftctx;k<=rightctx;k++) {
        for (i=0;i<Nmodes;++i) {
            for (j=0;j<Ndim;++j)
                fprintf(outputfile,"%e ",meany[k+leftctx][i][j]);
            fprintf(outputfile,"\n");
        }
        fprintf(outputfile,"\n");
    }
    fclose(outputfile);
    printf("Wrote distribution to %s\n",file); fflush(stdout);
    return(0);
}
