#include "trainclass.h"

// Read class info from the given file. Kinda long..
void read_classes(char *modelfile,
                  int *ngau, int *dim, 
                  float ***mean, float ***vars, float **c)
{
    float **tmean,**tvar,*tmixw;
    int   ndim, tngau, i, j;

    FILE *fp2;

    // Ignore all lines beginning with '#' and find how many classes
    printf("Reading models in %s\n",modelfile);
    if ((fp2=fopen(modelfile,"r")) == NULL)
        QUIT(("Unable to open %s for reading\n",modelfile))
    fscanf(fp2,"%d %d",&tngau,&ndim);
    *ngau = tngau; *dim = ndim;
    tmean = (float **) alloc2d(tngau,ndim,sizeof(float));
    tvar = (float **) alloc2d(tngau,ndim,sizeof(float));
    tmixw = (float *) calloc(tngau,sizeof(float));
    for (i=0; i<tngau; i++){
        fscanf(fp2,"%f",&tmixw[i]);
        for (j=0;j<ndim;j++) fscanf(fp2,"%f",&tmean[i][j]);
        for (j=0;j<ndim;j++) fscanf(fp2,"%f",&tvar[i][j]);
    }
    fclose(fp2);
    *mean = tmean; *vars = tvar; *c = tmixw;
    return;
}
