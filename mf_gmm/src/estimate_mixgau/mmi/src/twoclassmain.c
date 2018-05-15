#include <stdio.h>
#include <stdlib.h>
#include "header.h"

/*
 * DUMMY wrapper.
 * Intended for 2-class data, each of which is stored in an "infile".
 * The code expects exactly 2000 2-dimensional points per class
 */
int main(int argc, char **argv)
{
    float **data, **data1, **data2, ***mean, ***var, **mixw;
    FILE *fp, *fp2;
    int i, ngau;
    int class[2000];
    float classprior[2] = {0.5,0.5}; // priors
    char *infile1, *infile2;
    const int maxiter = 20;
    const int ndim = 2;
    const int nclass = 2;

    if (argc != 4)
        QUIT(("USAGE: %s infile1, infile2, ngau\n",argv[0]))

    data = (float **)alloc2d(2000,ndim,sizeof(float));
    data1 = data;
    data2 = &data[1000];
    for (i=0; i<1000; i++) class[i] = 0;
    for (i=1000; i<2000; i++) class[i] = 1;

    infile1 = argv[1];
    infile2 = argv[2];
    sscanf(argv[3],"%d",&ngau);

    fp = fopen(infile1, "r");
    for (i=0; i<1000; i++) fscanf(fp,"%f",&data1[i][0]);
    for (i=0; i<1000; i++) fscanf(fp,"%f",&data1[i][1]);
    fclose(fp);
    fp2 = fopen(infile2, "r");
    for (i=0; i<1000; i++) fscanf(fp2,"%f",&data2[i][0]);
    for (i=0; i<1000; i++) fscanf(fp2,"%f",&data2[i][1]);
    fclose(fp2);

    mean = (float***)alloc3d(nclass,ngau, ndim, sizeof(float));
    var = (float***)alloc3d(nclass,ngau, ndim, sizeof(float));
    mixw = (float**)alloc2d(nclass,ngau, sizeof(float));

    mgau_em(data1,1000,ndim,ngau,mean[0],var[0],mixw[0],maxiter);
    mgau_em(data2,1000,ndim,ngau,mean[1],var[1],mixw[1],maxiter);
    mgau_mmi(data,2000,ndim,nclass,class,ngau,mean,var,mixw,classprior,maxiter);

    store_distribution("DIST.cls1",ngau,ndim,mixw[0],mean[0],var[0]);    
    store_distribution("DIST.cls2",ngau,ndim,mixw[1],mean[1],var[1]);    

    return 0;
}
