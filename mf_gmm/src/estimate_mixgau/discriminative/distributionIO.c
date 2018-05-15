#include "discriminative.h"

void  read_distribution(char*,int*,int*,float**,float***,float***);
void  store_distribution(char*,int,int,float*,float**,float**);

void readalldistributions(char *ctllist, int *onclass,
                          int **onmodes, int *ondim,
                          float ***oc, float ****omean, float ****ovar)
{
    int i, j, k, nclass, *nmodes, ndim;
    float ***tmean, ***mean, ***tvar, ***var, **c, **tc;
    char fname[512];
    FILE *fp;

    // First find number of classes
    if ((fp = fopen(ctllist,"r")) == NULL)
        QUIT(("Unable to open %s\n",ctllist))
    for (nclass=0; fscanf(fp,"%s",fname) != EOF; nclass++); rewind(fp);
    *onclass = nclass;

    // read all into a temporary list
    tmean = (float ***)calloc(nclass,sizeof(float**));
    tvar = (float ***)calloc(nclass,sizeof(float**));
    tc = (float **)calloc(nclass,sizeof(float*));
    *onmodes = nmodes = (int *)calloc(nclass, sizeof(int));

    fscanf(fp,"%s",fname);
    read_distribution(fname,nmodes,ondim,tc,tmean,tvar);
    for (i=1; i<nclass; i++){
        fscanf(fp,"%s",fname);
        read_distribution(fname,nmodes+i,&ndim,tc+i,tmean+i,tvar+i);
        if (ndim != *ondim)
            QUIT(("All models do not have same dimension\n"))
    }
    fclose(fp);

    // Redo in alloc3d format. Assume all models have same nmodes
    *omean = mean = (float***) alloc3d(nclass,nmodes[0],ndim,sizeof(float));
    *ovar = var = (float***) alloc3d(nclass,nmodes[0],ndim,sizeof(float));
    *oc = c = (float**) alloc2d(nclass,nmodes[0],sizeof(float));

    for (i=0; i<nclass; i++) {
        for (j=0; j<nmodes[i]; j++) {
            c[i][j] = tc[i][j];
            for (k=0; k<ndim; k++) {
                mean[i][j][k] = tmean[i][j][k];
                var[i][j][k] = tvar[i][j][k];
            }
        }
        free2d((void**)tmean[i]); free2d((void**)tvar[i]); free(tc[i]);
    }
    free(tmean); free(tvar); free(tc);

    return;
}

void storealldistributions(char *ctllist, int nclass,
                           int *nmodes, int ndim,
                           float **c, float ***mean, float ***var)
{
    int i, tnclass;
    char fname[512];
    FILE *fp;

    // First find number of classes
    if ((fp = fopen(ctllist,"r")) == NULL)
        QUIT(("Unable to open %s\n",ctllist))
    for (tnclass=0; fscanf(fp,"%s",fname) != EOF; tnclass++); rewind(fp);
    if (tnclass != nclass) {
        printf("Only %d distributions given in %s\n",tnclass,ctllist);
	printf("Only %d of %d distributions will be written\n",tnclass,nclass);
	fflush(stdout);
    }

    for (i=0; i<tnclass; i++) {
        fscanf(fp,"%s",fname);
        store_distribution(fname,nmodes[i],ndim,c[i],mean[i],var[i]);
    }
    fclose(fp);

    return;
}

/* Routine to store the distributions computed */

void  store_distribution(char *file, int Nmodes, int Ndim, 
                         float *c, float **mean, float **var)
{
    int     i,j;
    FILE    *fopen(),*outputfile;

    if ((outputfile = fopen(file,"w")) == NULL) {
        printf("Unable to open %s for writing!\n",file);
        printf("Will store models in /tmp/DISTRIB instead\n");
        if ((outputfile = fopen("/tmp/DISTRIB","w")) == NULL) 
            QUIT(("Unable to open /tmp/DISTRIB either! No output!\n"))
    }
    
    fprintf(outputfile,"%d %d\n",Nmodes,Ndim);
    for (i=0;i<Nmodes;++i)
    {
        fprintf(outputfile,"%e\n",c[i]);
        for (j=0;j<Ndim;++j)
            fprintf(outputfile,"%e ",mean[i][j]);
        fprintf(outputfile,"\n");
        for (j=0;j<Ndim;++j)
            fprintf(outputfile,"%e ",var[i][j]);
        fprintf(outputfile,"\n");
    }
    fclose(outputfile);
    return;
}

/*========================================================================*/
/* Routine to read initial distributions, if any */

void  read_distribution(char *file, int *nmodes, int *ndim, 
                         float **c, float ***mean, float ***var)
{
    int     i,j;
    float *lc, **lmean, **lvar;
    FILE    *fopen(),*inputfile;

    if ((inputfile = fopen(file,"r")) == NULL) 
        QUIT(("Unable to open %s for reading\n",file))
    
    fscanf(inputfile,"%d %d",nmodes,ndim);
    *c = lc = (float *)calloc(*nmodes,sizeof(float));
    *mean = lmean = (float**) alloc2d(*nmodes,*ndim,sizeof(float));
    *var = lvar = (float**) alloc2d(*nmodes,*ndim,sizeof(float));

    for (i=0;i<*nmodes;++i)
    {
        fscanf(inputfile,"%f",&lc[i]);
        for (j=0;j<*ndim;++j)
            fscanf(inputfile,"%f",&lmean[i][j]);
        for (j=0;j<*ndim;++j)
            fscanf(inputfile,"%f",&lvar[i][j]);
    }
    fclose(inputfile);
    return;
}

void readpriors(char *priorsfile, int nclass, float **oclasspriors)
{
    int i;
    float *classpriors;
    FILE *fp;

    *oclasspriors = classpriors = (float *) calloc(nclass, sizeof(float));
    if (priorsfile != NULL) {
        if ((fp = fopen(priorsfile,"r")) == NULL)
            QUIT(("Unable to open %s for reading\n",priorsfile))

        for (i=0; fscanf(fp,"%f",&classpriors[i]) != EOF; i++);
        fclose(fp);

        if (i != nclass)
            QUIT(("%d priors required, read %d priors\n",nclass,i))
    }
    else
        for (i=0; i<nclass; i++) classpriors[i] = 1.0/(float)nclass;

    return;
}
