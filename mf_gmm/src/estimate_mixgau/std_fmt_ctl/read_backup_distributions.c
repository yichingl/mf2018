#include <stdio.h>
#include <stdlib.h>

int    read_backup_distribution(filen,mean,variance,pi,Nmodes,Ndim)
float  ***mean,***variance,**pi;
int    *Nmodes,Ndim;
char   *filen;
{
    int    j,k,imodes,idim;
    float  **varbuff,**meanbuff,*cbuff;
    char   **alloc2d();

    FILE   *fopen(),*dist;

    if ((dist = fopen(filen,"r")) == NULL)
    {
        printf("Unable to open distribution file\n");
        return(0);
    }

    fscanf(dist,"%d %d",&imodes,&idim);

    *Nmodes = imodes;

    if (idim != Ndim)
    {
	printf("ERROR in read_distribution:\n");
        printf("Distribution dimensionality != Dimensionality of data given\n");
	return(0);
    }

    if ((varbuff = (float **) alloc2d(imodes,idim,sizeof(float))) == NULL)
    {
        printf("Unable to allocate space for variances\n");
        return(0);
    }

    if ((meanbuff = (float **) alloc2d(imodes,idim,sizeof(float))) == NULL)
    {
        printf("Unable to allocate space for means\n");
        return(0);
    }

    if ((cbuff = (float *) malloc(imodes*sizeof(float))) == NULL)
    {
        printf("Unable to allocate space for priors\n");
        return(0);
    }

    for (k=0;k<imodes;++k)
    {
            fscanf(dist,"%f",&cbuff[k]);
            for (j=0;j<idim;++j)
                fscanf(dist,"%f ",&meanbuff[k][j]);
            for (j=0;j<idim;++j)
                 fscanf(dist,"%f ",&varbuff[k][j]);
    }
    fclose(dist);

    *variance = varbuff; *mean = meanbuff; *pi = cbuff;

    return(1);
}
