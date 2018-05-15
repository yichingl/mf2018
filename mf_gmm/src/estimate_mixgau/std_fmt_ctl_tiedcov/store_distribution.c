#include <stdio.h>
#include "header.h"


int  store_distribution(file,Nmodes,Ndim,c,mean,var)
int     Nmodes,Ndim;
float   **mean,**var,c[];
char    *file;
{
    int     i,j;

    FILE    *fopen(),*outputfile;

    if ((outputfile = fopen(file,"w")) == NULL) return(-1);
    
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
    printf("Wrote distribution to %s\n",file); fflush(stdout);
    return(0);
}
