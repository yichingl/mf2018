#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "header.h"

#define DESCRIPTION "\nThis program computes a mixture gaussian distribution for a set of vectors. \nThe data are assumed to be in an ascii file, of format\nf(0,0) f(0,1) ... f(0,D-1)\nf(1,0) f(1,1) .. f(1,D-1)\n...\nf(N-1,0).......f(N-1,D-1)\n where D is the dimensionality of the vectors and N is the number of vectors in the data set\nIf the input is not in the right format the program will die horribly.\n"

#define USAGE "USAGE:\n%s -i <input data file> \\\n\t-c <output codebook file> \\\n\t-o <output index file for inputs> \\\n\t-m <Total no of modes in distribution>\\\n\t-n <Max iterations [Default 100]>\n"

int readfile(char *infile, float ***buff, int *onvec, int *ondim)
{
    int i, j, nvec, ndim;
    float **lbuff;
    char line[10000];
    FILE *ifp;

    ifp = fopen(infile,"r");
    if (ifp == NULL) 
        QUIT(("Unable to open input file %s\n",infile));

    // First read a single line and calculate dimension
    fgets(line,10000,ifp);
    ndim = 0;
    if(strtok(line," \t")) {
        ndim++;  // Count the first word
        while(strtok(0," \t")) ndim++; // Count the rest
    }
    *ondim = ndim;
    rewind(ifp);

    // Count the no. of vectors;
    nvec = 0;
    while (fgets(line,10000,ifp) != NULL) nvec++;
    *onvec = nvec;
    rewind(ifp);

    // Allocate and read data
    *buff = lbuff = (float **) alloc2d(nvec, ndim, sizeof(float));
    for (i=0; i<nvec; i++)
        for (j=0; j<ndim; j++)
            fscanf(ifp,"%f",&lbuff[i][j]);

    fclose(ifp);

    return 0;
}

int savefiles(char *cbkfile, char *idxfile, float **mean, int Nmodes, int Ndim,
              int *bin, int nvecs)
{
    int i,j;
    FILE *fp;

    if (idxfile) {
        if ((fp = fopen(idxfile,"w")) == NULL)
            QUIT(("Error opening %s\n",idxfile))
        for (i=0; i<nvecs; i++) fprintf(fp,"%d\n",bin[i]);
        fclose(fp);
    }

    if (cbkfile) {
        if ((fp = fopen(cbkfile,"w")) == NULL)
            QUIT(("Error opening %s\n",cbkfile))
        for (i=0; i<Nmodes; i++) {
            for (j=0; j<Ndim; j++) fprintf(fp,"%e ",mean[i][j]);
            fprintf(fp,"\n");
        }
        fclose(fp);
    }
    return 0;
}

int main(int argc, char **argv)
{
    float   **vector, **mean;
    int     i, numvecs, Ndim, Nmodes, *bin, maxiter;
    char    *infile, *cbkfile=NULL, *idxfile=NULL;

    void    vector_quantize();


    maxiter = 100;
    /* Print out command, just to be nice */
    for (i=0; i<argc; i++) printf("%s ",argv[i]); printf("\n"); fflush(stdout);
    if (argc == 1)
    {
        printf(DESCRIPTION);
        QUIT ((USAGE,argv[0]));
    }
    for (i=1;i<argc;++i)
    {
        char sw;
        if (argv[i][0]!='-')
            QUIT ((USAGE,argv[0]));
        sw = argv[i][1];
        switch(sw)
        {
            case 'i':
            {
                infile = argv[++i];
                break;
            }
            case 'c': case 'C':
            {
                cbkfile = argv[++i];
                break;
            }
            case 'o': case 'O':
            {
                idxfile = argv[++i];
                break;
            }
            case 'm': case 'M':
            {
                Nmodes = atoi(argv[++i]);
                break;
            }
            case 'n': case 'N':
            {
                maxiter = atoi(argv[++i]);
                break;
            }
            default:
            {
                printf(DESCRIPTION);
                QUIT ((USAGE,argv[0]));
            }

        }
    }
    if (!cbkfile && !idxfile)
        QUIT(("No output files given. One of -c or -o must be specified\n"));

    readfile(infile, &vector, &numvecs, &Ndim);

    mean = (float **) alloc2d(Nmodes,Ndim,sizeof(float));

    if ((bin = (int *)calloc(numvecs, sizeof(int))) == NULL)
        QUIT(("Unable to allocate %d ints for bin array\n",numvecs));

    vector_quantize(mean,Nmodes,vector,numvecs,Ndim,bin,maxiter);

    savefiles(cbkfile, idxfile, mean, Nmodes, Ndim, bin, numvecs);

    free2d((void**)vector);
    free2d((void**)mean);
    free((void*)bin);

    return 0;
}
