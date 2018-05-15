#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include "header.h"

#define MAXLINESIZE 10000

int readcodebook(char *codefile, float ***Codes, int *Ncode, int *Ndim)
{
    int ncode, ndim, i, j;
    float **codes;
    FILE *fp = fopen(codefile,"r");
    if (fp == NULL) QUIT(("Unable to open %s for reading\n",codefile));

    fscanf(fp,"%d",&ncode); *Ncode = ncode;
    fscanf(fp,"%d",&ndim); *Ndim = ndim;

    *Codes = codes = (float **) alloc2d(ncode, ndim, sizeof(float));

    for (i=0;i<ndim;i++){
        for (j=0;j<ncode;j++) fscanf(fp,"%f",&codes[j][i]);
        fprintf(fp,"\n");
    }
    fclose (fp);

    return 1;
}


main(int argc, char **argv)
{
    float **vec, **extendedvec, **codes;
    int  *bow, Ndim, Nvec, Ncodes, Nfeat, i, fileno, width, cepdim;
    char *ctlfile,*outfile,*cbkfile,line[MAXLINESIZE];
    FILE *fctl,*fout;

    if (argc != 6) {
        printf("Usage: %s <ctl file> <width> <cepdim> <codebook file> <output file>>\n",argv[0]);
       exit(0);
    }
    ctlfile = argv[1]; cbkfile = argv[4]; outfile = argv[5];
    sscanf(argv[2],"%d",&width); sscanf(argv[3],"%d",&cepdim);

    readcodebook(cbkfile, &codes, &Ncodes, &Ndim);
    assert(Ndim == cepdim*width);
    bow = (int *) calloc(Ncodes, sizeof(int));
 
    if ((fctl=fopen(ctlfile,"r"))==NULL) QUIT(("Unable to read %s\n",ctlfile));
    if ((fout=fopen(outfile,"w"))==NULL) QUIT(("Unable to write %s\n",outfile));

    fileno = 0;
    while (fscanf(fctl,"%s",line) != EOF) {
        read_cep(line, &vec, &Nvec, cepdim);
        extend_vec(vec,cepdim,Nvec,width,&extendedvec,&Nfeat);
        fileno++; printf("%d\n", fileno); fflush(stdout);
        BagofWords(bow, extendedvec, Nfeat, codes, Ncodes, Ndim);
        for (i=0; i<Ncodes; i++) fprintf(fout,"%d ",bow[i]); fprintf(fout,"\n");
        free2d((void**) vec); free2d((void**) extendedvec);
    }
    fclose(fout); fclose(fctl);

    return 0;
}
