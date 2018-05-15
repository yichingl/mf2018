#include "discriminative.h"

// Read data for all classes
void  getfeatures(char*, float***, int*, int*, int);

void readallclassdata(char *ctllist, int dim, 
                      float *****odata, // nclass, nseqs, len, dim
                      int *onclass, 
                      int **onseqs,    // nclas
                      int ***olenseq)  // nclass, nseqs
{
    float ****data;
    int i, j, nclass, *nseqs, **len, fdim;
    char filen[512], filen2[512];
    FILE *fctl, *fp;

    if ((fctl = fopen(ctllist,"r")) == NULL)
        QUIT(("Unable to open %s for reading\n",ctllist))
    for (nclass=0; fscanf(fctl,"%s",filen) != EOF; nclass++); rewind(fctl);
    *onclass = nclass;
    *onseqs = nseqs = (int *)calloc(nclass,sizeof(int));
    *odata  = data = (float ****)calloc(nclass,sizeof(float***));
    *olenseq = len = (int **) calloc(nclass,sizeof(int *));

    for (i=0; i<nclass; i++) {
        fscanf(fctl,"%s",filen);
        if ((fp = fopen(filen,"r")) == NULL)
            QUIT(("Unable to open %s for reading\n",filen))
        for (j=0; fscanf(fp,"%s",filen2) != EOF; j++); rewind(fp);
        nseqs[i] = j;
        data[i] = (float ***)calloc(nseqs[i],sizeof(float **));
        len[i] = (int *)calloc(nseqs[i],sizeof(int));
        for (j=0; j<nseqs[i]; j++){
            fscanf(fp,"%s",filen2);
            getfeatures(filen2, data[i]+j, len[i]+j, &fdim, dim);
        }
        fclose(fp);
    }
    fclose(fctl);

    return;

}

void  getfeatures(char *cepfile, float ***feat, int *nfr, int *featdim, 
                  int cepdim)
{
    // Basically a stub routine
    *featdim = cepdim;
    if (read_cep(cepfile,feat,nfr,cepdim) == IO_ERR)
                QUIT(("Unable to read %s\n",cepfile))
    return;
}
