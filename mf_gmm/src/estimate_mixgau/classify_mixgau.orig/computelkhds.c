/*---------------------------------------------------------------------------
* Program that actually trains a gaussian mixture distribution using EM
* given a control file and the parameters (no. of mixture components needed)
* The mfc file is read, and a combination of mfc and delta-mfc features is
* obtained and the joint feature is used for the distribution estimation.
*----------------------------------------------------------------------------*/
#include <ayesha.h>

void  getfeatures(char*,float****,int**,int*,int*);
void  estimate_mixture_gau(float**,int,int,int,float**,float**,float*,int);
void read_classes(char*,char***,int*,int**,int*,float****,float****,float***);
float likelihood(float**,int,float*,float**,float**,int,int);
float gauss(float*,float*,float*,float,int);

int main(int argc, char **argv)
{
    float ***feat,**mixw,***mean,***vars,lkhd;
    int  nseg,*nframes,featdim,i,j,*ngau,nclass,dim;
    char *ctlfile,*outfile,*classfile,**class;
    FILE *fout;

    if (argc != 4) 
        QUIT(("USAGE:\n %s <ctlfile> <classfile> <outfile>\n",argv[0]))
    ctlfile = argv[1]; classfile = argv[2]; outfile=argv[3];

    if ((fout=fopen(outfile,"w")) == NULL)
        QUIT(("Unable to write %s\n",outfile))
    getfeatures(ctlfile,&feat,&nframes,&nseg,&featdim); 

    read_classes(classfile,&class,&nclass,&ngau,&dim,&mean,&vars,&mixw);
    if (dim != featdim)
        QUIT(("Mismatch in featdim %d and dist dim %d\n",featdim,dim))
    for (i=0;i<nseg;i++){
        for (j=0;j<nclass;j++){
          lkhd = likelihood(feat[i],nframes[i],mixw[j],
                                   mean[j],vars[j],ngau[j],dim);
          fprintf(fout,"%f ",lkhd);
        }
        fprintf(fout,"\n");
    }
    fclose(fout);
    return 0;
}


/*==========================================================================*/
/* Routine to read mfc file and get feature mfc+deltamfc */

#define CEPDIM 13
void  getfeatures(char *ctlfile, float ****feat, int **nfr, 
                                                    int *numseg, int *featdim)
{
    float **cep, **features=NULL, ***outfeat, cepmean[100];
    int i,j,k,beg,end,frmno,*nfrm,cepdim,outdim,totframes=0;
    int filenfrm,maxframes=0,nseg;
    char fname[256],prevfile[256],*dummy="/dev/null";
    FILE *fp;

    cepdim = (int)CEPDIM; outdim = cepdim;
    if ((fp=fopen(ctlfile,"r"))==NULL) QUIT(("Unable to open %s\n",ctlfile))

    nseg = 0;
    while(fscanf(fp,"%s %d %d",fname,&beg,&end)!=EOF) {
        totframes = end - beg + 1; 
        if (maxframes < totframes) maxframes = totframes; 
        nseg++;
    }
    rewind(fp);

    *numseg = nseg;
    *feat = outfeat = (float ***)alloc3d(nseg,maxframes,outdim,sizeof(float));
    *nfr = nfrm = (int *)calloc(nseg,sizeof(int));
    *featdim = outdim;

    strcpy(prevfile,dummy); nseg = 0;
    while(fscanf(fp,"%s %d %d",fname,&beg,&end)!=EOF) {
        if (strcmp(fname,prevfile) != 0) {
            printf("Reading from %s\n",fname); fflush(stdout);
            if (features != NULL) free2d((void **)features);
            if (read_cep(fname,&cep,&filenfrm,cepdim) == IO_ERR)
                QUIT(("Unable to read %s\n",fname))
            features = (float **)alloc2d(filenfrm,2*cepdim,sizeof(float));

            // for (i=2;i<filenfrm-2;i++){
            for (i=0;i<filenfrm;i++){
                for(j=0;j<cepdim;j++) features[i][j] = cep[i][j];
//                for(k=0;k<cepdim;k++,j++) {
//                    features[i][j] = 0.67*cep[i+1][k]+0.33*cep[i+2][k] -
//                             0.67*cep[i-1][k]-0.33*cep[i-2][k];
//                }
            }
//            for(j=0;j<cepdim;j++) {
//                features[0][j] = cep[0][j];
//                features[1][j] = cep[1][j];
//                features[filenfrm-1][j] = cep[filenfrm-1][j];
//                features[filenfrm-2][j] = cep[filenfrm-2][j];
//                features[0][j+cepdim]=features[2][j+cepdim];
//                features[1][j+cepdim]=features[2][j+cepdim];
//                features[filenfrm-1][j+cepdim]=features[filenfrm-3][j+cepdim];
//              features[filenfrm-2][j+cepdim]=features[filenfrm-3][j+cepdim];
 //           }
            free2d((void **)cep);
        }
        for (j=0;j<outdim;j++) cepmean[j] = 0;
        if (end < beg || end > filenfrm)
           QUIT(("Beg %d, End %d and nfr %d error %s\n",beg,end,filenfrm,fname))
        for (i=beg,k=0;i<=end;i++,k++){
            for (j=0;j<outdim;j++) outfeat[nseg][k][j] = features[i][j];
//            for (j=0;j<outdim;j++) cepmean[j] += features[i][j];
        }
//        for (j=0;j<outdim;j++) cepmean[j] /= (float)(end-beg+1);
//        for (i=beg,k=0;i<=end;i++,k++) {
//            for (j=0;j<cepdim;j++) outfeat[nseg][k][j] -= cepmean[j];
//        }
	nfrm[nseg] = end-beg+1; nseg++;
        strcpy(prevfile,fname);
    }
    if (features != NULL) free2d((void **)features);
    printf("%d of %d total frames obtained from cepfiles\n",filenfrm,totframes);

    return;
}

/*=========================================================================*/
float likelihood(float **x, int nframes, float *corprod, 
                 float **mean, float **invvar, int ngau, int Ndim)
{
    float lkhd,Tau[256],max,sum;
    int i,j;

    lkhd = 0;
    for (i=0;i<nframes;++i)
    {
        max = Tau[0] = gauss(x[i],mean[0],invvar[0],corprod[0],Ndim);
        for (j=1; j<ngau; j++){
             Tau[j] = gauss(x[i],mean[j],invvar[j],corprod[j],Ndim);
             if (Tau[j] > max) max = Tau[j];
        }
        sum = 0.0;
        for (j=0;j<ngau;++j) {
            if (Tau[j] > max-15) Tau[j] = exp(Tau[j] - max + 10);
            else Tau[j] = 0;
            sum += Tau[j];
        }
        lkhd += max + log(sum) - 10;
    }
    return(lkhd);
}

/*---------------------------------------------------------------------------*/

float gauss(float *x,float *mean,float *invvar,float corprod,int Ndim)
{
    auto int	j;
    auto float	power=0,density;

    for (j=0;j<Ndim;++j)
	power += (x[j]-mean[j])*(x[j]-mean[j])*invvar[j];
    power /= -2.0;
    density = corprod + power;
    return(density);
}
/*---------------------------------------------------------------------------*/

// Read class info from the given file. Kinda long..
void read_classes(char *classfile, char ***class, int *nclass, 
                  int **ngau, int *dim, 
                  float ****mean, float ****var, float ***mixw)
{
    float ***tmean,***tvar,**tmixw;
    int   ncls,*tngau,ndim,pdim=-1;
    char  **clsname;
    char  modelfile[256],line[1024];
    int i,j,clsid = 0;

    FILE *fp, *fp2;

    // Ignore all lines beginning with '#' and find how many classes
    if ((fp=fopen(classfile,"r")) == NULL)
        QUIT(("Unable to open %s for reading\n",classfile))
    fgets(line,256,fp);
    while(line[0]=='#') fgets(line,256,fp);
    for(ncls=1;fgets(line,256,fp)!=NULL;ncls++);
    rewind(fp);
    fgets(line,256,fp);
    while(line[0]=='#') fgets(line,256,fp);

    // Preallocate outermost pointer for multidimensional arrays
    clsname = (char **)alloc2d(ncls,100,sizeof(char));
    tngau = (int *)calloc(ncls,sizeof(int));
    tmean = (float ***)calloc(ncls,sizeof(float **));
    tvar = (float ***)calloc(ncls,sizeof(float **));
    tmixw = (float **)calloc(ncls,sizeof(float *));

    // Since we already have a line, use a do-until loop.
    do{
        sscanf(line,"%s %s",clsname[clsid],modelfile);
        printf("Reading models for %s\n",clsname[clsid]);
        if ((fp2=fopen(modelfile,"r")) == NULL)
            QUIT(("Unable to open %s for reading\n",modelfile))
        fscanf(fp2,"%d %d",&tngau[clsid],&ndim);
        if (pdim != -1 && pdim != ndim) 
            QUIT(("Dimension mismatch in some models\n"))
        pdim = ndim;
        tmean[clsid] = (float **) alloc2d(tngau[clsid],ndim,sizeof(float));
        tvar[clsid] = (float **) alloc2d(tngau[clsid],ndim,sizeof(float));
        tmixw[clsid] = (float *) calloc(tngau[clsid],sizeof(float));
        for (i=0; i<tngau[clsid]; i++){
            fscanf(fp2,"%f",&tmixw[clsid][i]);
            for (j=0;j<ndim;j++) fscanf(fp2,"%f",&tmean[clsid][i][j]);
            for (j=0;j<ndim;j++) fscanf(fp2,"%f",&tvar[clsid][i][j]);
            tmixw[clsid][i] = log((double)tmixw[clsid][i]);
            for (j=0;j<ndim;j++) {
                if (tvar[clsid][i][j] < 1.0e-20) tvar[clsid][i][j] = 1.0e-20;
                tmixw[clsid][i] -= 0.5*log((double)tvar[clsid][i][j]);
            }
            for (j=0;j<ndim;j++) tvar[clsid][i][j] = 1.0/tvar[clsid][i][j];
        }
        fclose(fp2);
        clsid++;
    } while(fgets(line,256,fp)!=NULL);
    fclose(fp);

    *mean = tmean; *var = tvar; *mixw = tmixw;
    *ngau = tngau; *nclass = ncls; *dim = ndim; *class = clsname;
    return;
}
