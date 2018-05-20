#include "trainclass.h"

#include <sys/stat.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/file.h>
#include <sys/uio.h>

#define TRUE 1
#define FALSE 0

/* Macro to byteswap an int variable.  x = ptr to variable */
#define SWAP_INT(x)   *(x) = ((0x000000ff & (*(x))>>24) | \
                                (0x0000ff00 & (*(x))>>8) | \
                                (0x00ff0000 & (*(x))<<8) | \
                                (0xff000000 & (*(x))<<24))
/* Macro to byteswap a float variable.  x = ptr to variable */
#define SWAP_FLOAT(x) SWAP_INT((int *) x)

int read_cep(char *file, float***cep, int *numframes, int cepsize)
{
    FILE *fp;
    int n_float;
    struct stat statbuf;
    int i, n, byterev, sf, ef;
    float **mfcbuf;

    if (stat(file, &statbuf) < 0) {
        printf("stat(%s) failed\n", file);
        return -1;
    }

    if ((fp = fopen(file, "rb")) == NULL) {
	printf("fopen(%s,rb) failed\n", file);
	return -1;
    }
    
    /* Read #floats in header */
    if (fread(&n_float, sizeof(int), 1, fp) != 1) {
	fclose (fp);
	return -1;
    }
    
    /* Check if n_float matches file size */
    byterev = FALSE;
    if ((n_float*sizeof(float) + 4) != statbuf.st_size) { //st_size size in bytes
    	n = n_float;
    	SWAP_INT(&n);

    	if ((n*sizeof(float) + 4) != statbuf.st_size) {
    	    printf("Header size field: %d(%08x); filesize: %d(%08x)\n",
    		    n_float,n_float,(int)statbuf.st_size,(int)statbuf.st_size);
    	    fclose (fp);
    	    return -1;
    	}
    	n_float = n;
    	byterev = TRUE;
    }
    if (n_float <= 0) {
    	printf("Header size field: %d\n",  n_float);
    	fclose (fp);
    	return -1;
    }
    
    /* n = #frames of input */
    n = n_float/cepsize;
    if (n * cepsize != n_float) {
    	printf("Header size field: %d; not multiple of %d\n", n_float, cepsize);
    	fclose (fp);
    	return -1;
    }
    sf = 0;
    ef = n;

    mfcbuf = (float **) alloc2d (n, cepsize, sizeof(float));
    
    /* Read mfc data and byteswap if necessary */
    n_float = n * cepsize;
    if (fread (mfcbuf[0], sizeof(float), n_float, fp) != n_float) {
	printf("Error reading mfc data\n");
	fclose (fp);
	return -1;
    }
    if (byterev) {
	for (i = 0; i < n_float; i++)
	    SWAP_FLOAT(&(mfcbuf[0][i]));
    }
    fclose (fp);

    *numframes = n;
    *cep = mfcbuf;
    return (1);
}

int readcep_part(char *file, int sf, int ef, float***cep, 
                                    int *numframes, int cepsize)
{
    FILE *fp;
    int n_float;
    struct stat statbuf;
    int i, n, byterev;
    float **mfcbuf;

    if (stat(file, &statbuf) < 0) {
        printf("stat(%s) failed\n", file);
        return -1;
    }

    if ((fp = fopen(file, "rb")) == NULL) {
	printf("fopen(%s,rb) failed\n", file);
	return -1;
    }
    
    /* Read #floats in header */
    if (fread(&n_float, sizeof(int), 1, fp) != 1) {
	fclose (fp);
	return -1;
    }
    
    /* Check if n_float matches file size */
    byterev = FALSE;
    if ((n_float*sizeof(float) + 4) != statbuf.st_size) {
	n = n_float;
	SWAP_INT(&n);

	if ((n*sizeof(float) + 4) != statbuf.st_size) {
	    printf("Header size field: %d(%08x); filesize: %d(%08x)\n",
		    n_float,n_float,(int)statbuf.st_size,(int)statbuf.st_size);
	    fclose (fp);
	    return -1;
	}

	n_float = n;
	byterev = TRUE;
    }
    if (n_float <= 0) {
	printf("Header size field: %d\n",  n_float);
	fclose (fp);
	return -1;
    }
    
    /* n = #frames of input */
    n = n_float/cepsize;
    if (n * cepsize != n_float) {
	printf("Header size field: %d; not multiple of %d\n", n_float, cepsize);
	fclose (fp);
	return -1;
    }

    mfcbuf = (float **) alloc2d (n, cepsize, sizeof(float));
    
    /* Read mfc data and byteswap if necessary */
    n_float = n * cepsize;
    if (fread (mfcbuf[0], sizeof(float), n_float, fp) != n_float) {
	printf("Error reading mfc data\n");
	fclose (fp);
	return -1;
    }
    if (byterev) {
	for (i = 0; i < n_float; i++)
	    SWAP_FLOAT(&(mfcbuf[0][i]));
    }
    fclose (fp);

    *numframes = ef - sf + 1;
    *cep = (float **) alloc2d(*numframes,cepsize,sizeof(float));
    { int j;
      for (i = 0; i < *numframes; i++) 
          for (j = 0; j < cepsize; j++)
              *cep[i][j] = mfcbuf[sf+i][j];
    }

    free2d((void **)mfcbuf);
    return (1);
}
void  getfeatures(char *ctlfile, float ***feat, int *nfr, int *featdim, int cepdim)
{
    float **cep, **features=NULL, **outfeat, cepmean[100];
    int i,j,k,beg,end,frmno,nfrm,outdim,totframes=0,begseg;
    char fname[256],prevfile[256],*dummy="/dev/null";
    FILE *fp;

    outdim = DELTA? 2*cepdim : cepdim;

    if ((fp=fopen(ctlfile,"r"))==NULL) QUIT(("Unable to open %s\n",ctlfile))
    while(fscanf(fp,"%s %d %d",fname,&beg,&end)!=EOF)
        totframes += end - beg + 1;
    rewind(fp);

    outfeat = (float **)alloc2d(totframes,outdim,sizeof(float));
    *featdim = outdim; *feat = outfeat;

    strcpy(prevfile,dummy); frmno = 0;
    while(fscanf(fp,"%s %d %d",fname,&beg,&end)!=EOF){
        if (beg < 1) beg = 1;
        if (strcmp(fname,prevfile) != 0){
            printf("Reading from %s\n",fname); fflush(stdout);
            read_cep(fname,&cep,&nfrm,cepdim) == IO_ERR;
}

int main(int argc, char **argv)
{
    float **feat,*mixw,**mean,**vars;
    int  nframes,featdim,nmodes,cepdim;
    char *ctlfile,*outfile;

    if (argc != 4 && argc != 5) QUIT(("USAGE:\n %s <ctlfile> <nmodes> <outfile> [featdim(def13)]\n",argv[0])){
    ctlfile = argv[1]; sscanf(argv[2],"%d",&nmodes);outfile=argv[3];}

    cepdim = argc == 5 ? atoi(argv[4]) : 13;

    getfeatures(ctlfile,&feat,&nframes,&featdim,cepdim);
    return(1);
}
