/******************************************************************************
*    Main routine for computing the distribution for CDCN.
*    Coded by Bhiksha Raj, June 94
******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "header.h"

#define DESCRIPTION "\nThis program computes stereo corrections for a GMM\n"

#define USAGE "USAGE:\n%s -c <control file> \\\n\t-x <indir1> \\\n\t-y <indir2> \\\n\t-o <output file> \\\n\t-l <historycontext> \\\n\t-r <futurecontext> \\\n\t-i the file containing the initial distribution to continue EM from>\n"

int main(int argc, char **argv)
{
    float   **vectorx,**vectory;
    float   **mean,**variance,*c;
    float   ***meany, *count;
    int     numvecs,Ndim,Nmodes,nvecy,totalvecs;
    int     i,j,k, leftctx, rightctx;
    char    *ctlfile=NULL,*indirx=NULL,*indiry=NULL,*inputgmm=NULL,*outfile=NULL;
    char    basefile[512],filename[512],sw;

    FILE    *fopen(),*listfile;

    /* Print out command, just to be nice */
    for (i=0; i<argc; i++) printf("%s ",argv[i]); printf("\n"); fflush(stdout);
    if (argc == 1) {
        printf(DESCRIPTION);
        QUIT ((USAGE,argv[0]));
    }
    for (i=1;i<argc;++i) {
        if (argv[i][0]!='-')
            QUIT ((USAGE,argv[0]));
        sw = argv[i][1];
        switch(sw) {
            case 'c': case 'C': {
                ctlfile = argv[++i];
                break;
            }
            case 'x': {
                indirx = argv[++i];
                break;
            }
            case 'y': {
                indiry = argv[++i];
                break;
            }
            case 'o': case 'O': {
                outfile = argv[++i];
                break;
            }
	    case 'i': {
                inputgmm = argv[++i];
		break;
            }
            case 'l': {
                sscanf(argv[++i],"%d",&leftctx);
                break;
            }
            case 'r': {
                sscanf(argv[++i],"%d",&rightctx);
                break;
            }
            default: {
                printf(DESCRIPTION);
                QUIT ((USAGE,argv[0]));
            }

        }
    }


    listfile = fopen(ctlfile,"r");
    if (listfile == NULL) 
        QUIT(("Unable to open control file %s\n",ctlfile));

    /* Read input distribution */
    if (read_input_distribution(inputgmm,&mean, &variance,&c,&Nmodes,&Ndim)) {
        QUIT(("Error reading %s\b",inputgmm))
    }
    meany = (float ***) alloc3d(leftctx+rightctx+1,Nmodes,Ndim,sizeof(float));
    count = (float *) calloc(Nmodes,sizeof(float));

    totalvecs = 0;
    while (fscanf(listfile,"%s\n",basefile) != EOF) {
        sprintf(filename,"%s/%s",indirx,basefile);
        if (read_cep(filename,&vectorx,&numvecs,Ndim) == -1) {
            QUIT(("Unable to read %s\n",filename));
        }
        sprintf(filename,"%s/%s",indiry,basefile);
        if (read_cep(filename,&vectory,&nvecy,Ndim) == -1) {
            QUIT(("Unable to read %s\n",filename));
        }
        if (nvecy != numvecs) {
            QUIT(("Clean and stereo data do not match\n"));
        }
        estimate_gmm_stereo_corrections(vectorx, vectory,
            numvecs, Ndim, Nmodes, mean, variance, c, rightctx, leftctx,
            meany, count);
        totalvecs += numvecs;

        free2d((void**)vectorx); free2d((void**)vectory);
    }
    printf("%d vectors in all\n",totalvecs);
    fclose(listfile);
	
    // Normalize
    for (i=-leftctx; i<=rightctx; i++) {
        for (j=0; j<Nmodes; j++) {
            for (k=0; k<Ndim; k++) meany[i+leftctx][j][k] /= count[j];
        }
    }
    store_distribution(outfile,Nmodes,Ndim,c,mean,variance,
                       leftctx,rightctx,meany);

    free (c);
    free2d((void**)mean);
    free2d((void**)variance);
    free3d((void***)meany); free(count);
    return 0;
}
