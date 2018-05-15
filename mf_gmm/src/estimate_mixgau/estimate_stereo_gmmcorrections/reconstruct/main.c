/******************************************************************************
*    Main routine for computing the distribution for CDCN.
*    Coded by Bhiksha Raj, June 94
******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "header.h"

#define DESCRIPTION "\nThis program computes stereo corrections to a noisy file using stereo GMMs\n"

#define USAGE "USAGE:\n%s -c <control file>\\\n\t-x <indir>\\\n\t-o <outdir>\\\n\t-i the file containing the distribution>\n"

int main(int argc, char **argv)
{
    float   **vectorx,***vectory;
    float   **mean,**variance,*c;
    float   ***meany;
    int     numvecs,Ndim,Nmodes;
    int     i, leftctx, rightctx;
    char    *ctlfile=NULL,*indirx=NULL,*inputgmm=NULL,*outdir=NULL;
    char    basefile[512],filename[512],sw;

    FILE    *fopen(),*listfile;

    /* Print out command, just to be nice */
    for (i=0; i<argc; i++) printf("%s ",argv[i]); printf("\n"); fflush(stdout);
    if (argc == 1) {
        printf(DESCRIPTION);
        QUIT ((USAGE,argv[0]));
    }
    for (i=1;i<argc;++i) {
        if (argv[i][0]!='-') QUIT ((USAGE,argv[0]));
        sw = argv[i][1];
        switch(sw) {
            case 'c': { ctlfile = argv[++i]; break; }
            case 'x': { indirx = argv[++i]; break; }
	    case 'i': { inputgmm = argv[++i]; break; }
            case 'o': { outdir = argv[++i]; break; }
            default: { printf(DESCRIPTION); QUIT ((USAGE,argv[0])); }
        }
    }

    listfile = fopen(ctlfile,"r");
    if (listfile == NULL) 
        QUIT(("Unable to open control file %s\n",ctlfile));

    /* Read input distribution */
    if (read_input_distribution(inputgmm,&mean,&variance,&c,&Nmodes,&Ndim,&leftctx,&rightctx,&meany)) {
        QUIT(("Error reading %s\b",inputgmm))
    }

    while (fscanf(listfile,"%s\n",basefile) != EOF) {
        sprintf(filename,"%s/%s",indirx,basefile);
#ifdef DEBUG
        printf("Reading %s\n",filename);
#endif
        if (read_cep(filename,&vectorx,&numvecs,Ndim) == -1) {
            QUIT(("Unable to read %s\n",filename));
        }
//        sprintf(filename,"%s/%s",outdir,basefile);
        sprintf(filename,"%s",outdir);

        vectory = (float***) alloc3d(leftctx+rightctx+1, numvecs,
                                     Ndim, sizeof(float)); 
        estimate_gmm_stereo_corrections(vectorx, vectory,
            numvecs, Ndim, Nmodes, mean, variance, c, rightctx, leftctx,
            meany);

        write_file(leftctx+rightctx+1,numvecs,Ndim,vectory,filename);
        free2d((void**)vectorx); free3d((void***)vectory);
    }
    fclose(listfile);
	
    free (c);
    free2d((void**)mean);
    free2d((void**)variance);
    free3d((void***)meany);
    return 0;
}
