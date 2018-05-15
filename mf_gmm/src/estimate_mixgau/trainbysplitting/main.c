/******************************************************************************
*    Main routine for computing Gaussian Mixture Distributions
*    Coded by Bhiksha Raj, June 94
******************************************************************************/
#include <stdlib.h>
#include <string.h>
#include "EMheader.h"

#define DESCRIPTION "\nThis program computes a mixture gaussian distribution for a set of speech files.\n"

#define USAGE "USAGE:\n%s \t-c <control file> \\\n\t-o <output file> \\\n\t-d <dimensionality of data> \\\n\t-m <Total no of modes in distribution>\\\n\t-z the file containing the initial distribution to continue EM from> \n"

int main(int argc, char **argv)
{
    REAL   **vector;
    float **buff; // Read routines read float
    REAL   **mean,**variance,*c;
    int     numspch,numvecs,Ndim,Nmodes,curnmodes,maxlength=0;
    int     i,j,length,initialize,superiter;
    char    ctlfile[512],outfile[512];
    char    filename[512],initialcbk[512], sw;

    FILE    *fopen(),*listfile;

    /*
     * by default I assume the programs is to be run from data
     * only with no initial set of Gaussians. i.e. go to vq and them 
     * em.
     */
    /* SET DEFAULTS */
    initialize = 1 ; 

    /* Print out command, just to be nice */
    for (i=0; i<argc; i++) printf("%s ",argv[i]); printf("\n"); fflush(stdout);
    if (argc == 1) {
        printf(DESCRIPTION);
        QUIT ((USAGE,argv[0]));
    }
    for (i=1;i<argc;++i) 
    {
        if (argv[i][0]!='-')
            QUIT ((USAGE,argv[0]));
        sw = argv[i][1];
        switch(sw) 
        {
            case 'c': case 'C': 
            {
                strcpy(ctlfile,argv[++i]);
                break;
            }
            case 'o': case 'O': 
            {
                strcpy(outfile,argv[++i]);
                break;
            }
            case 'd': case 'D': 
            {
                Ndim = atoi(argv[++i]);
                break;
            }
            case 'm': case 'M': 
            {
                Nmodes = atoi(argv[++i]);
                break;
            }
	    case 'z': case 'Z': 
            {
		strcpy(initialcbk,argv[++i]);
		initialize = 0 ; 
		break;
            }
            default:
            {
                printf(DESCRIPTION);
                QUIT ((USAGE,argv[0]));
            }

        }
    }


    listfile = fopen(ctlfile,"r");
    if (listfile == NULL) 
        QUIT(("Unable to open control file %s\n",ctlfile));

    numvecs = 0; maxlength = 0;
    while (fscanf(listfile,"%s\n",filename) != EOF)
    {
#ifdef DEBUG
        printf("Reading %s\n",filename);
#endif
        if (read_cep(filename, &buff, &length, Ndim) == -1)
        {
        /*
         * What happens on a datafile read error? We can either just go on to
         * the next cepfile, or stop the program.  Set switch in header.h.
         */
#ifdef IGNORE_READ_ERR
            printf("READ FAILURE; skipping cepfile\n");
            fflush(stdout);
#else
            QUIT(("Unable to read %s\n",filename));
#endif
        }
        else
        {
            maxlength += length;
	    free2d((void**)buff);
        }
    }
    printf("%d vectors in all\n",maxlength);
    rewind(listfile);

    if((vector = (REAL **) alloc2d(maxlength,Ndim,sizeof(REAL))) == NULL)
        QUIT(("Unable to Allocate Space for Vector Array\n"));

    while ((fscanf(listfile,"%s\n",filename) != EOF) && (numvecs <= maxlength))
    {
        if (read_cep(filename, &buff, &length, Ndim) == -1)
        {
        /*
         * What happens on a datafile read error? We can either just go on to
         * the next cepfile, or stop the program.  Set switch in header.h.
         */
#ifdef IGNORE_READ_ERR
            printf("READ FAILURE; skipping cepfile\n");
            fflush(stdout);
#else
            QUIT(("Unable to read %s\n",filename));
#endif
        }
        else
        {
#ifdef DEBUG
	    printf("Success in reading %s feature file\n", filename);
#endif
            for (i=0;i<length;i++)
	    {
		for (j=0;j<Ndim;++j)
		    vector[numvecs][j] = (REAL) buff[i][j];
                ++numvecs;
	    }
	    free2d((void**)buff);
        }
        if (numvecs > maxlength)
            QUIT(("**** Too many frames? Check your dimensionality!! ****\n"));
    }
    fclose(listfile);

    if (numvecs == 0)
        QUIT(("Only 0 vectors to compute a DISTRIBUTION!\nI am quitting!\n"));

    numspch = numvecs;

   
    /*
     * do this only if we are not requesting a restart from a previous 
     * temp statistics file .
     */ 
    if (initialize) 
    { 
       /*
        * allocate the mean and variance and c arrays.
        */ 
        c = (REAL *) malloc(Nmodes*sizeof(REAL));
        mean = (REAL **) alloc2d(Nmodes,Ndim,sizeof(REAL));
        if((variance = (REAL **) alloc2d(Nmodes,Ndim,sizeof(REAL))) == NULL)
            QUIT(("Unable to allocate space for variances\n"));
        curnmodes = 0;
    }
    else 
    {
       /* Read the distribution the input file */
	if (read_distribution(initialcbk,&Nmodes,&Ndim,&c,&mean,&variance)!=0) 
            QUIT(("Unable to read initial distribution\n"));
        curnmodes = Nmodes;
    }
	
    estimate_mixture_gaussian(vector,numspch,Ndim,Nmodes,mean,variance,c,
                              outfile,curnmodes);
    if (store_distribution(outfile,Nmodes,Ndim,c,mean,variance) != 0) 
    {
        printf("Unable to open %s to store distribution\n",outfile);
        printf("Superiter = %d\n",superiter);
    }

    free2d((void**)vector);
    free(c);
    free2d((void**)mean);
    free2d((void**)variance);
    return 0;
}
