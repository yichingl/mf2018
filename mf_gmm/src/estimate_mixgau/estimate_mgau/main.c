/******************************************************************************
*    Main routine for computing the distribution for CDCN.
*    Coded by Bhiksha Raj, June 94
******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "header.h"

#define DESCRIPTION "\nThis program computes a mixture gaussian distribution for a set of speech\ncepstral files. The dimensionality of the cepstral files is assumed to be 13\n"

#define USAGE "USAGE:\n%s -e <filename extension (def none)> \\\n\t-c <control file> \\\n\t-i <indir (def none)> \\\n\t-o <output file> \\\n\t-m <Total no of modes in distribution>\\\n\t-d <dimensionality of the data (default 13)> \\\n\t-n <For CMN> \\\n\t-l <To include delta features> \\\n\t-z the file containing the initial distribution to continue EM from> \\\n\t-t <file to dump temporary results to> \\\n\t-r <fraction of lowest-0th-term vectors to discard (def 0)>\n"

int computefeat(float ***buff, int CMS, int DELTA, float FILTERLOWEGY,
                int Ndim, int *length)
{
    int i, j;
    float **lbuff = *buff;

    if (FILTERLOWEGY > 0) { /* Remove lowest egy vectors within 5% of min egy */
        float min = 1e+30;
        float max = -1e+30;
        float c0thresh;
        int k;
        float **ltmp = *buff;

        for (i=0; i<*length; i++) {
            if (ltmp[i][0] < min) min = ltmp[i][0];
            if (ltmp[i][0] > max) max = ltmp[i][0];
        }
        c0thresh = min + (max - min)*FILTERLOWEGY;
        for (i=0,k=0; i<*length; i++) {
            if (ltmp[i][0] < c0thresh) continue;
            for (j=0; j<Ndim; j++) ltmp[k][j] = ltmp[i][j];
            k++;
        }
        *length = k;
    }

    if (CMS) {
        float *mean;
        mean = (float *)calloc(Ndim,sizeof(float));
        for (i=0; i<*length; i++)
            for (j=0; j<Ndim; j++) mean[j] += lbuff[i][j];
        for (j=0; j<Ndim; j++) mean[j] /= (float)(*length);
        for (i=0; i<*length; i++)
            for (j=0; j<Ndim; j++) lbuff[i][j] -= mean[j];
        free(mean);
    }
    if (DELTA) {
        float **tmp;
        *length -= 4;
        *buff = tmp = (float **) alloc2d(*length,2*Ndim,sizeof(float));
        for (i=0; i<*length; i++) {
            for (j=0; j<Ndim; j++) tmp[i][j] = lbuff[i+2][j];
            for (j=0; j<Ndim; j++)
                tmp[i][j+Ndim] = lbuff[i+4][j] - lbuff[i][j];
        }
        free2d((void**)lbuff);
    }

    return 0;
}

int main(int argc, char **argv)
{
    float   **vector,**buff;
    float   **mean,**variance,*c;
    float   atemp, FILTERLOWEGY=0;
    int     numspch,numvecs,Ndim,Nmodes,maxlength=0,CMS=0, DELTA=0;
    int     i,j,k,size,length,*nbin,*bin;
    int	    initialize,superiter;
    char    ctlfile[512],*infileext=NULL,outfile[512],tempfile[512];
    char    basefile[512],filename[512],initialcodebk[512], sw, *indir=NULL;

    FILE    *fopen(),*listfile;
    float   scale();
    char    **alloc2d();

    void    estimate_multi_modals(),vector_quantize();
    int    read_backup_distribution(),store_distribution();


    /*
     * by default I assume the programs is to be run from data
     * only with no initial set of Gaussians. i.e. go to vq and them 
     * em.
     */
    /* SET DEFAULTS */
    initialize = 1 ; 
    strcpy (tempfile,"CDCN.DIST.TEMP");
    Ndim = 13;

	/* Print out command, just to be nice */
	for (i=0; i<argc; i++) printf("%s ",argv[i]); printf("\n"); fflush(stdout);
    if (argc == 1)
    {
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
            case 'e': case 'E':
            {
                infileext = argv[++i];
                break;
            }
            case 'i':
            {
                indir = argv[++i];
                break;
            }
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
            case 'm': case 'M':
            {
                Nmodes = atoi(argv[++i]);
                break;
            }
            case 'd': case 'D':
            {
                Ndim = atoi(argv[++i]);
                break;
            }
            case 'n':
            {
                CMS = 1;
                break;
            }
            case 'l':
            {
                DELTA = 1;
                break;
            }
            case 'r': 
            {
                sscanf(argv[++i],"%f",&FILTERLOWEGY);
                break;
            }
	    case 'z': case 'Z':
	    {
		strcpy(initialcodebk,argv[++i]);
		/*
  		 * do not initialize from VQ but from this file
		 */

		/*
		 * notice that this is ONLY prepared for restarting the
		 * program from a previous EM files with the same number
		 * of Gaussians. i.e. if the machinbe breaks you can 
		 * contine EM from there
		 */
		initialize = 0 ; 
		break;
            }
	    case 't': case 'T':
	    {
               /*
                * Temporary file to store partially converged distribution
                * to
                */
		strcpy(tempfile,argv[++i]);
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
    while (fscanf(listfile,"%s\n",basefile) != EOF)
    {
        /*
         * If an extension has been specified for the files, append it
         */
        if (infileext) {
            if (indir) sprintf(filename,"%s/%s.%s",indir,basefile,infileext);
            else sprintf(filename,"%s.%s",basefile,infileext);
        }
        else {
            if (indir) sprintf(filename,"%s/%s",indir,basefile);
            else sprintf(filename,"%s",basefile);
        }

        if (read_cep(filename,&buff,&length,Ndim) == -1)
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
            if (!DELTA)
                maxlength += length;
            else
                maxlength += length - 4;
	    free2d((void**)buff);
        }
    }
    printf("%d vectors in all\n",maxlength);
    rewind(listfile);

    /*
     * For our particular case it is needed that the Maximum dimensionality
     * is the same as the actual dimension of the data. Else the read_cep
     * routine fails!
     */
    if (DELTA) 
        vector = (float **) alloc2d(maxlength,2*Ndim,sizeof(float));
    else
        vector = (float **) alloc2d(maxlength,Ndim,sizeof(float));
    if(!vector) QUIT(("Unable to Allocate Space for Vector Array\n"));

    while ((fscanf(listfile,"%s\n",basefile) != EOF) && (numvecs <= maxlength))
    {
        /*
         * If an extension has been specified for the files, append it
         */
        if (infileext) {
            if (indir) sprintf(filename,"%s/%s.%s",indir,basefile,infileext);
            else sprintf(filename,"%s.%s",basefile,infileext);
        }
        else {
            if (indir) sprintf(filename,"%s/%s",indir,basefile);
            else sprintf(filename,"%s",basefile);
        }
        if (read_cep(filename,&buff,&length,Ndim) == -1)
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
            computefeat(&buff,CMS,DELTA,FILTERLOWEGY,Ndim,&length);
            for (i=0;i<length;i++)
	    {
		for (j=0;j<Ndim;++j)
		    vector[numvecs][j] = buff[i][j];
                if (DELTA)
		    for (j=0;j<Ndim;++j)
		        vector[numvecs][j+Ndim] = buff[i][j+Ndim];
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
    if (DELTA) Ndim *= 2;

    /*
     * We Vector Quantize to obtain the initial values for the EM.
     * If this codebook already exists, we skip the VQ and directly
     * compute the variances and c[]s after obtaining the mean values
     * as the code words in the existing codebook
     */
   
    /*
     * do this only if we are not requesting a restart from a previous 
     * temp statistics file .
     */ 
    if (initialize) 
    { 
       /*
        * allocate the mean and variance and c arrays.
        */ 
        c = (float *) malloc(Nmodes*sizeof(float));
        mean = (float **) alloc2d(Nmodes,Ndim,sizeof(float));
        if((variance = (float **) alloc2d(Nmodes,Ndim,sizeof(float))) == NULL)
            QUIT(("Unable to allocate space for variances\n"));

       nbin = (int *)malloc(Nmodes*sizeof(int)); /* no of vectors in a mode */
       if (nbin==NULL)
           QUIT(("Unable to allocate %d space for nbin array\n",Nmodes));

      /*
       * The vector_quantize routine performs VQ with a mahalonobis metric
       * and returns the codes as the means and the wieghts as the variances
       * of the initial estimates of the modes, which will further be 
       * employed in EM. Note that the variances are purely diagonal
       * We initialize all initial c[] to be equal
       */

       size = numspch*sizeof(int);
       if ((bin = (int *)malloc(size)) == NULL)
           QUIT(("Unable to allocate %d space for bin array\n",size));


       vector_quantize(mean,Nmodes,vector,numspch,Ndim,bin);

       for (i=0;i<Nmodes;++i)
           c[i] = 1.0/(float)Nmodes;
       for (k=0;k<Nmodes;++k)
       {
           nbin[k] = 0;
           for (i=0;i<Ndim;++i)
               variance[k][i] = 0;
       }
       for (i=0;i<numspch;++i)
       {
           for (j=0;j<Ndim;++j)
           {
               atemp = (vector[i][j] - mean[bin[i]][j]);
               variance[bin[i]][j] += atemp*atemp;
           }
           ++nbin[bin[i]];
       }
       for (k=0;k<Nmodes;++k)
       {
           for(j=0;j<Ndim;++j)
               variance[k][j] /= nbin[k];
       }

       free(bin);             /* We do not need this array anymore */
       free(nbin);  /* Chappie not needed anymore */
    }
    else 
    {
       /* 
        * if initialize = 0  ===> I want to skip the VQ and go to EM
        * straight
        */
	if (!read_backup_distribution(initialcodebk,&mean, &variance,&c, 
              &Nmodes, Ndim))
            QUIT(("Unable to read initial distribution\n"));
    }
	
    for (superiter=0;superiter<3;++superiter)
    {
        estimate_multi_modals(vector,numspch,Ndim,Nmodes,mean,variance,c,tempfile,9);
        if (store_distribution(outfile,Nmodes,Ndim,c,mean,variance) != 0)
        {
            printf("Unable to open %s to store distribution\n",outfile);
            printf("Superiter = %d\n",superiter);
        }
    }

    free2d((void**)vector);
    free (c);
    free2d((void**)mean);
    free2d((void**)variance);
    return 0;
}
