#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "header.h"


/*------------------------------------------------------------------------------
    This function defines the distance metric being used and
    can be varied as per requirement. Metric used here euclidean.
    This particular function also checks if the distance exeeds a
    previously given minimum and returns if it does. Useful to
    speed up clustering..
------------------------------------------------------------------------------*/

float  prune_distance(float *x,float *y,int Ndim,float mindist)
{
    int    i;
    float    dist,t;

    t = *x++ - *y++;
    dist = t * t;
    for (i=1;i<Ndim;++i){
        t = *x++ - *y++;
        dist += t * t;
        if (dist > mindist) break;
    }
    return(dist);
}

/*------------------------------------------------------------------------------
    Function for generating initial codebook
------------------------------------------------------------------------------*/
void Initial_codebook(int Ncodes, float **codes, 
                      int Nvecs, float **vector, int Ndim, int *bin)
{
    int    i,j,nbin,idx;
    float  mindist,codedist;

    /*
    for (i=0;i<Ncodes;++i){
        for (j=0;j<Ndim;++j)
            codes[i][j] = vector[i][j];
        bin[i] = 1;
    }
    */

    /* Just pick Ncodes vectors randomly from input */
    int *ndone = (int *)calloc(Ncodes,sizeof(int));
    for (i=0;i<Ncodes;++i){
        int repeat = 1;
        idx = (int)(drand48()*(float)Nvecs);
        while (repeat) { /* Ensure we arent reusing vectors */
            repeat = 0;
            for (j=0; j<i; j++) if (ndone[j] == idx) repeat = 1;
            if (repeat) idx = (int)(drand48()*(float)Nvecs);
        }
        for (j=0;j<Ndim;++j)
            codes[i][j] = vector[idx][j];
        ndone[i] = idx;
    }
    free(ndone);

    for (i=0;i<Nvecs;++i){
        nbin = 0; mindist = 1.0e+32;
        mindist = prune_distance(vector[i],codes[0],Ndim,mindist);
        for (j=1;j<Ncodes;++j){
            codedist = prune_distance(vector[i],codes[j],Ndim,mindist);
            if (codedist < mindist){
                mindist = codedist;
                nbin = j;
            }
        }
        for (j=0;j<Ndim;++j){
            codes[nbin][j] = codes[nbin][j]*bin[nbin]+vector[i][j];
            codes[nbin][j] /= (bin[nbin]+1);
        }
        ++bin[nbin];
    }
    return;
}


/*------------------------------------------------------------------------------
    Function to cluster the training vectors around the codes
    and calculate the sum distortion
------------------------------------------------------------------------------*/

float  Distortion_and_cluster(int *bins, float **vector, int Nvecs,
                              float **codes, int Ncodes, int Ndim)
{
    int    i,j,nearest;
    float    distortion,codedist,mindist;

    distortion = 0;
    for (i=0;i<Nvecs;++i){
	mindist = 1.0e+32;
        nearest = 0;
        mindist = prune_distance(vector[i],codes[0],Ndim,mindist);
        for (j=1;j<Ncodes;++j){
            codedist = prune_distance(vector[i],codes[j],Ndim,mindist);
            if (codedist < mindist){
                mindist = codedist;
                nearest = j;
            }
        }
        bins[i] = nearest;
        distortion += mindist;
    }
    return(distortion);
}


/*------------------------------------------------------------------------------
    Function to update codebook according to the clustering information
    receive through the array bin
------------------------------------------------------------------------------*/

/* Function to obtain largest bin in case of an empty bin, to split */

int Largest(int *nfulptr,int *Full,int *Noincode)
{
    int i,j,ibig,ifull;

    ifull = 0; ibig = Full[ifull];
    for (i=1;i<*nfulptr;++i){
        j = Full[i];
        if (Noincode[j]>Noincode[ibig]){
            ibig = j;
            ifull = i;
        }
    }

    for (i=ifull;i<*nfulptr-1;++i)
        Full[i] = Full[i+1];
    --(*nfulptr);
    
    return(ibig);
}

void  Update_codebook(int *bin, float **vector, int Nvecs,
                      float **codes, int Ncodes, int Ndim)
{
    int    i,j,k,l,*Noincode;
    int    nfull,nmpt,*Mpt,*Full;

    if ((Noincode = (int *)malloc(Ncodes*sizeof(int))) == NULL)
        QUIT(("Unable to allocate space for array Noincode\n"));

    if ((Mpt = (int *)malloc(Ncodes*sizeof(int))) == NULL)
        QUIT(("Unable to allocate space for array Mpt\n"));

    if ((Full = (int *)malloc(Ncodes*sizeof(int))) == NULL)
        QUIT(("Unable to allocate space for array Full\n"));

    for (i=0;i<Ncodes;++i){
        Noincode[i] = 0;
        for (j=0;j<Ndim;++j)
            codes[i][j] = 0;
    }
    
    for (i=0;i<Nvecs;++i){
        j = bin[i];
        for (k=0;k<Ndim;++k)
            codes[j][k] += vector[i][k];
        ++Noincode[j];
    }

    nfull = 0;
    nmpt = 0;
    for (i=0;i<Ncodes;++i){
        if (Noincode[i] > 1){
            Full[nfull] = i;
            ++nfull;
            for (j=0;j<Ndim;++j)
                codes[i][j] /= (float)Noincode[i];
        }
        else{
            Mpt[nmpt] = i;
            ++nmpt;
        }
    }
    
    if (nmpt<=nfull){
        for (i=0;i<nmpt;++i){
            j = Mpt[i];
            l = Largest(&nfull,Full,Noincode);
            for (k=0;k<Ndim;++k){
                codes[j][k] = codes[l][k] * 0.90;
                codes[l][k] *= 1.10;
            }
        }
    }
    else
        QUIT(("CATASTROPHE, CATASTROPHE\nRerun with better initialization\n"));


    free(Noincode);
    free(Mpt);
    free(Full);
    return;
}

/*--------------------------------------------------------------------
 * The actual VQ routine
 *--------------------------------------------------------------------*/

void  vector_quantize(float **codes, int Ncodes, 
                      float **vector, int Nvecs, int Ndim, int *bin)    
{
    float    distortion,improvement,prevdist,threshold;
    int    iter=0;

    improvement = 100;
    threshold = 1e-06;


    Initial_codebook(Ncodes,codes,Nvecs,vector,Ndim,bin);
    prevdist = Distortion_and_cluster(bin,vector,Nvecs,codes,Ncodes,Ndim);

    printf("Initial Distortion = %f\n",prevdist);
    while (((improvement > threshold)||(improvement < 0))&&(iter<5))
    {
        Update_codebook(bin,vector,Nvecs,codes,Ncodes,Ndim);
        distortion = 
            Distortion_and_cluster(bin,vector,Nvecs,codes,Ncodes,Ndim);
        improvement = (prevdist-distortion)/distortion;
        prevdist = distortion;
        printf("VQ distortion = %f, improvement = %f after %d iterations\n",distortion,improvement,iter); fflush(stdout);
        ++iter;
    }
    printf("VQ exited after %d iterations with total distortion %f\n",iter,prevdist); fflush(stdout);

    return;
}



void initialize_mgau(float **vector, int nvec, int Ndim,
                     float **mean, float **var, float *c, int Nmodes)
{
    int i, j, k, *nbin, *bin;
    nbin = (int *)calloc(Nmodes,sizeof(int)); /* no of vectors in a mode */
    if (nbin==NULL)
        QUIT(("Unable to allocate %d space for nbin array\n",Nmodes));

   /*
   * The vector_quantize routine performs VQ with a mahalonobis metric
   * and returns the codes as the means and the wieghts as the variances
   * of the initial estimates of the modes, which will further be 
   * employed in EM. Note that the variances are purely diagonal
   * We initialize all initial c[] to be equal
   */

   if ((bin = (int *)calloc(nvec,sizeof(int))) == NULL)
       QUIT(("Unable to allocate %d floats for bin array\n",nvec));

   vector_quantize(mean,Nmodes,vector,nvec,Ndim,bin);

   for (i=0;i<Nmodes;++i)
       c[i] = 1.0/(float)Nmodes;
   for (k=0;k<Nmodes;++k) {
       nbin[k] = 0;
       for (i=0;i<Ndim;++i)
           var[k][i] = 0;
   }
   for (i=0;i<nvec;++i) {
       for (j=0;j<Ndim;++j) {
           float atemp = (vector[i][j] - mean[bin[i]][j]);
           var[bin[i]][j] += atemp*atemp;
       }
       ++nbin[bin[i]];
   }
   for (k=0;k<Nmodes;++k) {
       for(j=0;j<Ndim;++j)
           var[k][j] /= nbin[k];
   }

   free(bin); 
   free(nbin);
}
