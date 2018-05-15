#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "header.h"


/*------------------------------------------------------------------------------
    Function for generating initial codebook
    Will do so by random draws from a Gaussian
------------------------------------------------------------------------------*/
float **initial_codebook(int Ncodes, int Ndim)
{
    float **cbk;
    int i, j;

    cbk = (float **)alloc2d(Ncodes,Ndim,sizeof(float));
    if (cbk == NULL) return cbk;

    for (i=0; i<Ncodes; i++)
        for (j=0; j<Ndim; j++)
            cbk[i][j] = (float)gaussrand();

    for (i=0;i<Ncodes;++i){
        float length = 0;
        for (j=0; j<Ndim; j++) length += cbk[i][j];
        length = sqrt(length);
        for (j=0; j<Ndim; j++) cbk[i][j] /= length;
    }

    return cbk;
}


/*------------------------------------------------------------------------------
    This function defines the distance metric being used.
    We use the projected error:
    ||X - P(X,C)||^2
    where P(X,C) = <X,C>*C
    Minimizing this is identical to maximizing <X,C>
------------------------------------------------------------------------------*/

float  projection(float *x, float *y, int Ndim)
{
    int    i;
    float  t;

    t = x[0]*y[0];
    for (i=1;i<Ndim;++i){
        t += x[i]*y[i];
    }
    return(t);
}

float  projection_error(float *x, float *y, int Ndim, float project)
{
    int    i;
    float  d, err;

    err = 0;
    for (i=0; i<Ndim; ++i){
        d = x[i] - project*y[i];
        err += d*d;
    }
    return(err);
}

/*------------------------------------------------------------------------------
    Function to cluster the training vectors around the codes
    and calculate the sum distortion
------------------------------------------------------------------------------*/

float  project_and_cluster(int *bins, float *bestproj,
                              float **vector, int Nvecs,
                              float **codes, int Ncodes, int Ndim)
{
    int    i,j,nearest;
    float  distortion,codeproj,maxproj;

    distortion = 0;
    for (i=0;i<Nvecs;++i){
        nearest = 0;
        maxproj = projection(vector[i],codes[0],Ndim);
        for (j=1;j<Ncodes;++j){
            codeproj = projection(vector[i],codes[j],Ndim);
            if (codeproj > maxproj){
                maxproj = codeproj;
                nearest = j;
            }
        }
        bins[i] = nearest; bestproj[i] = maxproj;
        distortion += projection_error(vector[i],codes[nearest],Ndim,maxproj);
    }
    return(distortion);
}


/*------------------------------------------------------------------------------
    Function to update codebook according to the clustering information
    receive through the array bin
------------------------------------------------------------------------------*/

/* Function to obtain largest bin in case of an empty bin, to split */

int Largest(int *nfulptr, int *Full, int *Noincode)
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

void  Update_codebook(int *bin, float *proj,
                      float **vector, int Nvecs,
                      float **codes, int Ncodes, int Ndim)
{
    int    i,j,k,l,*Noincode;
    int    nfull,nmpt,*Mpt,*Full;
    int         Largest();

    if ((Noincode = (int *)calloc(Ncodes,sizeof(int))) == NULL)
        QUIT(("Unable to allocate space for array Noincode\n"));

    if ((Mpt = (int *)calloc(Ncodes,sizeof(int))) == NULL)
        QUIT(("Unable to allocate space for array Mpt\n"));

    if ((Full = (int *)calloc(Ncodes,sizeof(int))) == NULL)
        QUIT(("Unable to allocate space for array Full\n"));

    // for (i=0;i<Ncodes;++i) for (j=0;j<Ndim;++j) codes[i][j] = 0;
    
    for (i=0;i<Nvecs;++i){
        j = bin[i];
        for (k=0;k<Ndim;++k)
            codes[j][k] += vector[i][k] * proj[i];
        ++Noincode[j];
    }

    nfull = 0;
    nmpt = 0;
    for (i=0;i<Ncodes;++i){
        if (Noincode[i] > 1){
            Full[nfull++] = i;
            // for (j=0;j<Ndim;++j) codes[i][j] /= (float)Noincode[i];
        }
        else{
            Mpt[nmpt++] = i;
        }
    }
    
    if (nmpt<=nfull){
        for (i=0;i<nmpt;++i){
            j = Mpt[i];
            l = Largest(&nfull,Full,Noincode);
            for (k=0;k<Ndim;++k){
                codes[j][k] = codes[l][k] + 0.50;
                codes[l][k] -= 0.50;
            }
        }
    }
    else QUIT(("CATASTROPHE!\nRerun with better initialization\n"));

    for (i=0;i<Ncodes;++i){
        float length = 0;
        for (j=0; j<Ndim; j++) length += codes[i][j];
        length = sqrt(length);
        for (j=0; j<Ndim; j++) codes[i][j] /= length;
    }


    free(Noincode);
    free(Mpt);
    free(Full);
    return;
}


/*------------------------------------------------------------------
 * Main VQ routine
 * ---------------------------------------------------------------*/

void  vector_quantize(float **codes, int Ncodes, float **vector, int Nvecs,
                      int Ndim, int maxiter)
{
    float    distortion,improvement,prevdist,threshold, *proj;
    int    iter=0, *bin;

    proj = (float *) calloc(Nvecs, sizeof(float));
    bin = (int *) calloc(Nvecs, sizeof(float));

    improvement = 100;
    threshold = 1e-06;

    codes = initial_codebook(Ncodes,Ndim);
    prevdist = project_and_cluster(bin,proj,vector,Nvecs,codes,Ncodes,Ndim);

    printf("Initial Distortion = %f\n",prevdist);
    while (((improvement > threshold)||(improvement < 0))&&(iter<maxiter))
    {
        Update_codebook(bin,proj,vector,Nvecs,codes,Ncodes,Ndim);
        distortion = 
            project_and_cluster(bin,proj,vector,Nvecs,codes,Ncodes,Ndim);
        improvement = (prevdist-distortion)/distortion;
        prevdist = distortion;
        printf("VQ distortion = %f, improvement = %f after %d iterations\n",distortion,improvement,iter); fflush(stdout);
        ++iter;
    }
    printf("VQ exited after %d iterations with total distortion %f\n",iter,prevdist); fflush(stdout);

    free(bin); free(proj);
    return;
}
