/*-----------------------------------------------------------------------------
    Function for Vector Quantization and generation of codebook,
    and binning information.
    The codebook generated is useful mainly for initialization 
    of other algorthims as the LBG algorithm is stopped after 7 iterations.
    The initialization of the VQ is done by picking the first few vectors
    and binning around them using a K-means.

    Parameters Passed
      Input:
	vector[][] - Input data set, (float)
	Nvecs - Number of vectors in data set (int)
	Ndim - Dimension of input set (int)
	Ncodes - Number of codes needed in codebook(int)

      Out:
	codes[][] - The codebook (float)
	bin[] - Binning information - bin[i] tells us which code the vector
		vector[i] belongs to.

     Caveat: This routine can only be used if vector[][] and codes[][]
             have been allocated using alloc2d()
-----------------------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "header.h"

void  vector_quantize(codes,Ncodes,vector,Nvecs,Ndim,bin)    
float    **vector,**codes;
int    Ncodes,Nvecs,Ndim,*bin;
{
    float    distortion,improvement,prevdist,threshold;
    int    iter=0;
    extern void    Initial_codebook(),Update_codebook();
    extern float    Distortion_and_cluster();

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


/*------------------------------------------------------------------------------
    Function for generating initial codebook
------------------------------------------------------------------------------*/
void Initial_codebook(Ncodes,codes,Nvecs,vector,Ndim,bin)
    float    **codes,**vector;
    int      Ncodes,Nvecs,Ndim,*bin;
    {
    auto int    i,j,nbin;
    auto float  mindist,codedist;
    float       prune_distance();

    for (i=0;i<Ncodes;++i){
        for (j=0;j<Ndim;++j)
            codes[i][j] = vector[i][j];
        bin[i] = 1;
    }

    for (i=Ncodes;i<Nvecs;++i){
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

float  Distortion_and_cluster(bins,vector,Nvecs,codes,Ncodes,Ndim)
    float    **vector,**codes;
    int    Nvecs,Ncodes,*bins,Ndim;
    {
    auto int    i,j,nearest;
    auto float    distortion,codedist,mindist;
    float        prune_distance();

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

void  Update_codebook(bin,vector,Nvecs,codes,Ncodes,Ndim)
    float    **vector,**codes;
    int    Nvecs,Ncodes,Ndim,*bin;
    {
    auto int    i,j,k,l,*Noincode;
    auto int    nfull,nmpt,*Mpt,*Full;
    int         Largest();

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

/* Function to obtain largest bin in case of an empty bin, to split */

int Largest(nfulptr,Full,Noincode)
    int    *Full,*Noincode,*nfulptr;
    {
    auto int i,j,ibig,ifull;

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


/*------------------------------------------------------------------------------
    This function defines the distance metric being used and
    can be varied as per requirement. Metric used here euclidean.
    This particular function also checks if the distance exeeds a
    previously given minimum and returns if it does. Useful to
    speed up clustering..
------------------------------------------------------------------------------*/

float  prune_distance(x,y,Ndim,mindist)
    float    *x,*y,mindist;
    int    Ndim;
    {
    auto int    i;
    auto float    dist,t;

    t = *x++ - *y++;
    dist = t * t;
    for (i=1;i<Ndim;++i){
        t = *x++ - *y++;
        dist += t * t;
        if (dist > mindist) break;
    }
    return(dist);
}
