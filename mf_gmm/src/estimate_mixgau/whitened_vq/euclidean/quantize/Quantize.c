#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "header.h"

/*------------------------------------------------------------------------------
    This function defines the distance metric being used and
    can be varied as per requirement. Metric used here euclidean.
    This particular function also checks if the distance exeeds a
    previously given minimum and returns if it does. Useful to
    speed up clustering..
------------------------------------------------------------------------------*/

float  prune_distance(float *x, float *y, int Ndim, float mindist)
{
    int    i;
    float    dist;
    dist = (x[0]-y[0])*(x[0]-y[0]);
    for (i=1;i<Ndim;++i){
        dist += (x[i]-y[i])*(x[i]-y[i]);
        if (dist > mindist) break;
    }
    return(dist);
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
    for (i=0;i<Nvecs;++i) {
	mindist = 1.0e+32;
        nearest = 0;
        mindist = prune_distance(vector[i],codes[0],Ndim,mindist);
        for (j=1;j<Ncodes;++j){
            codedist = prune_distance(vector[i],codes[j],Ndim,mindist);
            if (codedist < mindist) {
                mindist = codedist;
                nearest = j;
            }
        }
        bins[i] = nearest;
        distortion += mindist;
    }
    return(distortion);
}

/*----------------------------------------------------------------------
  Compute Bag of Words
------------------------------------------------------------------------*/

int  BagofWords(int *bow, float **vector, int Nvecs,
                float **codes, int Ncodes, int Ndim)
{
    int i;
    int *bin = (int *) calloc(Nvecs, sizeof(int));

    Distortion_and_cluster(bin,vector,Nvecs,codes,Ncodes,Ndim);

    for (i=0; i<Ncodes; i++) bow[i] = 0;
    for (i=0; i<Nvecs; i++) bow[bin[i]]++;

    free(bin);

    return 0;
}
