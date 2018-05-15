#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "header.h"

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


/*------------------------------------------------------------------
 * Main BOW routine
 * ---------------------------------------------------------------*/

void bagofwords(float **codes, int Ncodes, float **vector,
                 int Nvecs, int Ndim, int *counts, float *projcounts)
{
    float  distortion, *projs;
    int    i, j, *bin;

    bin = (int *) calloc(Nvecs, sizeof(int));
    projs = (float *) calloc(Nvecs, sizeof(float));

    for (i=0; i<Ncodes; i++) {counts[i] = 0; projcounts[i] = 0;}
    distortion = project_and_cluster(bin,projs,vector,Nvecs,codes,Ncodes,Ndim);
    for (i=0; i<Nvecs; i++) {
        j = bin[i];
        counts[j]++;
        projcounts[j] += projs[i];
    }
    printf("Bag of words computed; total error = %f\n",distortion);

    free(bin); free(projs);
    return;
}
