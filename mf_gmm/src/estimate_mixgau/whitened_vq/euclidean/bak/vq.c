#ifdef MULTITHREAD
    #define _GNU_SOURCE
    #include <sched.h>
    #include <pthread.h>
#endif

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#ifdef SIMD
    #include <x86intrin.h>
#endif

#include "header.h"

#define MAXITER 15

#ifdef MULTITHREAD
typedef struct t_arg {
    int id;
    int count;
    int *bin;
    float **vector;
    float **codes;
    int Ncodes;
    int Ndim;
    float *distortion;
} t_args;
#endif

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

#ifdef SIMD
    int j;
    __m128 X, Y, Z, W;

    dist = 0;
    W = _mm_setr_ps(0.0f, 0.0f, 0.0f, 0.0f);
    for (i=0; i<Ndim-4; i+=4) {
        X = _mm_load_ps(&x[i]);
        Y = _mm_load_ps(&y[i]);
        Z = _mm_sub_ps(X, Y);
        X = _mm_mul_ps(Z, Z);
        W = _mm_add_ps(X, W);
    }
    for (j=0; j<4; j++) dist += W[j];

    X = _mm_load_ps(&x[i]);
    Y = _mm_load_ps(&y[i]);
    Z = _mm_sub_ps(X, Y);
    X = _mm_mul_ps(Z, Z);
    for (j=0;j<Ndim-i;++j){
        dist += X[j];
    }
#else
    dist = (x[0]-y[0])*(x[0]-y[0]);
    for (i=1;i<Ndim;++i){
        dist += (x[i]-y[i])*(x[i]-y[i]);
        if (dist > mindist) break;
    }
#endif
    return(dist);
}



/*------------------------------------------------------------------------------
    Function for generating initial codebook
------------------------------------------------------------------------------*/
void Initial_codebook(int Ncodes, float **codes, int Nvecs,
                      float **vector, int Ndim, int *bin)
{
    int    i,j,nbin;
    float  mindist,codedist;

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

#ifdef MULTITHREAD
void  *Distortion_and_cluster(void *threadarg)
{
    int *bins, Ncodes, Ndim, Nvecs, thrid;
    float **vector, **codes;
    t_args  *args = (t_args *) threadarg;

    thrid = args->id;
    vector = args->vector;
    codes = args->codes;
    bins = args->bin;
    Nvecs = args->count;
    Ncodes = args->Ncodes;
    Ndim = args->Ndim;
#else
float  Distortion_and_cluster(int *bins, float **vector, int Nvecs,
                              float **codes, int Ncodes, int Ndim)
{
#endif
    int    i, j, nearest;
    float   distortion, codedist, mindist;

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
#ifdef MULTITHREAD
    (args->distortion)[0] = distortion;
    pthread_exit((void *)threadarg);
#else
    return(distortion);
#endif
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

#ifdef MULTITHREAD
    #define NUM_THREADS 8
#endif

void  vector_quantize(float **codes, int Ncodes, float **vector,
                      int Nvecs, int Ndim, int maxiter)    
{
    float  totaldistortion, improvement, prevdist, threshold;
    int    iter=0, *bin;
    clock_t start, end;
    double cpu_time_used;

#ifdef MULTITHREAD
    int    ithr, startvec, count, rc;
    float distortion[NUM_THREADS];
    pthread_t thread[NUM_THREADS];
    t_args *thrargs;
    pthread_attr_t attr;
    void *status;

    thrargs = (t_args*) calloc(NUM_THREADS, sizeof(t_args));
#endif

    bin = (int *) calloc(Nvecs, sizeof(int));
    if (!bin) QUIT(("Unable to alloc %d ints for bin array\n",Nvecs))
    if (maxiter <= 0) maxiter = MAXITER;

    improvement = 100;
    threshold = 1e-06;

    Initial_codebook(Ncodes,codes,Nvecs,vector,Ndim,bin);
    start = clock(); prevdist = 1e+30;
#ifdef MULTITHREAD
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
#endif
    while (((improvement > threshold)||(improvement < 0))&&(iter<maxiter)) {
        char *cbkfile = "Codebook.tmp";
#ifdef MULTITHREAD
        for (ithr=0; ithr<NUM_THREADS; ithr++) {
            startvec = (int)((float)Nvecs/(float)NUM_THREADS * ithr);
            count = (int)((float)Nvecs/(float)NUM_THREADS * (ithr+1)) - startvec;

            thrargs[ithr].id = ithr;
            thrargs[ithr].count = count;
            thrargs[ithr].bin = bin+startvec;
            thrargs[ithr].vector = vector+startvec;
            thrargs[ithr].codes = codes;
            thrargs[ithr].Ncodes = Ncodes;
            thrargs[ithr].Ndim = Ndim;
            thrargs[ithr].distortion = &distortion[ithr];

            rc = pthread_create(&thread[ithr],NULL,Distortion_and_cluster,(void*)&thrargs[ithr]);
            if (rc) QUIT(("ERROR: return code from pthread_create() is %d\n",rc))
        }
        for (totaldistortion = 0, ithr=0; ithr<NUM_THREADS; ithr++) {
            rc = pthread_join(thread[ithr], &status);
            if (rc) QUIT(("ERROR: return code from pthread_join() is %d\n",rc))
            totaldistortion += distortion[ithr];
        }
#else
        totaldistortion =
            Distortion_and_cluster(bin,vector,Nvecs,codes,Ncodes,Ndim);
#endif
        Update_codebook(bin,vector,Nvecs,codes,Ncodes,Ndim);
        improvement = (prevdist-totaldistortion)/totaldistortion;
        prevdist = totaldistortion;
        end = clock();
        cpu_time_used = ((double)(end-start)) / CLOCKS_PER_SEC;
        printf("VQ distortion = %f, improvement = %f after %d iterations in %f secs\n",totaldistortion,improvement,iter,cpu_time_used); fflush(stdout);
        savefile(cbkfile, codes, Ndim, -1, Ncodes);
        ++iter;
    }
    printf("VQ exited after %d iterations with total distortion %f\n",iter,prevdist); fflush(stdout);
    free(bin);
#ifdef MULTITHREAD
    pthread_attr_destroy(&attr);
    // pthread_exit(NULL);
#endif

    return;
}
