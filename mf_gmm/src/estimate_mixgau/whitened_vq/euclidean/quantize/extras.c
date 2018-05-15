#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <sys/stat.h>
#include <sys/file.h>

#include "header.h"

char          **alloc2d(int dim1, int dim2, int size)
{
	size_t   i,nelem;
	char	 *p, **pp;

	nelem = (size_t)dim1 * (size_t)dim2;
	p = calloc(nelem, (size_t) size);
	if (p == NULL) return (NULL);

	pp = (char **) calloc((size_t) dim1, (size_t) sizeof(char *));
	if (pp == NULL) {
		free(p);
		return (NULL);
	}
	for (i = 0; i < (size_t) dim1; i++)
		pp[i] = p + (i * (size_t) dim2 * (size_t) size);

	return (pp);
}

void free2d(void **p)
{
   if (p!=NULL) {
      free(p[0]);
      free(p);
   }
   return;
}

#ifndef FALSE
    #define FALSE   0
    #define TRUE   1
#endif
/* Macro to byteswap an int variable.  x = ptr to variable */
#define SWAP_INT(x)   *(x) = ((0x000000ff & (*(x))>>24) | \
                                (0x0000ff00 & (*(x))>>8) | \
                                (0x00ff0000 & (*(x))<<8) | \
                                (0xff000000 & (*(x))<<24))
/* Macro to byteswap a float variable.  x = ptr to variable */
#define SWAP_FLOAT(x) SWAP_INT((int *) x)

int read_cep(char *file, float***cep, int *numframes, int cepsize)
{
    FILE *fp;
    int n_float;
    struct stat statbuf;
    int i, n, byterev;
    float **mfcbuf;

    if (stat(file, &statbuf) < 0) {
        printf("stat_retry(%s) failed\n", file);
        return -1;
    }

    if ((fp = fopen(file, "rb")) == NULL) {
	printf("fopen(%s,rb) failed\n", file);
	return -1;
    }
    
    /* Read #floats in header */
    if (fread(&n_float, sizeof(int), 1, fp) != 1) {
	fclose (fp);
	return -1;
    }
    
    /* Check if n_float matches file size */
    byterev = FALSE;
    if ((n_float*sizeof(float) + 4) != statbuf.st_size) {
	n = n_float;
	SWAP_INT(&n);

	if ((n*sizeof(float) + 4) != statbuf.st_size) {
	    printf("Header size field: %d(%08x); filesize: %d(%08x)\n",
		    n_float,n_float,(int)statbuf.st_size,(int)statbuf.st_size);
	    fclose (fp);
	    return -1;
	}

	n_float = n;
	byterev = TRUE;
    }
    if (n_float <= 0) {
	printf("Header size field: %d\n",  n_float);
	fclose (fp);
	return -1;
    }
    
    /* n = #frames of input */
    n = n_float/cepsize;
    if (n * cepsize != n_float) {
	printf("Header size field: %d; not multiple of %d\n", n_float, cepsize);
	fclose (fp);
	return -1;
    }

    mfcbuf = (float **) alloc2d (n, cepsize, sizeof(float));
    if (mfcbuf == NULL) {
	printf("Unable to alloc %d floats to read %s\n",n,file);
	return -1;
    }
    
    /* Read mfc data and byteswap if necessary */
    n_float = n * cepsize;
    if (fread (mfcbuf[0], sizeof(float), n_float, fp) != n_float) {
	printf("Error reading mfc data\n");
	fclose (fp);
	return -1;
    }
    if (byterev) {
	for (i = 0; i < n_float; i++)
	    SWAP_FLOAT(&(mfcbuf[0][i]));
    }
    fclose (fp);

    *numframes = n;
    *cep = mfcbuf;
    return (1);
}


/*-------------------------------------- */
int extend_vec(float **vec, int cepdim, int Nvec, int width,
               float ***extendedvec, int *Nfeat)
{
    int i, j, k, l, wdim, nfeat;
    float **evec;

    if (Nvec < width) QUIT(("Insufficient frames to create extended feat\n"))

    *Nfeat = nfeat = Nvec - width  + 1;
    wdim = cepdim*width;

    *extendedvec = evec = (float **)alloc2d(nfeat, wdim, sizeof(float));
    for (i=0; i<nfeat; i++) {
        l = 0;
        for (j=0; j<width; j++) {
            for (k=0; k<cepdim; k++) {
                evec[i][l++] = vec[i+j][k];
            }
        }
    }

    return 1;
}

