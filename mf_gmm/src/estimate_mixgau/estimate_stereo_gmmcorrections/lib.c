#include <stdio.h>
#include <stdlib.h>

char          **alloc2d(int dim1,	/* "x" dimension */ 
			int dim2,	/* "y" dimension */
			int size	/* number of bytes each entry takes */
			)
{
	int             i;		/* loop control variable */
	unsigned        nelem;		/* total number of elements */

	char	       *p,		/* pointer to matrix memory */
		      **pp;		/* pointer to matrix mem table */

	nelem = (unsigned) dim1 * dim2;

	p = calloc(nelem, (unsigned) size);
	if (p == NULL) return (NULL);

	pp = (char **) calloc((unsigned) dim1, (unsigned) sizeof(char *));
	if (pp == NULL) {
		free(p);
		return (NULL);
	}

	for (i = 0; i < dim1; i++)
		pp[i] = p + (i * dim2 * size);

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

char      ***alloc3d( int       dim1,        /* "x" dimension */
             	      int 	dim2,        /* "y" dimension */
             	      int 	dim3,        /* "z" dimension */
             	      int 	size         /* number of bytes per entry */
		    )
{
    int     j,i;     	/* loop control variables */

    char    *p;
    char    **pp;
    char    ***ppp;        /* pointer to matrix mem table */

    // Allocate for all the elements
    if ((p = (char *)calloc(dim1*dim2*dim3,size)) == NULL)
        return(NULL);

    // Allocate for first level pointers
    if ((pp = (char **)calloc(dim1*dim2,sizeof(char *))) == NULL)
    {
        free(p);
        return(NULL);
    }

    // Allocate for the 3-D level pointers
    if ((ppp = (char ***)calloc(dim1,sizeof(char **))) == NULL)
    {
        free(pp); free(p);
        return(NULL);
    }

    // Fill the table with locations to where each row begins
    for (i = 0, j = 0; i < dim1*dim2; i++, j += dim3)
    {
        pp[i] = p + j*size;
    }
    for (i = 0, j = 0; i < dim1; i++, j += dim2)
    {
        ppp[i] = pp + j;
    }

    return(ppp);
}


void  free3d(void ***X)
{
    // The alloc3d routine simplies the freeing to the following:
    free(X[0][0]);
    free(X[0]);
    free(X);
}

#include <sys/types.h>
#include <sys/stat.h> 

/* Macro to byteswap an int variable.  x = ptr to variable */
#define SWAP_INT(x)   *(x) = ((0x000000ff & (*(x))>>24) | \
                                (0x0000ff00 & (*(x))>>8) | \
                                (0x00ff0000 & (*(x))<<8) | \
                                (0xff000000 & (*(x))<<24))
/* Macro to byteswap a float variable.  x = ptr to variable */
#define SWAP_FLOAT(x) SWAP_INT((int *) x)

#define TRUE 1
#define FALSE 0


int read_cep(char *file, float***cep, int *numframes, int cepsize)
{
    FILE *fp;
    int n_float;
    struct stat statbuf;
    int i, n, byterev, sf, ef;
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
    sf = 0;
    ef = n;

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
    return (0);
}
