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

	/*
	 * Compute total number of elements needed for the two-dimensional
	 * matrix
	 */
	nelem = (unsigned) dim1 * dim2;

	/*
	 * Allocate the memory needed for the matrix
	 */
	p = calloc(nelem, (unsigned) size);

	/*
	 * If the allocation were not successful, return a NULL pointer
	 */
	if (p == NULL) return (NULL);

	/*
	 * Now allocate a table for referencing the matrix memory
	 */
	pp = (char **) calloc((unsigned) dim1, (unsigned) sizeof(char *));

	/*
	 * If the allocation were not successful, return a NULL pointer
	 * and free the previously allocated memory
	 */
	if (pp == NULL) {
		free(p);
		return (NULL);
	}

	/*
	 * Fill the table with locations to where each row begins
	 */
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


void  free3d(char ***X)
{
    // The alloc3d routine simplies the freeing to the following:
    free(X[0][0]);
    free(X[0]);
    free(X);
}
