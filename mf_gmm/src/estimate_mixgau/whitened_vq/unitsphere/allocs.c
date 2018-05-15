#include <stdio.h>
#include <stdlib.h>

char  **alloc2d(int dim1,  /* "x" dimension */ 
                int dim2,  /* "y" dimension */
                int size   /* number of bytes each entry takes */
               )
{
    size_t        i, nelem; /* total number of elements */
    char           *p,     /* pointer to matrix memory */
              **pp;        /* pointer to matrix mem table */

    /* Total number of elements needed for the two-dimensional matrix */
    nelem = (size_t)dim1 * (size_t)dim2;

    p = calloc(nelem, (size_t) size);
    if (p == NULL) return (NULL);

    /* Now allocate a table for referencing the matrix memory */
    pp = (char **) calloc((unsigned) dim1, (unsigned) sizeof(char *));
    if (pp == NULL) {
        free(p);
        return (NULL);
    }

    /* Fill the table with locations to where each row begins */
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
