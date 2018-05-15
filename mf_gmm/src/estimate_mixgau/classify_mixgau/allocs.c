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
