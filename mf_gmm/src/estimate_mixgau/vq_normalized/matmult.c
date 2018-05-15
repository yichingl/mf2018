#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_cblas.h>
#include "header.h"

float **matmult(float **a, float **b, int rowsa, int colsa, int colsb)
{
    float **c = (float**) alloc2d(colsb, rowsa, sizeof(float));
    cblas_sgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
                rowsa, colsb, colsa,
                1, a[0], rowsa, b[0], colsa, 0, c[0], rowsa);
    return c;
}


int main(int argc, char **argv)
{
    float **a, **b, **c;
    int i, j, k, rowsa, colsa, colsb;

    rowsa = 3; colsa = 4; colsb = 5;

    a = (float**) alloc2d(colsa, rowsa, sizeof(float));
    b = (float**) alloc2d(colsb, colsa, sizeof(float));

    for (i=0,k=0; i<colsa; i++) for (j=0; j<rowsa; j++) a[i][j] = (float)k++;
    for (i=0; i<colsb; i++) for (j=0; j<colsa; j++) b[i][j] = (float)k++;

    c = matmult(a,b,rowsa,colsa,colsb);

    printf("A :\n");
    for (i=0; i<rowsa; i++) {
        for (j=0; j<colsa; j++) printf("%6.2f ",a[j][i]);
        printf("\n");
    }
    printf("\nB :\n");
    for (i=0; i<colsa; i++) {
        for (j=0; j<colsb; j++) printf("%6.2f ",b[j][i]);
        printf("\n");
    }
    printf("\nC = A*B :\n");
    for (i=0; i<rowsa; i++) {
        for (j=0; j<colsb; j++) printf("%6.2f ",c[j][i]);
        printf("\n");
    }
    return 0;
}
