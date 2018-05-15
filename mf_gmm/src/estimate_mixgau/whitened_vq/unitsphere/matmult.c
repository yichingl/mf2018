#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <gsl/gsl_cblas.h>

char **alloc2d(); 
void free2d();

float **matmult(float **a, int rowsa, int colsa,
                float **b, int rowsb, int colsb,
                float **c)
{
    assert(colsa==rowsb);
    if (c == NULL) c = (float**) alloc2d(colsb, rowsa, sizeof(float));
    cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                rowsa, colsb, colsa,
                1, a[0], colsa, b[0], colsb, 0, c[0], colsb);
    return c;
}

float **matmulttxpA(float **a, int rowsa, int colsa,
                    float **b, int rowsb, int colsb,
                    float **c)
{
    assert(rowsa==rowsb);
    if (c==NULL) c = (float **) alloc2d(colsa, colsb, sizeof(float));

    cblas_sgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                colsa, colsb, rowsa,
                1, a[0], colsa, b[0], colsb, 0, c[0], colsb);
    return c;
}

float **matmulttxpB(float **a, int rowsa, int colsa,
                    float **b, int rowsb, int colsb,
                    float **c)
{
    assert(colsa==colsb);
    if (c==NULL) c = (float **) alloc2d(rowsa, rowsb, sizeof(float));

    cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
                rowsa, rowsb, colsa,
                1, a[0], colsa, b[0], colsb, 0, c[0], rowsb);
    return c;
}

int main(int argc, char **argv)
{
    float **a, **b, **c, *d, *e, **f;
    int i, j, k, rowsa, rowsb, colsa, colsb;

    rowsa = 3; colsa = rowsb = 4; colsb = 5;

    a = (float**) alloc2d(rowsa, colsa, sizeof(float));
    f = (float**) alloc2d(colsa, rowsa, sizeof(float));
    b = (float**) alloc2d(rowsb, colsb, sizeof(float));
    c = (float**) alloc2d(rowsa, colsb, sizeof(float));
    d = (float*) calloc(colsa, sizeof(float));
    e = (float*) calloc(colsb, sizeof(float));

    for (i=0,k=0; i<rowsa; i++) for (j=0; j<colsa; j++) {f[j][i] = a[i][j] = (float)k++;};
    for (i=0; i<rowsb; i++) for (j=0; j<colsb; j++) b[i][j] = (float)k++;
    for (j=0; j<colsa; j++) d[j] = 2*(float)j + 1;

    matmult(a,rowsa,colsa,b,rowsb,colsb,c);

    printf("A :\n");
    for (i=0; i<rowsa; i++) {
        for (j=0; j<colsa; j++) printf("%6.2f ",a[i][j]);
        printf("\n");
    }
    printf("F = A' :\n");
    for (i=0; i<colsa; i++) {
        for (j=0; j<rowsa; j++) printf("%6.2f ",f[i][j]);
        printf("\n");
    }
    printf("\nB :\n");
    for (i=0; i<rowsb; i++) {
        for (j=0; j<colsb; j++) printf("%6.2f ",b[i][j]);
        printf("\n");
    }
    printf("\nC = A*B :\n");
    for (i=0; i<rowsa; i++) {
        for (j=0; j<colsb; j++) printf("%6.2f ",c[i][j]);
        printf("\n");
    }

    free2d((void**)c);
    c = matmulttxpA(f,colsa,rowsa,b,rowsb,colsb,NULL);
    printf("\n\nC = F'*B :\n");
    for (i=0; i<rowsa; i++) {
        for (j=0; j<colsb; j++) printf("%6.2f ",c[i][j]);
        printf("\n");
    }
    printf("\n");

    c = matmulttxpA(a,rowsa,colsa,a,rowsa,colsa,NULL);
    printf("\n\nC = A'*A :\n");
    for (i=0; i<colsa; i++) {
        for (j=0; j<colsa; j++) printf("%6.2f ",c[i][j]);
        printf("\n");
    }

    free2d((void**)c);
    c = matmulttxpB(b,rowsb,colsb,b,rowsb,colsb,NULL); /* colsa = rowsb */
    printf("\n\nC = B*B' :\n");
    for (i=0; i<rowsb; i++) {
        for (j=0; j<rowsb; j++) printf("%6.2f ",c[i][j]);
        printf("\n");
    }

    printf("\nD :\n");
    for (i=0; i<colsa; i++) printf("%6.2f ",d[i]); printf("\n");
    matmult(&d,1,colsa,b,rowsb,colsb,&e);
    printf("\nD*B :\n");
    for (i=0; i<colsb; i++) printf("%6.2f ",e[i]); printf("\n");

    return 0;
}
