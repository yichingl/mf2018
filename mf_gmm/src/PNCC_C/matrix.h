/*
 * Header file for matrix operation functions
 */

#define SCALAR double

typedef SCALAR *vector, **matrix;

vector new_vector(int);
matrix new_matrix(int, int);
double innerproduct(int, vector, vector);
double house(int, vector);
void   free_vector(vector);
void   free_matrix(matrix);
void   vecprint(vector, int, int, char*);
void   matprint(matrix, int, int, char*);
void   tridiagonalize(int, matrix, vector, vector);
void   error(char*);

double matinv(int, matrix, matrix);
double gjmatinv(int, matrix);
void   invr(int, matrix);

int eigen(int, matrix, vector, vector);
