#ifndef MATRIX_MAT_H
#define MATRIX_MAT_H


void matuniform(double min, double max, int a, int b, double A[a][b]);
void vecuniform(double min, double max, int a, double v[a]);
void matcopy(int a, int b, double A[a][b], double B[a][b]);
double dot(int a, double u[a], double v[a]);
void printmat(int a1, int a2, double A[a1][a2]);
void printrowvec(int a, double u[a]);

#endif