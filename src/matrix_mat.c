#include <stdlib.h> // rand(), RAND_MAX
#include <stdio.h> // printf(), fflush(), stdout

#include "matrix_mat.h"


void matuniform(double min, double max, int a, int b, double A[a][b]) {
    double h = max - min;
    for (int i = 0; i < a; i++) {
        for (int j = 0; j < b; j++) {
            A[i][j] = min + h*rand()/RAND_MAX;
        }
    }
    return;
}

void vecuniform(double min, double max, int a, double v[a]) {
    double h = max - min;
    for (int i = 0; i < a; i++) {
        v[i] = min + h*rand()/RAND_MAX;
    }
    return;
}

void matcopy(int a, int b, double A[a][b], double B[a][b]) {
    for (int i = 0; i < a; i++) {
        for (int j = 0; j < b; j++) {
            B[i][j] = A[i][j];
        }
    }
    return;
}

double dot(int a, double u[a], double v[a]) {
    double result = 0;
    for (int i = 0; i < a; i++) {
        result += (double) u[i] * v[i];
    }
    return result;
}

void printmat(int a1, int a2, double A[a1][a2]) {
    printf("{");
    for (int i = 0; i < a1; i++) {
        printf("{");
        for (int j = 0; j < a2; j++) {
            printf("%f", A[i][j]);
            if (j < a2-1) {
                printf(", ");
            }
        }
        printf("}");
        if (i < a1-1) {
            printf(",\n");
        } 
    }
    printf("}\n");
    fflush(stdout);
    return;
}

void printrowvec(int a, double u[a]) {
    printf("{");
    for (int i = 0; i < a; i++) {
        printf("%f", u[i]);
        if (i < a-1) {
            printf(", ");
        } 
    }
    printf("}\n");
    fflush(stdout);
    return;
}
