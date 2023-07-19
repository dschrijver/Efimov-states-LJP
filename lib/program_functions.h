#ifndef PROGRAM_FUNCTIONS_H
#define PROGRAM_FUNCTIONS_H

#include "globals.h"

double psi_train(double r[3][3]);
double V(double r[3][3], double sigma);
double p(double r[3][3]);
double A(double r[3][3]);
double A2(double r[3][3]);
void Ow(double r[3][3]);
void AOw(double r[3][3]);
double E(double r[3][3], double sigma);
void OwE(double r[3][3], double sigma);
void x_to_r(double x[3][3], double r[3][3]);
void r_to_x(double r[3][3], double x[3][3]);
double x_to_R(double x[3][3]);
double x_to_R2(double x[3][3]);

#endif