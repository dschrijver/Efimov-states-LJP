#ifndef CALC_GRADIENTS_H
#define CALC_GRADIENTS_H


#include "globals.h"

double K_gradient(double h, double* raccept, int* domain_bumps);
double E_gradient(double h, double sigma, double* raccept);

#endif