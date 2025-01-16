#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H


#include "globals.h"

double monte_carlo_prep(double h, double* raccept, int* domain_bumps, double* resA);
double monte_carlo_minE(double h, double sigma, double* raccept);

#endif 