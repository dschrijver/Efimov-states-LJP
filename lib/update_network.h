#ifndef UPDATE_NETWORK_H
#define UPDATE_NETWORK_H


#include "globals.h"

double prep_network();
double minimize_energy(double h);
double sweep(double h);
int adam_K(int t);
int adam_E(int t);

#endif