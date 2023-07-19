#ifndef NETWORK_DERIVATIVES_H
#define NETWORK_DERIVATIVES_H


#include "globals.h"

void dpsi_dw(double r[3][3]);
void dpsi_dr(double r[3][3]);
void d2psi_dr2(double r[3][3]);

#endif