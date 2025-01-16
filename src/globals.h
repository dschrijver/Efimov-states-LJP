#ifndef GLOBALS_H
#define GLOBALS_H


#include "../constants.h"

// Global on main process, declared and defined in main.c
extern int pipes[CHILD_PROCESSES * 2][2];

extern double w1[N_in][N_hid];
extern double b1[N_hid];
extern double w2[N_hid];

extern double u_in[N_in];
extern double u_hid[N_hid];

extern double dpsidw1[N_in][N_hid];
extern double dpsidb1[N_hid];
extern double dpsidw2[N_hid];

extern double Oww1[N_in][N_hid];
extern double Owb1[N_hid];
extern double Oww2[N_hid];

extern double sumOww1[N_in][N_hid];
extern double sumOwb1[N_hid];
extern double sumOww2[N_hid];

extern double AOww1[N_in][N_hid];
extern double AOwb1[N_hid];
extern double AOww2[N_hid];

extern double sumAOww1[N_in][N_hid];
extern double sumAOwb1[N_hid];
extern double sumAOww2[N_hid];

extern double dKw1[N_in][N_hid];
extern double dKb1[N_hid];
extern double dKw2[N_hid];

extern double OwEw1[N_in][N_hid];
extern double OwEb1[N_hid];
extern double OwEw2[N_hid];

extern double sumOwEw1[N_in][N_hid];
extern double sumOwEb1[N_hid];
extern double sumOwEw2[N_hid];

extern double dEw1[N_in][N_hid];
extern double dEb1[N_hid];
extern double dEw2[N_hid];

extern double mtw1[N_in][N_hid];
extern double mtb1[N_hid];
extern double mtw2[N_hid];

extern double vtw1[N_in][N_hid];
extern double vtb1[N_hid];
extern double vtw2[N_hid];

extern double C[2][3], D[2][3];

extern double dpsidr_list[2][3];
extern double d2psidr2_list[2][3];

#endif