#ifndef MAIN_H
#define MAIN_H


#include "constants.h"

// Global on main process, declared and defined in main.c.
// Also global on child processes, but isolated. 
int pipes[CHILD_PROCESSES * 2][2];

double w1[N_in][N_hid];
double b1[N_hid] = {0};
double w2[N_hid];

double u_in[N_in];
double u_hid[N_hid];

double dpsidw1[N_in][N_hid];
double dpsidb1[N_hid];
double dpsidw2[N_hid];

double Oww1[N_in][N_hid];
double Owb1[N_hid];
double Oww2[N_hid];

double sumOww1[N_in][N_hid];
double sumOwb1[N_hid];
double sumOww2[N_hid];

double AOww1[N_in][N_hid];
double AOwb1[N_hid];
double AOww2[N_hid];

double sumAOww1[N_in][N_hid];
double sumAOwb1[N_hid];
double sumAOww2[N_hid];

double dKw1[N_in][N_hid];
double dKb1[N_hid];
double dKw2[N_hid];

double OwEw1[N_in][N_hid];
double OwEb1[N_hid];
double OwEw2[N_hid];

double sumOwEw1[N_in][N_hid];
double sumOwEb1[N_hid];
double sumOwEw2[N_hid];

double dEw1[N_in][N_hid];
double dEb1[N_hid];
double dEw2[N_hid];

double mtw1[N_in][N_hid];
double mtb1[N_hid];
double mtw2[N_hid];

double vtw1[N_in][N_hid];
double vtb1[N_hid];
double vtw2[N_hid];

double C[2][3], D[2][3];

double dpsidr_list[2][3];
double d2psidr2_list[2][3];

int main(int argc, char *argv[]);

#endif