#ifndef CHILD_PROCESSES_H
#define CHILD_PROCESSES_H


#include "globals.h"

void child_process(int id);
void eval_K(int n_samples, int id);
void eval_E(int n_samples, int id, double sigma);

#endif