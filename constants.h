#ifndef CONSTANTS_H
#define CONSTANTS_H


// I usually choose as many processes as my CPU has cores, but
// that does not mean that that is the optimal choice. 
#define TOTAL_PROCESSES         (int) 6
#define CHILD_PROCESSES         (int) (TOTAL_PROCESSES - 1)

// Amountd of input neurons and hidden neurons
#define N_in                    (int) 3 // Don't change this one!
#define N_hid                   (int) 30 // This one can be changed freely

// General settings
#define N_prep_updates          (int) 3000
#define N_minimize_updates      (int) 5000
#define N_sweep_updates         (int) 20000
#define start_dynamic_stepsize  (int) 100
#define R_max                   (double) 15.0
#define h_start                 (double) 3.0
#define sigma_min               (double) 0.628094308 // a = +0.5
#define sigma_end               (double) 0.89663     // a = -100
// #define sigma_end               (double) 0.920174361 // a = -8.0
#define h_max                   (double) 100.0 
#define R2_max                  (double) (R_max*R_max)

// Monte Carlo settings
#define neq                     (int) 10000
#define nskip                   (int) 12
#define nsize                   (int) 10000

// Cut-off of potential at short range
#define V_max                   (double) 100.0

// Dynamic stepsize
// The constants a and b need to be recalculated when a
// different ideal acceptance rate is chosen. 
#define a_ds                    (double) 0.548446959 
#define b_ds                    (double) 0.004860716
#define P_ideal                 (double) 0.32 

// Adam scheme settings
#define alpha                   (double) 0.001
#define beta1                   (double) 0.9
#define beta2                   (double) 0.999
#define epsilon                 (double) 10e-8

#endif