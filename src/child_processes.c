#include <stdio.h> // printf(), fflush(), stdout
#include <unistd.h> // read(), write(), close()

#include "program_functions.h" // x_to_r(), A(), A2(), E(), Ow(), AOw(), OwE()
#include "child_processes.h"


void child_process(int id) {
    int n_samples, mode;
    double sigma;

    // Close pipes that are not needed
    for (int i = 0; i < CHILD_PROCESSES; i++) {
        if (i == id) {
            close(pipes[2*i][1]);
            close(pipes[2*i+1][0]);
        } else {
            close(pipes[2*i][0]);
            close(pipes[2*i][1]);
            close(pipes[2*i+1][0]);
            close(pipes[2*i+1][1]);
        }
    }

    printf("Started child process %d\n", id);
    fflush(stdout);
    
    while (1) {
        // Wait until Monte Carlo integration has started
        read(pipes[2*id][0], &n_samples, sizeof(int));
        // Preparation: 0, minimizing energy: 1
        read(pipes[2*id][0], &mode, sizeof(int));
        // Update network parameters
        read(pipes[2*id][0], w1, N_in*N_hid*sizeof(double));
        read(pipes[2*id][0], b1, N_hid*sizeof(double));
        read(pipes[2*id][0], w2, N_hid*sizeof(double));

        if (!mode) {
            eval_K(n_samples, id);
        } else {
            read(pipes[2*id][0], &sigma, sizeof(double));
            eval_E(n_samples, id, sigma);
        }
    }
}

void eval_K(int n_samples, int id) {
    double x[3][3], r[3][3];

    // Zero the sums
    double sumA = 0.0, sumA2 = 0.0;
    for (int i = 0; i < N_hid; i++) {
        sumOwb1[i] = 0.0;
        sumOww2[i] = 0.0;
        sumAOwb1[i] = 0.0;
        sumAOww2[i] = 0.0;
        for (int j = 0; j < N_in; j++) {
            sumOww1[j][i] = 0.0;
            sumAOww1[j][i] = 0.0;
        }
    }

    // Receive samples and evaluate the integrands
    for (int i = 0; i < n_samples; i++) {
        read(pipes[2*id][0], &x, 9*sizeof(double));
        x_to_r(x, r);
        sumA += A(r);
        sumA2 += A2(r);
        Ow(r);
        AOw(r);
        for (int i = 0; i < N_hid; i++) {
            sumOwb1[i] += Owb1[i];
            sumOww2[i] += Oww2[i];
            sumAOwb1[i] += AOwb1[i];
            sumAOww2[i] += AOww2[i];
            for (int j = 0; j < N_in; j++) {
                sumOww1[j][i] += Oww1[j][i];
                sumAOww1[j][i] += AOww1[j][i];
            }
        }
    }

    // Send the results back to main process
    write(pipes[2*id+1][1], &sumA, sizeof(double));
    write(pipes[2*id+1][1], &sumA2, sizeof(double));

    write(pipes[2*id+1][1], sumOww1, N_in*N_hid*sizeof(double));
    write(pipes[2*id+1][1], sumOwb1, N_hid*sizeof(double));
    write(pipes[2*id+1][1], sumOww2, N_hid*sizeof(double));

    write(pipes[2*id+1][1], sumAOww1, N_in*N_hid*sizeof(double));
    write(pipes[2*id+1][1], sumAOwb1, N_hid*sizeof(double));
    write(pipes[2*id+1][1], sumAOww2, N_hid*sizeof(double));

    return;
}

void eval_E(int n_samples, int id, double sigma) {
    double x[3][3], r[3][3];
    
    // Zero the sums
    double sumE = 0.0;
    for (int i = 0; i < N_hid; i++) {
        sumOwb1[i] = 0.0;
        sumOww2[i] = 0.0;
        sumOwEb1[i] = 0.0;
        sumOwEw2[i] = 0.0;
        for (int j = 0; j < N_in; j++) {
            sumOww1[j][i] = 0.0;
            sumOwEw1[j][i] = 0.0;
        }
    }

    // Receive samples and evaluate the integrands
    for (int i = 0; i < n_samples; i++) {
        read(pipes[2*id][0], &x, 9*sizeof(double));
        x_to_r(x, r);
        sumE += E(r, sigma);
        Ow(r);
        OwE(r, sigma);
        for (int i = 0; i < N_hid; i++) {
            sumOwb1[i] += Owb1[i];
            sumOww2[i] += Oww2[i];
            sumOwEb1[i] += OwEb1[i];
            sumOwEw2[i] += OwEw2[i];
            for (int j = 0; j < N_in; j++) {
                sumOww1[j][i] += Oww1[j][i];
                sumOwEw1[j][i] += OwEw1[j][i];
            }
        }
    }

    // Send the results back to main process
    write(pipes[2*id+1][1], &sumE, sizeof(double));

    write(pipes[2*id+1][1], sumOww1, N_in*N_hid*sizeof(double));
    write(pipes[2*id+1][1], sumOwb1, N_hid*sizeof(double));
    write(pipes[2*id+1][1], sumOww2, N_hid*sizeof(double));

    write(pipes[2*id+1][1], sumOwEw1, N_in*N_hid*sizeof(double));
    write(pipes[2*id+1][1], sumOwEb1, N_hid*sizeof(double));
    write(pipes[2*id+1][1], sumOwEw2, N_hid*sizeof(double));
    
    return;
}