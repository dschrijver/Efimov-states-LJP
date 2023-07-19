#include <stdlib.h> // rand(), RAND_MAX
#include <unistd.h> // read(), write()

#include "program_functions.h" // x_to_r(), x_to_R2(), p()
#include "matrix_mat.h" // matcopy()
#include "monte_carlo.h"


double monte_carlo_prep(double h, double* raccept, int* domain_bumps, double* resA) {
    int ntotal = nsize*nskip;

    // Initial position in space is at the origin
    double x[3][3] = {0}, r[3][3];
    x_to_r(x, r);
    double w = p(r);

    int naccept = 0;
    *domain_bumps = 0;
    double xnew[3][3];
    double wnew, R2, eta;

    // Burn-in
    for (int i = 0; i < neq; i++) {
        // Propose a new step
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                eta = (double)rand()/(double)RAND_MAX;
                xnew[j][k] = x[j][k] + (double)2.0*h*(eta-(double)0.5);
            }
        }
        // Check if the step does not leave the domain
        R2 = x_to_R2(xnew);
        if (R2 < R2_max) {
            x_to_r(xnew, r);
            wnew = p(r);
            eta = (double)rand()/(double)RAND_MAX;
            // Metropolis sampling condition
            if (wnew > (w*eta)) {
                w = wnew;
                matcopy(3, 3, xnew, x);
            } 
        }
    }

    // Initiate child processes
    int n_samples, mode = 0, division, remainder;
    division = nsize / CHILD_PROCESSES;
    remainder = nsize % CHILD_PROCESSES;
    for (int i = 0; i < CHILD_PROCESSES; i++) {
        if (remainder > 0) {
            n_samples = division + 1;
            remainder--;
        } else {
            n_samples = division;
        }
        write(pipes[2*i][1], &n_samples, sizeof(int));
        write(pipes[2*i][1], &mode, sizeof(int));
        // Send new network parameters
        write(pipes[2*i][1], w1, N_in*N_hid*sizeof(double));
        write(pipes[2*i][1], b1, N_hid*sizeof(double));
        write(pipes[2*i][1], w2, N_hid*sizeof(double));
    }

    // Actual sampling loop
    for (int i = 0; i < ntotal; i++) {
        // Propose a new step
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                eta = (double)rand()/(double)RAND_MAX;
                xnew[j][k] = x[j][k] + (double)2.0*h*(eta-(double)0.5);
            }
        }
        // Check if the step does not leave the domain
        R2 = x_to_R2(xnew);
        if (R2 < R2_max) {
            x_to_r(xnew, r);
            wnew = p(r);
            eta = (double)rand()/(double)RAND_MAX;
            // Metropolis sampling condition
            if (wnew > (w*eta)) {
                w = wnew;
                matcopy(3, 3, xnew, x);
                naccept += 1;
            } 
        } else {
            *domain_bumps += 1;
        }
        if (i % nskip == 0) {
            write(pipes[2*((i/nskip)%CHILD_PROCESSES)][1], x, 9*sizeof(double));
        }
    }

    // Zero the sums and receive the results
    double child_sum, sumA = 0.0, sumA2 = 0.0;
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

    for (int i = 0; i < CHILD_PROCESSES; i++) {
        read(pipes[2*i+1][0], &child_sum, sizeof(double));
        sumA += child_sum;
        read(pipes[2*i+1][0], &child_sum, sizeof(double));
        sumA2 += child_sum;
        read(pipes[2*i+1][0], Oww1, N_in*N_hid*sizeof(double));
        read(pipes[2*i+1][0], Owb1, N_hid*sizeof(double));
        read(pipes[2*i+1][0], Oww2, N_hid*sizeof(double));
        read(pipes[2*i+1][0], AOww1, N_in*N_hid*sizeof(double));
        read(pipes[2*i+1][0], AOwb1, N_hid*sizeof(double));
        read(pipes[2*i+1][0], AOww2, N_hid*sizeof(double));
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

    // Average the sums
    double avg_mul = (double)1.0/(double)nsize;
    *resA = sumA * avg_mul;
    double resA2 = sumA2 * avg_mul;
    for (int i = 0; i < N_hid; i++) {
        sumOwb1[i] *= avg_mul;
        sumOww2[i] *= avg_mul;
        sumAOwb1[i] *= avg_mul;
        sumAOww2[i] *= avg_mul;
        for (int j = 0; j < N_in; j++) {
            sumOww1[j][i] *= avg_mul;
            sumAOww1[j][i] *= avg_mul;
        }
    }

    *raccept = (double)naccept/(double)ntotal;
    return (double)((*resA)*(*resA)/resA2);
}

double monte_carlo_minE(double h, double sigma, double* raccept) {
    int ntotal = nsize*nskip;

    // Initial position in space is at the origin
    double x[3][3] = {0}, r[3][3];
    x_to_r(x, r);
    double w = p(r);

    int naccept = 0;
    double xnew[3][3];
    double wnew, eta;

    // Burn-in
    for (int i = 0; i < neq; i++) {
        // Propose a new step
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                eta = (double)rand()/(double)RAND_MAX;
                xnew[j][k] = x[j][k] + (double)2.0*h*(eta-(double)0.5);
            }
        }
        x_to_r(xnew, r);
        wnew = p(r);
        eta = (double)rand()/(double)RAND_MAX;
        // Metropolis sampling condition
        if (wnew > (w*eta)) {
            w = wnew;
            matcopy(3, 3, xnew, x);
        } 
    }

    // Initiate child processes
    int n_samples, mode = 1, division, remainder;
    division = nsize / CHILD_PROCESSES;
    remainder = nsize % CHILD_PROCESSES;
    for (int i = 0; i < CHILD_PROCESSES; i++) {
        if (remainder > 0) {
            n_samples = division + 1;
            remainder--;
        } else {
            n_samples = division;
        }
        write(pipes[2*i][1], &n_samples, sizeof(int));
        write(pipes[2*i][1], &mode, sizeof(int));
        // Send new network parameters
        write(pipes[2*i][1], w1, N_in*N_hid*sizeof(double));
        write(pipes[2*i][1], b1, N_hid*sizeof(double));
        write(pipes[2*i][1], w2, N_hid*sizeof(double));

        write(pipes[2*i][1], &sigma, sizeof(double));
    }

    // Actual sampling loop
    for (int i = 0; i < ntotal; i++) {
        // Propose a new step
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                eta = (double)rand()/(double)RAND_MAX;
                xnew[j][k] = x[j][k] + (double)2.0*h*(eta-(double)0.5);
            }
        }
        x_to_r(xnew, r);
        wnew = p(r);
        eta = (double)rand()/(double)RAND_MAX;
        // Metropolis sampling condition
        if (wnew > (w*eta)) {
            w = wnew;
            matcopy(3, 3, xnew, x);
            naccept += 1;
        } 
        if (i % nskip == 0) {
            write(pipes[2*((i/nskip)%CHILD_PROCESSES)][1], x, 9*sizeof(double));
        }
    }

    // Zero the sums and receive the results
    double child_sum, sumE = 0.0;
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

    for (int i = 0; i < CHILD_PROCESSES; i++) {
        read(pipes[2*i+1][0], &child_sum, sizeof(double));
        sumE += child_sum;
        read(pipes[2*i+1][0], Oww1, N_in*N_hid*sizeof(double));
        read(pipes[2*i+1][0], Owb1, N_hid*sizeof(double));
        read(pipes[2*i+1][0], Oww2, N_hid*sizeof(double));

        read(pipes[2*i+1][0], OwEw1, N_in*N_hid*sizeof(double));
        read(pipes[2*i+1][0], OwEb1, N_hid*sizeof(double));
        read(pipes[2*i+1][0], OwEw2, N_hid*sizeof(double));
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

    // Average the sums
    double avg_mul = (double)1.0/(double)nsize;
    for (int i = 0; i < N_hid; i++) {
        sumOwb1[i] *= avg_mul;
        sumOww2[i] *= avg_mul;
        sumOwEb1[i] *= avg_mul;
        sumOwEw2[i] *= avg_mul;
        for (int j = 0; j < N_in; j++) {
            sumOww1[j][i] *= avg_mul;
            sumOwEw1[j][i] *= avg_mul;
        }
    }

    *raccept = (double)naccept/(double)ntotal;
    return sumE*avg_mul;
}