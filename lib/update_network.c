#include <stdio.h> // printf(), fflush(), stdout
#include <time.h> // clock(), clock_t, CLOCKS_PER_SEC
#include <math.h> // log(), pow(), sqrt()

#include "calc_gradients.h" // K_gradient(), E_gradient()
#include "update_network.h"


double prep_network() {
    // Zero Adam scheme parameters
    int t = 0;
    for (int i = 0; i < N_hid; i++) {
        mtb1[i] = 0.0;
        vtb1[i] = 0.0;
        mtw2[i] = 0.0;
        vtw2[i] = 0.0;
        for (int j = 0; j < N_in; j++) {
            mtw1[j][i] = 0.0;
            vtw1[j][i] = 0.0;
        }
    }

    double h = h_start;
    double K, raccept;
    int domain_bumps;

    // Main update loop
    for (int i = 0; i < N_prep_updates; i++) {
        clock_t begin = clock();
        // Calculate the gradient
        K = K_gradient(h, &raccept, &domain_bumps);
        // Update the step size
        if (i > start_dynamic_stepsize) {
            h *= log((double) a_ds * P_ideal + b_ds) / log((double) a_ds * raccept + b_ds);
        }
        // Update the network parameters
        t = adam_K(t);
        clock_t end = clock();
        double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
        printf("\rK = %f, domain bumps: %d, step %d/%d, time taken: %.4f         ", K, domain_bumps, i+1, N_prep_updates, time_spent);
        fflush(stdout);
    }

    return h;
}

double minimize_energy(double h) {
    // Zero Adam scheme parameters
    int t = 0;
    for (int i = 0; i < N_hid; i++) {
        mtb1[i] = 0.0;
        vtb1[i] = 0.0;
        mtw2[i] = 0.0;
        vtw2[i] = 0.0;
        for (int j = 0; j < N_in; j++) {
            mtw1[j][i] = 0.0;
            vtw1[j][i] = 0.0;
        }
    }

    double sigma = sigma_min;
    double E, raccept;

    // Main update loop
    for (int i = 0; i < N_minimize_updates; i++) {
        clock_t begin = clock();
        // Calculate the gradient
        E = E_gradient(h, sigma, &raccept);
        // Update the step size
        h *= log((double) a_ds * P_ideal + b_ds) / log((double) a_ds * raccept + b_ds);
        // Update the network parameters
        t = adam_E(t);
        clock_t end = clock();
        double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
        printf("\rE = %f, step %d/%d, time taken: %.4f         ", E, i+1, N_minimize_updates, time_spent);
        fflush(stdout);
    }

    return h;
}

double sweep(double h) {
    // Zero Adam scheme parameters
    int t = 0;
    for (int i = 0; i < N_hid; i++) {
        mtb1[i] = 0.0;
        vtb1[i] = 0.0;
        mtw2[i] = 0.0;
        vtw2[i] = 0.0;
        for (int j = 0; j < N_in; j++) {
            mtw1[j][i] = 0.0;
            vtw1[j][i] = 0.0;
        }
    }

    double sigma = sigma_min;
    double sigma_step = (sigma_end - sigma_min)/((double)N_sweep_updates-(double)1.0);
    double E, raccept;

    // Main update loop
    for (int i = 0; i < N_sweep_updates; i++) {
        clock_t begin = clock();
        // Calculate the gradient
        E = E_gradient(h, sigma, &raccept);
        // Update the step size
        h *= log((double) a_ds * P_ideal + b_ds) / log((double) a_ds * raccept + b_ds);
        // Update the network parameters
        t = adam_E(t);
        clock_t end = clock();
        double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
        printf("\rE = %f, sigma = %.5f, step %d/%d, time taken: %.4f         ", E, sigma, i+1, N_sweep_updates, time_spent);
        fflush(stdout);
        if (h > h_max) {
            break;
        }
        sigma += sigma_step;
    }

    return h;
}

int adam_K(int t) {
    t += 1;

    double mul1 = (double)1.0 / ((double)1.0 - pow(beta1, (double)t));
    double mul2 = (double)1.0 / ((double)1.0 - pow(beta2, (double)t));
    double mthat, vthat;

    for (int i = 0; i < N_hid; i++) {
        mtb1[i] = beta1 * mtb1[i] + ((double)1.0 - beta1) * dKb1[i];
        vtb1[i] = beta2 * vtb1[i] + ((double)1.0 - beta2) * dKb1[i] * dKb1[i];
        mthat = mtb1[i] * mul1;
        vthat = vtb1[i] * mul2;
        b1[i] -= alpha * mthat / (sqrt(vthat) + epsilon);

        mtw2[i] = beta1 * mtw2[i] + ((double)1.0 - beta1) * dKw2[i];
        vtw2[i] = beta2 * vtw2[i] + ((double)1.0 - beta2) * dKw2[i] * dKw2[i];
        mthat = mtw2[i] * mul1;
        vthat = vtw2[i] * mul2;
        w2[i] -= alpha * mthat / (sqrt(vthat) + epsilon);

        for (int j = 0; j < N_in; j++) { 
            mtw1[j][i] = beta1 * mtw1[j][i] + ((double)1.0 - beta1) * dKw1[j][i];
            vtw1[j][i] = beta2 * vtw1[j][i] + ((double)1.0 - beta2) * dKw1[j][i] * dKw1[j][i];
            mthat = mtw1[j][i] * mul1;
            vthat = vtw1[j][i] * mul2;
            w1[j][i] -= alpha * mthat / (sqrt(vthat) + epsilon);
        }
    }

    return t;
}

int adam_E(int t) {
    t += 1;

    double mul1 = (double)1.0 / ((double)1.0 - pow(beta1, (double)t));
    double mul2 = (double)1.0 / ((double)1.0 - pow(beta2, (double)t));
    double mthat, vthat;

    for (int i = 0; i < N_hid; i++) {
        mtb1[i] = beta1 * mtb1[i] + ((double)1.0 - beta1) * dEb1[i];
        vtb1[i] = beta2 * vtb1[i] + ((double)1.0 - beta2) * dEb1[i] * dEb1[i];
        mthat = mtb1[i] * mul1;
        vthat = vtb1[i] * mul2;
        b1[i] -= alpha * mthat / (sqrt(vthat) + epsilon); 

        mtw2[i] = beta1 * mtw2[i] + ((double)1.0 - beta1) * dEw2[i];
        vtw2[i] = beta2 * vtw2[i] + ((double)1.0 - beta2) * dEw2[i] * dEw2[i];
        mthat = mtw2[i] * mul1;
        vthat = vtw2[i] * mul2;
        w2[i] -= alpha * mthat / (sqrt(vthat) + epsilon);

        for (int j = 0; j < N_in; j++) { 
            mtw1[j][i] = beta1 * mtw1[j][i] + ((double)1.0 - beta1) * dEw1[j][i];
            vtw1[j][i] = beta2 * vtw1[j][i] + ((double)1.0 - beta2) * dEw1[j][i] * dEw1[j][i];
            mthat = mtw1[j][i] * mul1;
            vthat = vtw1[j][i] * mul2;
            w1[j][i] -= alpha * mthat / (sqrt(vthat) + epsilon);
        }
    }

    return t;
}