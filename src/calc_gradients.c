#include "monte_carlo.h"
#include "calc_gradients.h"


double K_gradient(double h, double* raccept, int* domain_bumps) {
    double resA;
    double K = monte_carlo_prep(h, raccept, domain_bumps, &resA);

    for (int i = 0; i < N_hid; i++) {
        dKb1[i] = (double)-2.0*K*(sumAOwb1[i] / resA - sumOwb1[i]);
        dKw2[i] = (double)-2.0*K*(sumAOww2[i] / resA - sumOww2[i]);
        for (int j = 0; j < N_in; j++) {
            dKw1[j][i] = (double)-2.0*K*(sumAOww1[j][i] / resA - sumOww1[j][i]);
        }
    }

    return K;
}

double E_gradient(double h, double sigma, double* raccept) {
    double E_res = monte_carlo_minE(h, sigma, raccept);

    for (int i = 0; i < N_hid; i++) {
        dEb1[i] = (double)2.0*(sumOwEb1[i] - E_res*sumOwb1[i]);
        dEw2[i] = (double)2.0*(sumOwEw2[i] - E_res*sumOww2[i]);
        for (int j = 0; j < N_in; j++) {
            dEw1[j][i] = (double)2.0*(sumOwEw1[j][i] - E_res*sumOww1[j][i]);
        }
    }

    return E_res;
}