#include <math.h> // exp(), erf(), sqrt(), pow()

#include "network_calc.h" // calc()
#include "network_derivatives.h" // dpsi_dw(), d2psi_dr2()
#include "matrix_mat.h" // dot()
#include "program_functions.h"


double psi_train(double r[3][3]) {
    double x[3][3];
    r_to_x(r, x);
    double R = x_to_R(x);
    double res1 = exp(-(R-(double)1.5)*(R-(double)1.5)/((double)1.1));
    double res2 = (double)1.0 + erf((double)2.12*(R - (double)1.5));
    return res1*res2;
}

double V(double r[3][3], double sigma) {
    double sigma4 = pow(sigma, (double)4.0);
    double sum = 0.0, dist1 = 0.0, dist2 = 0.0, dist3 = 0.0;
    double res, mul = sqrt((double)3.0);
    for (int i = 0; i < 3; i++) {
        dist1 += r[2][i]*r[2][i];
        res = (double)0.5*(mul*r[1][i]+r[2][i]);
        dist2 += res*res;
        res = (double)0.5*(mul*r[1][i]-r[2][i]);
        dist3 += res*res;
    }
    sum += sigma4/(pow(dist1, (double)5.0)) - (double)1.0/(dist1*dist1*dist1);
    sum += sigma4/(pow(dist2, (double)5.0)) - (double)1.0/(dist2*dist2*dist2);
    sum += sigma4/(pow(dist3, (double)5.0)) - (double)1.0/(dist3*dist3*dist3);
    sum *= (double)16.0;
    if (sum > V_max) {
        sum = V_max;
    } 
    return sum;
}

double p(double r[3][3]) {
    double result = calc(r);
    return result*result;
}

double A(double r[3][3]) {
    return psi_train(r) / calc(r);
}

double A2(double r[3][3]) {
    double result = A(r);
    return result*result;
}

void Ow(double r[3][3]) {
    dpsi_dw(r);
    double inv_psi_val = (double)1.0/calc(r);
    for (int i = 0; i < N_hid; i++) {
        Owb1[i] = dpsidb1[i]*inv_psi_val;
        Oww2[i] = dpsidw2[i]*inv_psi_val;
        for (int j = 0; j < N_in; j++) {
            Oww1[j][i] = dpsidw1[j][i]*inv_psi_val;
        }
    }
    return;
}

void AOw(double r[3][3]) {
    dpsi_dw(r);
    double inv_psi_val = (double)1.0/calc(r);
    double A_res = A(r);
    for (int i = 0; i < N_hid; i++) {
        AOwb1[i] = dpsidb1[i]*inv_psi_val*A_res;    
        AOww2[i] = dpsidw2[i]*inv_psi_val*A_res;
        for (int j = 0; j < N_in; j++) {
            AOww1[j][i] = dpsidw1[j][i]*inv_psi_val*A_res;
        }
    }
    return;
}

double E(double r[3][3], double sigma) {
    d2psi_dr2(r);
    double psi = calc(r);
    double H_kin = 0.0, V_int;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 3; j++) {
            H_kin -= d2psidr2_list[i][j];
        }
    }
    H_kin /= psi;
    V_int = V(r, sigma);
    return H_kin + V_int;
}

void OwE(double r[3][3], double sigma) {
    Ow(r);
    double E_res = E(r, sigma);
    for (int i = 0; i < N_hid; i++) {
        OwEb1[i] = E_res * Owb1[i];
        OwEw2[i] = E_res * Oww2[i];
        for (int j = 0; j < N_in; j++) {
            OwEw1[j][i] = E_res * Oww1[j][i];
        }
    }
    return;
}

void x_to_r(double x[3][3], double r[3][3]) {
    for (int j = 0; j < 3; j++) {
        r[0][j] = -sqrt((double)1.0/(double)6.0)*(x[0][j] + x[1][j] + x[2][j]);
        r[1][j] = sqrt((double)4.0/(double)3.0)*(x[2][j] - (double)0.5*(x[0][j] + x[1][j]));
        r[2][j] = x[1][j]-x[0][j];
    }
    return;
}

void r_to_x(double r[3][3], double x[3][3]) {
    double mul1 = sqrt((double)6.0), mul2 = sqrt((double)3.0);
    for (int j = 0; j < 3; j++) {
        x[0][j] = (double)1.0/(double)6.0*((double)-2.0*mul1*r[0][j] - mul2*r[1][j] - (double)3.0*r[2][j]);
        x[1][j] = (double)1.0/(double)6.0*((double)-2.0*mul1*r[0][j] - mul2*r[1][j] + (double)3.0*r[2][j]);
        x[2][j] = (double)1.0/(double)3.0*(-mul1*r[0][j] + mul2*r[1][j]);
    }
    return;
}

double x_to_R(double x[3][3]) {
    double r1[3], rho1[3];
    for (int j = 0; j < 3; j++) {
        r1[j] = x[1][j] - x[2][j];
        rho1[j] = x[0][j] - (x[1][j]+x[2][j])/((double)2);
    }
    return sqrt(dot(3, r1, r1) + (double)4.0/(double)3.0*dot(3, rho1, rho1));
}

double x_to_R2(double x[3][3]) {
    double r1[3], rho1[3];
    for (int j = 0; j < 3; j++) {
        r1[j] = x[1][j] - x[2][j];
        rho1[j] = x[0][j] - (x[1][j]+x[2][j])/((double)2);
    }
    return dot(3, r1, r1) + (double)4.0/(double)3.0*dot(3, rho1, rho1);
}
