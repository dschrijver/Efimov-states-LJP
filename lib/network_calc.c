#include <math.h> // exp(), tanh()

#include "network_calc.h"


double calc(double r[3][3]) {
    double u_out = 0;
    xi(r); 
    for (int i = 0; i < N_hid; i++) {
        u_hid[i] = b1[i];
        for (int j = 0; j < N_in; j++) {
            u_hid[i] += w1[j][i] * u_in[j];
        }
        u_out += w2[i] * tanh(u_hid[i]);
    }
    return exp(u_out);
}

void xi(double r[3][3]) {
    double res1 = 0.0, res2 = 0.0, subres2, res3 = 0.0, subres3;
    for (int i = 0; i < 3; i++) {
        res1 += (double) (r[2][i]*r[2][i]);
        subres2 = r[1][i]*sqrt((double)3.0) + r[2][i];
        res2 += (double) (subres2*subres2);
        subres3 = subres2 - (double)2.0*r[2][i];
        res3 += (double) (subres3*subres3);
    }
    u_in[0] = sqrt(res1);
    u_in[1] = (double)0.5 * sqrt(res2);
    u_in[2] = (double)0.5 * sqrt(res3);
    return;
}

void Ci(double r[3][3], int l) {
    xi(r);
    double mul = sqrt((double)3.0);
    C[0][0] = (double) 0.0;
    C[1][0] = r[2][l] / u_in[0];
    C[1][1] = (mul*r[1][l] + r[2][l])/((double)4.0*u_in[1]);
    C[0][1] = mul * C[1][1];
    C[1][2] = -(mul*r[1][l] - r[2][l])/((double)4.0*u_in[2]);
    C[0][2] = -mul * C[1][2];
    return;
}

void Di(double r[3][3], int l) {
    xi(r);
    double mul = sqrt((double)3.0);
    double res;
    D[0][0] = (double) 0.0;
    D[1][0] = -r[2][l]*r[2][l]/(u_in[0]*u_in[0]*u_in[0]) + (double)1.0/u_in[0];
    res = mul*r[1][l] + r[2][l];
    D[1][1] = -(res*res)/((double)16.0*u_in[1]*u_in[1]*u_in[1]) + (double)0.25/u_in[1];
    D[0][1] = (double)3.0 * D[1][1];
    res = mul*r[1][l] - r[2][l];
    D[1][2] = -(res*res)/((double)16.0*u_in[2]*u_in[2]*u_in[2]) + (double)0.25/u_in[2];
    D[0][2] = (double)3.0 * D[1][2];
    return;
}