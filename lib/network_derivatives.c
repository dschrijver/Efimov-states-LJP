#include <math.h> // tanh()

#include "network_calc.h" // calc(), Ci(), Di()
#include "network_derivatives.h"


void dpsi_dw(double r[3][3]) {
    double psi = calc(r);
    double fu_hid, fderiv;

    for (int i = 0; i < N_hid; i++) {
        fu_hid = tanh(u_hid[i]);
        fderiv = (double)1.0-fu_hid*fu_hid;
        dpsidw2[i] = psi * fu_hid;
        dpsidb1[i] = psi * w2[i] * fderiv;
        for (int j = 0; j < N_in; j++) {
            dpsidw1[j][i] = psi * w2[i] * fderiv * u_in[j];
        }
    }
    
    return;
}

void dpsi_dr(double r[3][3]) {
    double psi = calc(r);
    double fu_hid, fderiv, subres;

    for (int l = 0; l < 3; l++) {
        Ci(r, l);
        for (int k = 0; k < 2; k++) {
            dpsidr_list[k][l] = 0.0;
            for (int i = 0; i < N_hid; i++) {
                subres = 0.0;
                for (int j = 0; j < N_in; j++) {
                    subres += w1[j][i] * C[k][j];
                }
                fu_hid = tanh(u_hid[i]);
                fderiv = (double)1.0-fu_hid*fu_hid;
                subres *= w2[i]*fderiv;
                dpsidr_list[k][l] += subres;
            }
            dpsidr_list[k][l] *= psi;
        }
    }

    return;
}

void d2psi_dr2(double r[3][3]) {
    double psi = calc(r);
    dpsi_dr(r);

    double fu_hid, fderiv, f2deriv, sumC, sumD;

    for (int l = 0; l < 3; l++) {
        Ci(r, l);
        Di(r, l);
        for (int k = 0; k < 2; k++) {
            d2psidr2_list[k][l] = 0.0;
            for (int i = 0; i < N_hid; i++) {
                sumC = 0.0;
                sumD = 0.0;
                for (int j = 0; j < N_in; j++) {
                    sumC += w1[j][i] * C[k][j];
                    sumD += w1[j][i] * D[k][j];
                }
                fu_hid = tanh(u_hid[i]);
                fderiv = (double)1.0-fu_hid*fu_hid;
                f2deriv = (double)-2.0*fu_hid + (double)2.0*(fu_hid*fu_hid*fu_hid);
                d2psidr2_list[k][l] += w2[i]*(f2deriv*sumC*sumC + fderiv*sumD);
            }
            d2psidr2_list[k][l] = (double)1.0/psi*(dpsidr_list[k][l]*dpsidr_list[k][l]) + psi*d2psidr2_list[k][l];
        }
    }

    return;
}
