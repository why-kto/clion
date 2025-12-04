#include "BasisFunctions.h"

const double GAUSS_POINTS[3] = { -sqrt(3.0/5.0), 0.0, sqrt(3.0/5.0) };
const double GAUSS_WEIGHTS[3] = { 5.0/9.0, 8.0/9.0, 5.0/9.0 };

double q1d_m1(double t) { return 0.5 * t * (t - 1.0); }
double q1d_0 (double t) { return 1.0 - t * t; }
double q1d_p1(double t) { return 0.5 * t * (t + 1.0); }

double q1d_m1_d(double t) { return t - 0.5; }
double q1d_0_d (double t) { return -2.0 * t; }
double q1d_p1_d(double t) { return t + 0.5; }

double lin1d_m1(double t) { return 0.5 * (1.0 - t); }
double lin1d_p1(double t) { return 0.5 * (1.0 + t); }

double calc_basis_9node(int k, double ksi, double eta) {
    switch (k) {
        case 0: return q1d_m1(ksi)*q1d_m1(eta);
        case 1: return q1d_m1(ksi)*q1d_0(eta);
        case 2: return q1d_m1(ksi)*q1d_p1(eta);
        case 3: return q1d_0(ksi) *q1d_m1(eta);
        case 4: return q1d_0(ksi) *q1d_0(eta);
        case 5: return q1d_0(ksi) *q1d_p1(eta);
        case 6: return q1d_p1(ksi)*q1d_m1(eta);
        case 7: return q1d_p1(ksi)*q1d_0(eta);
        case 8: return q1d_p1(ksi)*q1d_p1(eta);
        default: return 0.0;
    }
}

double calc_grad_ksi_9node(int k, double ksi, double eta) {
    switch (k) {
        case 0: return q1d_m1_d(ksi)*q1d_m1(eta);
        case 1: return q1d_m1_d(ksi)*q1d_0(eta);
        case 2: return q1d_m1_d(ksi)*q1d_p1(eta);
        case 3: return q1d_0_d(ksi) *q1d_m1(eta);
        case 4: return q1d_0_d(ksi) *q1d_0(eta);
        case 5: return q1d_0_d(ksi) *q1d_p1(eta);
        case 6: return q1d_p1_d(ksi)*q1d_m1(eta);
        case 7: return q1d_p1_d(ksi)*q1d_0(eta);
        case 8: return q1d_p1_d(ksi)*q1d_p1(eta);
        default: return 0.0;
    }
}

double calc_grad_eta_9node(int k, double ksi, double eta) {
    switch (k) {
        case 0: return q1d_m1(ksi)*q1d_m1_d(eta);
        case 1: return q1d_m1(ksi)*q1d_0_d(eta);
        case 2: return q1d_m1(ksi)*q1d_p1_d(eta);
        case 3: return q1d_0(ksi) *q1d_m1_d(eta);
        case 4: return q1d_0(ksi) *q1d_0_d(eta);
        case 5: return q1d_0(ksi) *q1d_p1_d(eta);
        case 6: return q1d_p1(ksi)*q1d_m1_d(eta);
        case 7: return q1d_p1(ksi)*q1d_0_d(eta);
        case 8: return q1d_p1(ksi)*q1d_p1_d(eta);
        default: return 0.0;
    }
}

double calc_basis_4node(int k, double ksi, double eta) {
    switch (k) {
        case 0: return lin1d_m1(ksi)*lin1d_m1(eta);
        case 1: return lin1d_m1(ksi)*lin1d_p1(eta);
        case 2: return lin1d_p1(ksi)*lin1d_m1(eta);
        case 3: return lin1d_p1(ksi)*lin1d_p1(eta);
        default: return 0.0;
    }
}