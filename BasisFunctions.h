#pragma once
#include "Global.h"

extern const double GAUSS_POINTS[3];
extern const double GAUSS_WEIGHTS[3];

double calc_basis_9node(int k, double ksi, double eta);
double calc_grad_ksi_9node(int k, double ksi, double eta);
double calc_grad_eta_9node(int k, double ksi, double eta);
double calc_basis_4node(int k, double ksi, double eta);