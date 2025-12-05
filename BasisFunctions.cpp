#include "BasisFunctions.h"

// ТОЧКИ И ВЕСА ГАУССА (3x3)

const double GAUSS_POINTS[3]  = { -sqrt(3.0/5.0), 0.0,     sqrt(3.0/5.0) }; // координаты узлов
const double GAUSS_WEIGHTS[3] = {  5.0/9.0,       8.0/9.0, 5.0/9.0 };       // веса узлов


// ВСПОМОГАТЕЛЬНЫЕ ОДНОМЕРНЫЕ ФУНКЦИИ (1D) (нужны для построения 2D функций)

// квадратичные функции
double q1d_m1(double t) { return 0.5 * t * (t - 1.0); } // если t = -1, эта функция дает 1, иначе 0
double q1d_0 (double t) { return 1.0 - t * t; }         // если t =  0, эта функция дает 1, иначе 0
double q1d_p1(double t) { return 0.5 * t * (t + 1.0); } // если t = +1, эта функция дает 1, иначе 0

// производные квадратичных функций (нужны для градиентов)
double q1d_m1_d(double t) { return t - 0.5; }
double q1d_0_d (double t) { return -2.0 * t; }
double q1d_p1_d(double t) { return t + 0.5; }

// линейные функции (для Лямбды)
double lin1d_m1(double t) { return 0.5 * (1.0 - t); } // если t = -1, эта функция дает 1, иначе 0
double lin1d_p1(double t) { return 0.5 * (1.0 + t); } // если t = +1, эта функция дает 1, иначе 0


// ДВУМЕРНЫЕ ФУНКЦИИ ФОРМЫ (9 УЗЛОВ, k - локальный номер узла (0..8), суть в перемножении одномерных)

double calc_basis_9node(int k, double ksi, double eta) {
    switch (k) {
        // нижний ряд (eta = -1)
        case 0: return q1d_m1(ksi)*q1d_m1(eta); // левый нижний
        case 1: return q1d_m1(ksi)*q1d_0(eta);  // нижний центр
        case 2: return q1d_m1(ksi)*q1d_p1(eta); // правый нижний

        // средний ряд (eta = 0)
        case 3: return q1d_0(ksi) *q1d_m1(eta); // левый центр
        case 4: return q1d_0(ksi) *q1d_0(eta);  // центр
        case 5: return q1d_0(ksi) *q1d_p1(eta); // правый центр

        // верхний ряд (eta = +1)
        case 6: return q1d_p1(ksi)*q1d_m1(eta); // левый верхний
        case 7: return q1d_p1(ksi)*q1d_0(eta);  // верхний центр
        case 8: return q1d_p1(ksi)*q1d_p1(eta); // правый верхний
        default: return 0.0;
    }
}


// ГРАДИЕНТЫ

// здесь берем производную только от ksi, а eta оставляем как есть
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

// здесь берем производную только от eta, а ksi оставляем как есть
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


// ДВУМЕРНЫЕ ФУНКЦИИ ДЛЯ МАТЕРИАЛА (4 УЗЛА, используются ТОЛЬКО для интерполяции Лямбды (теплопроводности) и Гаммы)

double calc_basis_4node(int k, double ksi, double eta) {
    switch (k) {
        case 0: return lin1d_m1(ksi)*lin1d_m1(eta); // левый нижний
        case 1: return lin1d_m1(ksi)*lin1d_p1(eta); // левый верхний
        case 2: return lin1d_p1(ksi)*lin1d_m1(eta); // правый нижний
        case 3: return lin1d_p1(ksi)*lin1d_p1(eta); // правый верхний
        default: return 0.0;
    }
}