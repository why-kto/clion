#pragma once
#include "Global.h"

// КВАДРАТУРЫ ГАУССА

extern const double GAUSS_POINTS[3];  // координаты узлов
extern const double GAUSS_WEIGHTS[3]; // веса узлов
// extern означает буквально "cами массивы с числами лежат в .cpp файле, а здесь только обещаем, что они существуют"
// это нужно, чтобы при подключении .h файла в разные места не создавались дубликаты массивов.


// ПРОТОТИПЫ ФУНКЦИЙ

// для температуры (Биквадратичные, 9 узлов)
// возвращает значение функции формы для узла k в точке (ksi, eta)
double calc_basis_9node(int k, double ksi, double eta);

// возвращает градиенты (наклоны) по осям ksi и eta
// это нужно для Матрицы Жесткости (перемножение градиентов)
double calc_grad_ksi_9node(int k, double ksi, double eta);
double calc_grad_eta_9node(int k, double ksi, double eta);

// для материала (коэффициент диффузии аппроксимируется билинейно, поэтому 4 узла)
double calc_basis_4node(int k, double ksi, double eta);