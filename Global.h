#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib> // для exit, system

using namespace std;

#define ERR_FUNC_CRASH(err_text) { cout << err_text << endl; system("pause"); exit(2); }
#define ERR_FUNC_CRASH_FILE(err_text) { cout << err_text << '\"' << filename << '\"' << endl; system("pause"); exit(1); }

struct Element { int nodes[9]; };

struct CSRMatrix {
    double* di;
    double* gg;
    int*    ig;
    int*    jg;

    // Вспомогательный метод для очистки, если нужно будет
    void clean() {
        delete[] di; delete[] gg; delete[] ig; delete[] jg;
    }
};

struct program_configuration {
    int Nx, Ny;
    double Lx, Ly;
    int maxiter;
    double eps;
    double TestId;
};

struct boundary_conditions_parameters {
    int left, up, right, down;
    double left_value, top_value, right_value, bottom_value;
    double robin_left_value, robin_top_value, robin_right_value, robin_bottom_value;
};

struct WorkVectors {
    double* r;
    double* z;
    double* p;
};

enum BoundarySide {
    SIDE_LEFT,
    SIDE_RIGHT,
    SIDE_BOTTOM,
    SIDE_TOP
};