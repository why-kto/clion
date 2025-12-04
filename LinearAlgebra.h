#include "Global.h"

void allocWorkspace(program_configuration& cfg, WorkVectors& extra, double*& x, double*& f);

// Базовые векторные операции
double scal(int N, double* vec1, double* vec2); // у тебя была такая функция, хотя dotProduct дублирует её
double dotProduct(double*& v1, double*& v2, int N);
double vectorNorm(double*& v, int N);
void addScaledVector(int N, double alpha, double* x, double* y); // y += alpha*x
void addScaledVector(int N, double alpha, double* x);            // x *= alpha

// Операции с CSR матрицей
void mulSymCSRMatVec(program_configuration& cfg, CSRMatrix& A, double*& x, double*& f);
void buildICFactor(program_configuration& cfg, CSRMatrix& A, CSRMatrix& S);
void applyICPreconditioner(program_configuration& cfg, CSRMatrix& S, double* r, double* z, double* y);
void computeResidualAxMinusb(program_configuration& cfg, CSRMatrix& A, double* x, double* f, double* r);

// Основной решатель
int solve_MSG_IC(program_configuration& cfg, CSRMatrix& A, CSRMatrix& S, double*& f, double*& x, WorkVectors& extra);