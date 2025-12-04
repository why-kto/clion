#include "Global.h"
#include <cstring>

//222

void allocWorkspace(program_configuration& cfg, WorkVectors& extra, double*& x, double*& f) {
    int NodesX = 2 * cfg.Nx + 1;
    int NodesY = 2 * cfg.Ny + 1;
    int N = NodesX * NodesY;

    x = new double[N];
    f = new double[N];

    extra.r = new double[N];
    extra.z = new double[N];
    extra.p = new double[N];
}

double scal(int N, double* vec1, double* vec2) {
    double result = 0.0;
    for (int i = 0; i < N; i++) result += vec1[i] * vec2[i];
    return result;
}

double dotProduct(double*& v1, double*& v2, int N) {
    double sum = 0.0;
    for(int i = 0; i < N; i++) sum += v1[i] * v2[i];
    return sum;
}

double vectorNorm(double*& v, int N) {
    return sqrt(dotProduct(v, v, N));
}

void addScaledVector(int N, double alpha, double* x, double* y) {
    for (int i = 0; i < N; ++i)
        y[i] += alpha * x[i];
}

void addScaledVector(int N, double alpha, double* x) {
    for (int i = 0; i < N; ++i)
        x[i] *= alpha;
}

void mulSymCSRMatVec(program_configuration& cfg, CSRMatrix& A, double*& x, double*& f) {
    int NodesX = 2 * cfg.Nx + 1;
    int NodesY = 2 * cfg.Ny + 1;
    int N = NodesX * NodesY;

    double* di = A.di;
    double* gg = A.gg;
    int* ig = A.ig;
    int* jg = A.jg;

    for (int i = 0; i < N; ++i) {
        f[i] = di[i] * x[i];
    }

    for (int row = 0; row < N; row++)
    {
        for (int col = ig[row]; col < ig[row + 1]; ++col)
        {
            double value = gg[col];
            int j = jg[col];

            f[row] += value * x[j];
            f[j] += value * x[row];
        }
    }
}

void buildICFactor(program_configuration& cfg, CSRMatrix& A, CSRMatrix& S) {
    int NodesX = 2 * cfg.Nx + 1;
    int NodesY = 2 * cfg.Ny + 1;
    int N = NodesX * NodesY;

    S.ig = A.ig;
    S.jg = A.jg;
    int M = A.ig[N];

    S.di = new double[N];
    S.gg = new double[M];
    for(int i=0; i<N; i++) S.di[i] = 0.0;
    for(int i=0; i<M; i++) S.gg[i] = 0.0;

    double* adi = A.di;
    double* agg = A.gg;
    double* sdi = S.di;
    double* sgg = S.gg;
    int*    ig  = A.ig;
    int*    jg  = A.jg;

    for (int i = 0; i < N; ++i) {

        double sum_d = 0.0; // вклад диагонали

        // обновляем внедиагональные элементы + сразу считаем вклад диагоналей
        for (int k = ig[i]; k < ig[i + 1]; ++k) {

            int j = jg[k]; // столбец j < i для элемента [i,j]

            double sum = 0.0;

            // проходим по строкам i и j и ищем совпадающие столбцы m < j
            int ki = ig[i]; // указатель по строке i
            int kj = ig[j]; // указатель по строке j

            int end_i = ig[i + 1];
            int end_j = ig[j + 1];

            while (true) {
                int col_i = jg[ki];
                int col_j = jg[kj];

                if (col_i == col_j) {  // найден общий столбец m = col_i
                    sum += sgg[ki] * sgg[kj]; ++ki; ++kj;
                    if (ki >= end_i || kj >= end_j) break;
                }
                else if (col_i < col_j) {
                    ++ki; if (ki >= end_i) break;
                }
                else {
                    ++kj; if (kj >= end_j) break;
                }
            }

            // формула: S_ij = (A_ij - sum) / S_jj
            sgg[k] = (agg[k] - sum) / sdi[j];

            double sij = sgg[k]; //вклад диагонали
            sum_d += sij * sij;  //вклад диагонали
        }

        // обновляем диагональ
        double diag_val = adi[i] - sum_d;
        if (diag_val <= 0.0) ERR_FUNC_CRASH("readed matrix is not SPD"); //защита от отрицательного и нулевого
        sdi[i] = sqrt(diag_val);

    }

}

void applyICPreconditioner(program_configuration& cfg, CSRMatrix& S, double* r, double* z, double* y) {
    int NodesX = 2 * cfg.Nx + 1;
    int NodesY = 2 * cfg.Ny + 1;
    int N = NodesX * NodesY;

    double* di = S.di;
    double* gg = S.gg;
    int*    ig = S.ig;
    int*    jg = S.jg;

    // 1) прямой ход: S * y = r
    for (int i = 0; i < N; ++i) {
        double sum = 0.0;
        int k_start = ig[i], k_end = ig[i + 1];
        for (int k = k_start; k < k_end; ++k) {
            int j = jg[k];
            sum += gg[k] * y[j];
        }
        y[i] = (r[i] - sum) / di[i];
    }

    // 2) обратный ход: S^T * z = y
    for (int i = N - 1; i >= 0; --i) {
        double zi = y[i] / di[i];
        z[i] = zi;

        int k_start = ig[i], k_end = ig[i + 1];
        for (int k = k_start; k < k_end; ++k) {
            int j = jg[k];
            y[j] -= gg[k] * zi;
        }
    }
}

void computeResidualAxMinusb(program_configuration& cfg, CSRMatrix& A, double* x, double* f, double* r){
    int NodesX = 2 * cfg.Nx + 1;
    int NodesY = 2 * cfg.Ny + 1;
    int N = NodesX * NodesY;

    double* di = A.di;
    double* gg = A.gg;
    int*    ig = A.ig;
    int*    jg = A.jg;

    memcpy(r, f, N * sizeof(f[0]));

    for (int i = 0; i < N; ++i) {
        r[i] -= di[i] * x[i];
    }

    for (int row = 0; row < N; ++row) {
        for (int col = ig[row]; col < ig[row + 1]; ++col) {
            int    j   = jg[col];
            double value = gg[col];

            r[row] -= value * x[j];
            r[j]   -= value * x[row];
        }
    }
}

int solve_MSG_IC(program_configuration& cfg, CSRMatrix& A, CSRMatrix& S, double*& f, double*& x, WorkVectors& extra) {
    int NodesX = 2 * cfg.Nx + 1;
    int NodesY = 2 * cfg.Ny + 1;
    int N = NodesX * NodesY;

    double* r = extra.r; // r_k — невязка (r_k = f - A x_k)
    double* z = extra.z; // z_k — направление спуска
    double* p = extra.p; // p_k — направление
    double* Ap = f;      // Ap_k = A * p_k (используем f как рабочий буфер)

    double norm_f = vectorNorm(f, N); // Норма правой части ||f|| (считаем до порчи f)

    computeResidualAxMinusb(cfg, A, x, f, r); // r0 = f - A*x0

    applyICPreconditioner(cfg, S, r, z, f); // z0 = M^{-1} r0  (IC)

    memcpy(p, z, N * sizeof(z[0])); // p0 = z0
    double rz = dotProduct(r, z, N); // (r0, z0)

    double rel = vectorNorm(r, N) / norm_f;  // относительная невязка ||r0|| / ||f||


    if (rel < cfg.eps) return 0; // если начальное приближение уже хорошее — выходим

    //----------итерационный цикл предобусловленного МСГ----------
    int iter = 1;
    // условие останова в шапке: либо превысили maxiter,
    // либо относительная невязка ещё не меньше eps (формула (3.9))
    for (; iter < cfg.maxiter && rel >= cfg.eps; ++iter) {

        // 1) Ap_k = A * p_k
        mulSymCSRMatVec(cfg, A, p, Ap);

        // 2) alpha_k = (r_k, z_k) / (A p_k, p_k)
        double Ap_p = dotProduct(Ap, p, N);   // (A p_k, p_k)
        double alpha = rz / Ap_p;

        // 3) x_{k+1} = x_k + alpha_k * p_k
        addScaledVector(N, alpha, p, x);         // x += alpha * p

        // 4) r_{k+1} = r_k - alpha_k * A z_k (формула (3.6))
        addScaledVector(N, -alpha, Ap, r);       // r -= alpha * Ap

        // относительная невязка ||r_{k+1}|| / ||f||
        rel = vectorNorm(r, N) / norm_f;

        // 5) z_{k+1} = M^{-1} r_{k+1}  (IC)
        applyICPreconditioner(cfg, S, r, z, f);

        // 6) (r_{k+1}, z_{k+1})
        double rz_new = dotProduct(r, z, N);

        // 7) beta_{k+1} = (r_{k+1}, z_{k+1}) / (r_k, z_k)
        double beta = rz_new / rz;

        // 8) p_{k+1} = z_{k+1} + beta * p_k
        addScaledVector(N, beta, p);         // p = beta * p_k
        addScaledVector(N, 1.0, z, p);       // p = z_{k+1} + beta * p_k

        rz = rz_new;
    }

    return iter;
}