#include "Global.h"
#include "BasisFunctions.h"
#include <vector>
#include <set>

#include "LinearAlgebra.h"

double get_f(double x, double y, int test_id) {

    if (test_id == 1 || test_id == 2) return -2.0;
    if (test_id == 3) return -4.0;

    return 1.0;
}

double get_exact_u(double x, double y, int test_id) {
    if (test_id == 1) return x * x;
    if (test_id == 2) return y * y;

    // ДОБАВИЛИ ЭТО:
    if (test_id == 3) return x * x + y * y;

    return 0.0;
}

void build_local_b(double* local_b, program_configuration& cfg, int ex, int ey) {
    double h_x = cfg.Lx / cfg.Nx;
    double h_y = cfg.Ly / cfg.Ny;

    double x_start = ex * h_x;
    double y_start = ey * h_y;

    for(int k=0; k<9; k++) local_b[k] = 0.0;

    double det_J = (h_x * h_y) / 4.0;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {

            double ksi = GAUSS_POINTS[j];
            double eta = GAUSS_POINTS[i];
            double weight = GAUSS_WEIGHTS[i] * GAUSS_WEIGHTS[j];

            double global_x = x_start + (ksi + 1.0) * h_x / 2.0;
            double global_y = y_start + (eta + 1.0) * h_y / 2.0;
            double f_val = get_f(global_x, global_y, cfg.TestId);

            for (int k = 0; k < 9; k++) {

                double psi = calc_basis_9node(k, ksi, eta);

                local_b[k] += f_val * psi * weight * det_J;
            }
        }
    }
}

void element_fill(Element *Elements_all, program_configuration& cfg) {
    int Nx = cfg.Nx;
    int Ny = cfg.Ny;
    int NodesX = 2 * Nx + 1;

    for (int ey = 0; ey < Ny; ey++) {
        for (int ex = 0; ex < Nx; ex++) {

            int elem_idx = ey * Nx + ex;

            for (int local_x = 0; local_x < 3; local_x++) {
                for (int local_y = 0; local_y < 3; local_y++) {

                    int global_x = ex * 2 + local_x;
                    int global_y = ey * 2 + local_y;

                    int global_node_id = global_y * NodesX + global_x;

                    int local_id = local_x * 3 + local_y;

                    Elements_all[elem_idx].nodes[local_id] = global_node_id;
                }
            }
        }
    }
}

Element& get_element_nodes(Element *&Elements_all, program_configuration& cfg, int ex, int ey) {
    return Elements_all[ey * cfg.Nx + ex];
}

void fill_lambda(double* GlobalLambda, program_configuration& cfg) {
    int NodesX = 2 * cfg.Nx + 1;
    int NodesY = 2 * cfg.Ny + 1;

    for (int j = 0; j < NodesY; j++) {
        for (int i = 0; i < NodesX; i++) {
            int idx = j * NodesX + i;

            if (i < NodesX / 2) {
                GlobalLambda[idx] = 1.0;
            } else {
                GlobalLambda[idx] = 10.0;
            }
        }
    }
}

void build_local_matrix(double **&mtrx, double* lambda_values, double* gamma_values, program_configuration& cfg) {
    double h_x = cfg.Lx / cfg.Nx;
    double h_y = cfg.Ly / cfg.Ny;

    double det_J = (h_x*h_y)/4.0;
    double mult_x = 2.0/h_x;
    double mult_y = 2.0/h_y;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) { //здесь фиксируем точку над которой дрон
            double ksi_pos = GAUSS_POINTS[j];
            double eta_pos = GAUSS_POINTS[i];

            //материал в текущей точке дрона
            double lambda_sum = 0.0;
            double gamma_sum  = 0.0;
            for (int k = 0; k < 4; k++) {
                double weight_func = calc_basis_4node(k, ksi_pos, eta_pos);
                lambda_sum += lambda_values[k] * weight_func;
                gamma_sum  += gamma_values[k]  * weight_func;
            }

            //наклоны в текущей точке дрона
            double grad_x[9], grad_y[9];
            double psi[9];
            for (int k = 0; k < 9; k++) {
                grad_x[k] = calc_grad_ksi_9node(k, ksi_pos, eta_pos)*mult_x;
                grad_y[k] = calc_grad_eta_9node(k, ksi_pos, eta_pos)*mult_y;
                psi[k]    = calc_basis_9node(k, ksi_pos, eta_pos);
            }

            double weight = GAUSS_WEIGHTS[i] * GAUSS_WEIGHTS[j];
            for (int u = 0; u < 9; u++) {
                for (int v = 0; v < 9; v++) {
                    double stiffness = (grad_x[u] * grad_x[v] + grad_y[u] * grad_y[v]) * lambda_sum;

                    double mass = psi[u] * psi[v] * gamma_sum;

                    mtrx[u][v] += (stiffness + mass) * weight * det_J;
                }
            }

        }
    }

}

void assembly(CSRMatrix &A, Element *Elements_all, program_configuration& cfg, double* GlobalLambda, double* GlobalGamma) {
    int Nx = cfg.Nx;
    int Ny = cfg.Ny;

    double** local_matrix = new double*[9];
    for(int i=0; i<9; i++) local_matrix[i] = new double[9];

    for (int ey = 0; ey < Ny; ey++) {
        for (int ex = 0; ex < Nx; ex++) {

            for(int i=0; i<9; i++)
                for(int j=0; j<9; j++) local_matrix[i][j] = 0.0;

            Element elem = get_element_nodes(Elements_all, cfg, ex, ey);

            double lambda_local[4];
            double gamma_local[4];

            int corners[4] = {0, 2, 6, 8};

            for(int k=0; k<4; k++) {
                int global_node = elem.nodes[corners[k]];
                lambda_local[k] = GlobalLambda[global_node];
                gamma_local[k]  = GlobalGamma[global_node];
            }

            build_local_matrix(local_matrix, lambda_local, gamma_local, cfg);

            for (int u = 0; u < 9; u++) {
                for (int v = 0; v < 9; v++) {

                    int row = elem.nodes[u];
                    int col = elem.nodes[v];
                    double val = local_matrix[u][v];

                    if (row == col)
                        A.di[row] += val;

                    if (col < row) {

                        int start = A.ig[row];
                        int end   = A.ig[row+1];

                        for (int k = start; k < end; k++)
                            if (A.jg[k] == col) {
                                A.gg[k] += val;
                                break;
                            }

                    }
                }
            }
        }
    }

    for(int i=0; i<9; i++) delete[] local_matrix[i];
    delete[] local_matrix;
}

void assembly_b(double* Global_B, Element *Elements_all, program_configuration& cfg) {
    int Nx = cfg.Nx;
    int Ny = cfg.Ny;

    double local_b[9];

    for (int ey = 0; ey < Ny; ey++) {
        for (int ex = 0; ex < Nx; ex++) {

            Element elem = get_element_nodes(Elements_all, cfg, ex, ey);

            build_local_b(local_b, cfg, ex, ey);

            for (int k = 0; k < 9; k++) {
                int global_id = elem.nodes[k];
                Global_B[global_id] += local_b[k];
            }
        }
    }
}

void apply_neumann(double* b, Element *Elements_all, program_configuration& cfg, BoundarySide side, double theta) {
    double h_x = cfg.Lx / cfg.Nx;
    double h_y = cfg.Ly / cfg.Ny;

    // Настраиваем параметры в зависимости от стороны
    int start_ex, end_ex, start_ey, end_ey;
    double det_J_1D;
    int local_nodes[3]; // Какие локальные узлы лежат на грани
    bool integrate_by_ksi; // true - бежим по ksi, false - по eta
    double fixed_coord;    // Значение фиксированной координаты

    if (side == SIDE_BOTTOM) {
        start_ex = 0; end_ex = cfg.Nx; start_ey = 0; end_ey = 1;
        det_J_1D = h_x / 2.0;
        local_nodes[0]=0; local_nodes[1]=3; local_nodes[2]=6;
        integrate_by_ksi = true; fixed_coord = -1.0; // eta fixed
    }
    else if (side == SIDE_TOP) {
        start_ex = 0; end_ex = cfg.Nx; start_ey = cfg.Ny - 1; end_ey = cfg.Ny;
        det_J_1D = h_x / 2.0;
        local_nodes[0]=2; local_nodes[1]=5; local_nodes[2]=8;
        integrate_by_ksi = true; fixed_coord = 1.0; // eta fixed
    }
    else if (side == SIDE_LEFT) {
        start_ex = 0; end_ex = 1; start_ey = 0; end_ey = cfg.Ny;
        det_J_1D = h_y / 2.0;
        local_nodes[0]=0; local_nodes[1]=1; local_nodes[2]=2;
        integrate_by_ksi = false; fixed_coord = -1.0; // ksi fixed
    }
    else { // SIDE_RIGHT
        start_ex = cfg.Nx - 1; end_ex = cfg.Nx; start_ey = 0; end_ey = cfg.Ny;
        det_J_1D = h_y / 2.0;
        local_nodes[0]=6; local_nodes[1]=7; local_nodes[2]=8;
        integrate_by_ksi = false; fixed_coord = 1.0; // ksi fixed
    }

    // Универсальный цикл по граничным элементам
    for (int ey = start_ey; ey < end_ey; ey++) {
        for (int ex = start_ex; ex < end_ex; ex++) {

            Element elem = get_element_nodes(Elements_all, cfg, ex, ey);

            // Интегрируем по 3 точкам Гаусса
            for (int j = 0; j < 3; j++) {
                double var_coord = GAUSS_POINTS[j];
                double weight = GAUSS_WEIGHTS[j];

                double ksi = integrate_by_ksi ? var_coord : fixed_coord;
                double eta = integrate_by_ksi ? fixed_coord : var_coord;

                for (int k = 0; k < 3; k++) {
                    int loc_id = local_nodes[k];
                    int glob_id = elem.nodes[loc_id];
                    double psi = calc_basis_9node(loc_id, ksi, eta);

                    b[glob_id] += theta * psi * weight * det_J_1D;
                }
            }
        }
    }
}

void apply_robin(CSRMatrix &A, double* b, Element *Elements_all, program_configuration& cfg, BoundarySide side, double beta, double u_env) {
    double h_x = cfg.Lx / cfg.Nx;
    double h_y = cfg.Ly / cfg.Ny;

    int start_ex, end_ex, start_ey, end_ey;
    double det_J_1D;
    int local_nodes[3];
    bool integrate_by_ksi;
    double fixed_coord;

    // ТА ЖЕ САМАЯ НАСТРОЙКА, что и в Неймане (копипаст логики выбора стороны)
    if (side == SIDE_BOTTOM) {
        start_ex = 0; end_ex = cfg.Nx; start_ey = 0; end_ey = 1;
        det_J_1D = h_x / 2.0;
        local_nodes[0]=0; local_nodes[1]=3; local_nodes[2]=6;
        integrate_by_ksi = true; fixed_coord = -1.0;
    } else if (side == SIDE_TOP) {
        start_ex = 0; end_ex = cfg.Nx; start_ey = cfg.Ny - 1; end_ey = cfg.Ny;
        det_J_1D = h_x / 2.0;
        local_nodes[0]=2; local_nodes[1]=5; local_nodes[2]=8;
        integrate_by_ksi = true; fixed_coord = 1.0;
    } else if (side == SIDE_LEFT) {
        start_ex = 0; end_ex = 1; start_ey = 0; end_ey = cfg.Ny;
        det_J_1D = h_y / 2.0;
        local_nodes[0]=0; local_nodes[1]=1; local_nodes[2]=2;
        integrate_by_ksi = false; fixed_coord = -1.0;
    } else { // RIGHT
        start_ex = cfg.Nx - 1; end_ex = cfg.Nx; start_ey = 0; end_ey = cfg.Ny;
        det_J_1D = h_y / 2.0;
        local_nodes[0]=6; local_nodes[1]=7; local_nodes[2]=8;
        integrate_by_ksi = false; fixed_coord = 1.0;
    }

    for (int ey = start_ey; ey < end_ey; ey++) {
        for (int ex = start_ex; ex < end_ex; ex++) {
            Element elem = get_element_nodes(Elements_all, cfg, ex, ey);

            for (int gp = 0; gp < 3; gp++) {
                double var = GAUSS_POINTS[gp];
                double weight = GAUSS_WEIGHTS[gp];
                double ksi = integrate_by_ksi ? var : fixed_coord;
                double eta = integrate_by_ksi ? fixed_coord : var;

                // 1. Вектор b
                for (int k = 0; k < 3; k++) {
                    int glob_id = elem.nodes[local_nodes[k]];
                    double psi = calc_basis_9node(local_nodes[k], ksi, eta);
                    b[glob_id] += beta * u_env * psi * weight * det_J_1D;
                }

                // 2. Матрица A
                for (int u_idx = 0; u_idx < 3; u_idx++) {
                    for (int v_idx = 0; v_idx < 3; v_idx++) {
                        int loc_u = local_nodes[u_idx];
                        int loc_v = local_nodes[v_idx];
                        int row = elem.nodes[loc_u];
                        int col = elem.nodes[loc_v];

                        double val = beta * calc_basis_9node(loc_u, ksi, eta) * calc_basis_9node(loc_v, ksi, eta) * weight * det_J_1D;

                        if (row == col) A.di[row] += val;
                        else if (col < row) {
                            for (int k = A.ig[row]; k < A.ig[row+1]; k++) {
                                if (A.jg[k] == col) { A.gg[k] += val; break; }
                            }
                        }
                    }
                }
            }
        }
    }
}

void apply_dirichlet(CSRMatrix &A, double* b, program_configuration& cfg, BoundarySide side, double val) {
    int NodesX = 2 * cfg.Nx + 1;
    int NodesY = 2 * cfg.Ny + 1;
    double Penalty = 1.0e10; // (или 1.0e14)

    // Вычисляем шаг сетки (расстояние между соседними узлами)
    double hx_node = cfg.Lx / (NodesX - 1);
    double hy_node = cfg.Ly / (NodesY - 1);

    // Лямбда-функция для фиксации
    auto fix = [&](int node_idx, double current_x, double current_y) {

        double value_to_set = val; // По умолчанию берем то, что передали (константу)

        // АВТОМАТИКА: Если включен Тест 3, считаем точное значение сами
        if (cfg.TestId == 3) {
            value_to_set = get_exact_u(current_x, current_y, cfg.TestId);
        }

        A.di[node_idx] += Penalty;
        b[node_idx]    += Penalty * value_to_set;
    };

    if (side == SIDE_BOTTOM) {
        for (int i = 0; i < NodesX; i++)
            fix(i, i * hx_node, 0.0); // Передаем координаты: (x, 0)
    }
    else if (side == SIDE_TOP) {
        int start_node = (NodesY - 1) * NodesX;
        for (int i = 0; i < NodesX; i++)
            fix(start_node + i, i * hx_node, cfg.Ly); // Координаты: (x, Ly)
    }
    else if (side == SIDE_LEFT) {
        for (int j = 0; j < NodesY; j++)
            fix(j * NodesX, 0.0, j * hy_node); // Координаты: (0, y)
    }
    else if (side == SIDE_RIGHT) {
        for (int j = 0; j < NodesY; j++)
            fix(j * NodesX + NodesX - 1, cfg.Lx, j * hy_node); // Координаты: (Lx, y)
    }
}

void build_matrix_portrait(CSRMatrix &A, Element *Elements_all, program_configuration& cfg) {
    int Nx = cfg.Nx;
    int Ny = cfg.Ny;

    int NodesX = 2 * Nx + 1;
    int NodesY = 2 * Ny + 1;
    int TotalNodes = NodesX * NodesY;
    vector<set<int>> adj(TotalNodes);

    for (int ey = 0; ey < Ny; ey++) {
        for (int ex = 0; ex < Nx; ex++) {

            Element elem = get_element_nodes(Elements_all, cfg, ex, ey);

            for (int i = 0; i < 9; i++) {
                for (int j = 0; j < 9; j++) {

                    int row_node = elem.nodes[i]; //в какую строку пишем
                    int col_node = elem.nodes[j]; //какой элемент в строке

                    if (col_node < row_node) {
                        adj[row_node].insert(col_node);
                    }
                }
            }
        }
    }

    A.ig = new int[TotalNodes + 1];
    A.ig[0] = 0;
    for (int i = 1; i < TotalNodes + 1; i++)
        A.ig[i] = A.ig[i - 1] + adj[i - 1].size();

    A.jg = new int[A.ig[TotalNodes]];

    int j_global = 0;
    for (int i = 0; i < TotalNodes; i++) {
        for (int j : adj[i]) {
            A.jg[j_global] = j;
            j_global++;
        }
    }

    int size_gg = A.ig[TotalNodes];
    A.gg = new double[size_gg];
    A.di = new double[TotalNodes];

    for (int i = 0; i < size_gg; i++) A.gg[i] = 0.0;
    for (int i = 0; i < TotalNodes; i++) A.di[i] = 0.0;
}

void apply_boundary_conditions(CSRMatrix &A, double *&f, Element *&Elements_all, program_configuration &cfg, boundary_conditions_parameters &params) {
    if (params.left == 1)      apply_dirichlet(A, f, cfg, SIDE_LEFT, params.left_value);
    else if (params.left == 2) apply_neumann(f, Elements_all, cfg, SIDE_LEFT, params.left_value);
    else                       apply_robin(A, f, Elements_all, cfg, SIDE_LEFT, params.left_value, params.robin_left_value);

    if (params.right == 1)      apply_dirichlet(A, f, cfg, SIDE_RIGHT, params.right_value);
    else if (params.right == 2) apply_neumann(f, Elements_all, cfg, SIDE_RIGHT, params.right_value);
    else                        apply_robin(A, f, Elements_all, cfg, SIDE_RIGHT, params.right_value, params.robin_right_value);

    if (params.up == 1)      apply_dirichlet(A, f, cfg, SIDE_TOP, params.top_value);
    else if (params.up == 2) apply_neumann(f, Elements_all, cfg, SIDE_TOP, params.top_value);
    else                     apply_robin(A, f, Elements_all, cfg, SIDE_TOP, params.top_value, params.robin_top_value);

    if (params.down == 1)      apply_dirichlet(A, f, cfg, SIDE_BOTTOM, params.bottom_value);
    else if (params.down == 2) apply_neumann(f, Elements_all, cfg, SIDE_BOTTOM, params.bottom_value);
    else                       apply_robin(A, f, Elements_all, cfg, SIDE_BOTTOM, params.bottom_value, params.robin_bottom_value);
}