#include "IO.h"

void load_program_config(program_configuration& params, string filename) {
    ifstream in(filename);
    if(!in.is_open()) { ERR_FUNC_CRASH_FILE("Error opening file ") };
    if (!(in >> params.Nx >> params.Ny >> params.Lx >> params.Ly >> params.maxiter >> params.eps >> params.TestId)) { ERR_FUNC_CRASH_FILE("Error reading file ") };
    in.close();
}

void configure_boundaries(boundary_conditions_parameters& params, program_configuration& cfg) {
    // Сброс в дефолтные значения (на всякий случай)
    params.left = params.right = params.up = params.down = 1; // Дирихле
    params.left_value = params.right_value = params.top_value = params.bottom_value = 0.0;
    params.robin_left_value = params.robin_right_value = params.robin_top_value = params.robin_bottom_value = 0.0;

    if (cfg.TestId == 1) { // u = x^2
        params.left = 1;   params.left_value = 0.0;
        params.right = 1;  params.right_value = cfg.Lx * cfg.Lx;
        params.up = 2;     params.top_value = 0.0; // Нейман 0
        params.down = 2;   params.bottom_value = 0.0;
    }
    else if (cfg.TestId == 2) { // u = y^2
        params.left = 2;   params.left_value = 0.0;
        params.right = 2;  params.right_value = 0.0;
        params.up = 1;     params.top_value = cfg.Ly * cfg.Ly;
        params.down = 1;   params.bottom_value = 0.0;
    }
    else if (cfg.TestId == 3 || cfg.TestId == 4) { // Полином x^2+y^2 или Тест сходимости
        // Везде Дирихле 0.0 (реальные значения считаются внутри apply_dirichlet по формуле)
        params.left = 1; params.right = 1; params.up = 1; params.down = 1;
    }
    else {
        // Если это не тест, а основная задача - читаем из файла
        load_boundary_conditions(params, "boundary_conditions.txt");
    }
}

void load_boundary_conditions(boundary_conditions_parameters& params, string filename) {
    ifstream in(filename);
    if(!in.is_open()) { ERR_FUNC_CRASH_FILE("Error opening file ") }

    if (!(in >> params.left >> params.left_value)) { ERR_FUNC_CRASH_FILE("Error reading file ") }
    if (params.left == 3) if (!(in >> params.robin_left_value))  { ERR_FUNC_CRASH_FILE("Error reading file ") }

    if (!(in >> params.up >> params.top_value)) { ERR_FUNC_CRASH_FILE("Error reading file ") }
    if (params.up == 3) if (!(in >> params.robin_top_value))  { ERR_FUNC_CRASH_FILE("Error reading file ") }

    if (!(in >> params.right >> params.right_value)) { ERR_FUNC_CRASH_FILE("Error reading file ") }
    if (params.right == 3) if (!(in >> params.robin_right_value))  { ERR_FUNC_CRASH_FILE("Error reading file ") }

    if (!(in >> params.down >> params.bottom_value))  { ERR_FUNC_CRASH_FILE("Error reading file ") }
    if (params.down == 3) if (!(in >> params.robin_bottom_value))  { ERR_FUNC_CRASH_FILE("Error reading file ") }

    in.close();
}