#include "IO.h"

void load_program_config(program_configuration& params, string filename) {
    ifstream in(filename);
    if(!in.is_open()) { ERR_FUNC_CRASH_FILE("Error opening file ") };
    if (!(in >> params.Nx >> params.Ny >> params.Lx >> params.Ly >> params.maxiter >> params.eps >> params.TestId)) { ERR_FUNC_CRASH_FILE("Error reading file ") };
    in.close();
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