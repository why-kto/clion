#include "Global.h"
#include <fstream>

void load_program_config(program_configuration& params, string filename);
void configure_boundaries(boundary_conditions_parameters& params, program_configuration& cfg);
void load_boundary_conditions(boundary_conditions_parameters& params, string filename);