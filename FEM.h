#include "Global.h"
#include "BasisFunctions.h"
#include "LinearAlgebra.h" // Нужно для CSRMatrix
#include <set>

// Работа с сеткой
void element_fill(Element *Elements_all, program_configuration& cfg);
Element& get_element_nodes(Element *&Elements_all, program_configuration& cfg, int ex, int ey);
void fill_lambda(double* GlobalLambda, program_configuration& cfg);

// Сборка
void build_matrix_portrait(CSRMatrix &A, Element *Elements_all, program_configuration& cfg);
void assembly(CSRMatrix &A, Element *Elements_all, program_configuration& cfg, double* GlobalLambda, double* GlobalGamma);
void assembly_b(double* Global_B, Element *Elements_all, program_configuration& cfg);

// Краевые условия
void apply_boundary_conditions(CSRMatrix &A, double *&f, Element *&Elements_all, program_configuration &cfg, boundary_conditions_parameters &params);