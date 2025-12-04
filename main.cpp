#include <iomanip>

#include "Global.h"
#include "IO.h"
#include "FEM.h"

int main() {

    program_configuration cfg;
    boundary_conditions_parameters params;

    load_program_config(cfg, "kuslau.txt");

    if (cfg.TestId == 1) {
        cout << ">>> MODE: TEST 1 (u = x^2) <<<" << endl;

        // Слева (1) и Справа (3) - Дирихле (значения x^2)
        // Сверху (2) и Снизу (4) - Нейман 0 (изоляция)
        params.left = 1;   params.left_value = 0.0;  // x=0 -> u=0
        params.right = 1;  params.right_value = cfg.Lx * cfg.Lx; // x=Lx -> u=Lx^2
        params.up = 2;     params.top_value = 0.0;
        params.down = 2;   params.bottom_value = 0.0;
    }
    else if (cfg.TestId == 2) {
        cout << ">>> MODE: TEST 2 (u = y^2) <<<" << endl;
        // Настраиваем границы:
        // Слева (1) и Справа (3) - Нейман 0
        // Сверху (2) и Снизу (4) - Дирихле (значения y^2)
        params.left = 2;   params.left_value = 0.0;
        params.right = 2;  params.right_value = 0.0;
        params.up = 1;     params.top_value = cfg.Ly * cfg.Ly; // y=Ly -> u=Ly^2
        params.down = 1;   params.bottom_value = 0.0; // y=0 -> u=0
    }
    else if (cfg.TestId == 3) {
        cout << ">>> MODE: TEST 3 (u = x^2 + y^2) <<<" << endl;
        // Все границы - Дирихле (Тип 1).
        // Значения (0.0) тут не важны, так как внутри apply_dirichlet
        // сработает наша новая логика и посчитает формулу.
        params.left = 1;   params.left_value = 0.0;
        params.right = 1;  params.right_value = 0.0;
        params.up = 1;     params.top_value = 0.0;
        params.down = 1;   params.bottom_value = 0.0;
    }
    else if (cfg.TestId == 4) {
        cout << ">>> MODE: TEST 4 (Convergence Study) <<<" << endl;
        // Все границы - Дирихле u=0
        params.left = 1;   params.left_value = 0.0;
        params.right = 1;  params.right_value = 0.0;
        params.up = 1;     params.top_value = 0.0;
        params.down = 1;   params.bottom_value = 0.0;
    }
    else load_boundary_conditions(params, "boundary_conditions.txt");

    int Nx = cfg.Nx;
    int Ny = cfg.Ny;
    int NodesX = 2 * Nx + 1;
    int NodesY = 2 * Ny + 1;
    int TotalNodes = NodesX * NodesY;
    int TotalElements = Nx * Ny;

    cout << "FEM Analysis Start" << endl;
    cout << "Grid: " << Nx << "x" << Ny << " Elements (" << TotalNodes << " Nodes)" << endl;

    // ------------------------------------------
    // 2. ГЕОМЕТРИЯ (Сетка)
    // ------------------------------------------
    Element* Elements_all = new Element[TotalElements];
    element_fill(Elements_all, cfg);

    // ------------------------------------------
    // 3. ПОРТРЕТ МАТРИЦЫ (Скелет)
    // ------------------------------------------
    CSRMatrix A;
    build_matrix_portrait(A, Elements_all, cfg);

    // ------------------------------------------
    // 4. ПАМЯТЬ ПОД ВЕКТОРА (Твоя функция)
    // ------------------------------------------
    WorkVectors extra;
    double* x = nullptr;
    double* f = nullptr; // Это будет наша правая часть (вектор b)

    allocWorkspace(cfg, extra, x, f);

    // ------------------------------------------
    // 5. СБОРКА СИСТЕМЫ (Физика)
    // ------------------------------------------


    double* GlobalLambda = new double[TotalNodes];

    // if (cfg.TestId > 0) {
    //     // Для тестов Лямбда всегда 1.0
    //     for(int i=0; i<TotalNodes; i++) GlobalLambda[i] = 1.0;
    // } else {
    //     // Для обычного режима - твоя функция с распределением (1.0 и 10.0)
    //     fill_lambda(GlobalLambda, cfg);
    // }
    for(int i=0; i<TotalNodes; i++) GlobalLambda[i] = 1.0;

    double* GlobalGamma = new double[TotalNodes];

    // if (cfg.TestId > 0 && cfg.TestId < 4) {
    //     // Для тестов Гамма всегда 0.0 кроме теста 4
    //     for(int i=0; i<TotalNodes; i++) GlobalGamma[i] = 0.0;
    // } else if (cfg.TestId == 4) {
    //     //тест 4
    //     for(int i=0; i<TotalNodes; i++) GlobalGamma[i] = 1.0;
    // }
    for(int i=0; i<TotalNodes; i++) GlobalGamma[i] = 0.0;

    assembly(A, Elements_all, cfg, GlobalLambda, GlobalGamma);

    // Сборка Вектора f (правая часть)
    // Важно: сначала обнулить f, так как assembly_b делает +=
    for(int i=0; i<TotalNodes; i++) f[i] = 0.0;
    assembly_b(f, Elements_all, cfg);

    // ------------------------------------------
    // 6. КРАЕВЫЕ УСЛОВИЯ (Закрепляем края)
    // ------------------------------------------

    apply_boundary_conditions(A, f, Elements_all, cfg, params);

    // ------------------------------------------
    // 7. РЕШЕНИЕ (Солвер)
    // ------------------------------------------

    // Матрица S для предобуславливателя
    CSRMatrix S;
    // Строим неполное разложение Холецкого
    buildICFactor(cfg, A, S);

    // Начальное приближение x0 = 0
    for(int i=0; i<TotalNodes; i++) x[i] = 0.0;

    cout << "Solving system..." << endl;
    int iters = solve_MSG_IC(cfg, A, S, f, x, extra);

    cout << "Done in " << iters << " iterations." << endl;

    // ------------------------------------------
    // 8. ВЫВОД РЕЗУЛЬТАТА
    // ------------------------------------------
    cout << "\nTemperature Field:" << endl;
    cout.precision(2); // 2 знака после запятой
    cout << fixed;

    // Выводим сеткой, чтобы было наглядно (снизу вверх)
    for (int j = NodesY - 1; j >= 0; j--) { // Строки (Y)
        for (int i = 0; i < NodesX; i++) {  // Столбцы (X)
            int node_idx = j * NodesX + i;
            // Вывод числа с выравниванием
            if (abs(x[node_idx]) < 1e-10) cout << "0.00\t";
            else cout << x[node_idx] << "\t";
        }
        cout << endl;
    }

    if (cfg.TestId == 4) cout << endl << scientific << setprecision(15) << abs(x[NodesY * (NodesX / 2) + NodesY / 2]) - 1.0 << endl;
    system("pause");
    return 0;
}