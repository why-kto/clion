#include <iomanip>  // библиотека для вывода таблицы решения (setw, precision)
#include "Global.h" // наши структуры данных (Element, CSRMatrix, Config...)
#include "IO.h"     // функции ввода-вывода (чтение файлов)
#include "FEM.h"    // математическое ядро (сборка, базисные функции, физика)

int main() {

    program_configuration cfg;              // создаем структуру для хранения настроек
    load_program_config(cfg, "kuslau.txt"); // загружаем из файла настройки

    int Nx = cfg.Nx;                  // количество элементов по оси x (читаем из файла)
    int Ny = cfg.Ny;                  // количество элементов по оси y (читаем из файла)
    int NodesX = 2 * Nx + 1;          // количество узлов по оси x; формула: 2 * кол-во элементов по оси x + 1 (для 9-ти узловых элемнтов)
    int NodesY = 2 * Ny + 1;          // количество узлов по оси y; формула: 2 * кол-во элементов по оси y + 1 (для 9-ти узловых элемнтов)
    int TotalNodes = NodesX * NodesY; // общее количество узлов в сетке; формула: количество узлов по оси x * количество узлов по оси y
    int TotalElements = Nx * Ny;      // общее количество элементов в сетке; формула: количество элементов по оси x * количество элементов по оси y

    // выводим информацию, чтобы убедиться, что конфиг считался верно
    cout << "FEM Analysis Start" << endl;
    cout << "Mode: " << (cfg.TestId > 0 ? "TEST" : "MAIN TASK") << " #" << cfg.TestId << endl;
    cout << "Grid: " << Nx << "x" << Ny << " Elements (" << TotalNodes << " Nodes)" << endl;

    Element* Elements_all = new Element[TotalElements]; // выделяем память под массив элементов
    element_fill(Elements_all, cfg);                    // заполняем этот массив (функция связывает глобальные номера узлов с каждым конкретным элементом)

    CSRMatrix A;                                 // создаем структуру для хранения разреженной матрицы
    build_matrix_portrait(A, Elements_all, cfg); // запоняем массивы ig и jg (какой узел с каким соседит) + выделяем память для di и gg

    WorkVectors extra;                // создаем структуру (3 указателя r, z, p) для вспомогательных векторов (для решателя)
    double* x = nullptr;              // создаем указатель на вектор решения (температура)
    double* f = nullptr;              // создаем указатель на вектор правой части (источники тепла)
    allocWorkspace(cfg, extra, x, f); // выделяем память под все вышеперечисленное (r, z, p, x, f)

    double* GlobalLambda = new double[TotalNodes]; // создаем указатель на массив коэффициентов диффузии (теплопроводности) для каждого узла и сразу выделяем память
    fill_lambda(GlobalLambda, cfg);                // заполняем лямбду (массив коэффициентов диффузии) с помощью функции

    double* GlobalGamma = new double[TotalNodes];  // создаем указатель на массив коэффициентов реакции (теплоотдачи/поглощения) для каждого узла и сразу выделяем память
    fill_gamma(GlobalGamma, cfg);                  // заполняем гамму (массив коэффициентов реакции) с помощью функции

    // 1) цикл по всем элементам
    // 2) подсчет локальных матриц (интегрируем градиенты * лямбду + функции * гамму)
    // 3) суммируем вклады в глобальную матрицу A
    assembly(A, Elements_all, cfg, GlobalLambda, GlobalGamma);

    // собираем вектор правой части
    for(int i=0; i<TotalNodes; i++) f[i] = 0.0; // зануляем вектор так как assembly_b делает += (накопление)
    assembly_b(f, Elements_all, cfg);           // интегрируем функцию источника (правую часть уравнения) и добавляем в вектор f

    boundary_conditions_parameters params;      // создаем структуру для краевых условий
    configure_boundaries(params, cfg);          // в зависимости от номера теста заполняем структуру
    apply_boundary_conditions(A, f, Elements_all, cfg, params); // применяем условия

    CSRMatrix S;              // создаем структуру для хранения разреженной матрицы
    buildICFactor(cfg, A, S); // делаем неполное разложение

    for(int i=0; i<TotalNodes; i++) x[i] = 0.0; // начальное приближение для солвера

    cout << "Solving system..." << endl;                   // выводим информацию о том что долшли до решения
    int iters = solve_MSG_IC(cfg, A, S, f, x, extra);      // решаем
    cout << "Done in " << iters << " iterations." << endl; // выводим количество иттераций за которое решение сошлось

    cout << endl << "Temperature Field:" << endl; // вспомогательный вывод
    cout.precision(2);                            // ставим 2 знака после запятой для красоты
    cout << fixed;                                // фиксированный формат

    // выводим массив решения x как двумерную сетку.
    for (int j = NodesY - 1; j >= 0; j--) {  // строки (Y)
        for (int i = 0; i < NodesX; i++) {   // столбцы (X)
            int node_idx = j * NodesX + i;   // формула пересчета координат в индекс
            double val = x[node_idx];        // вывод
            if (abs(val) < 1e-10) val = 0.0; // убираем "минус ноль" (-0.00 иногда вылезает из-за погрешности float)

            cout << setw(8) << val;          // выравниваем колонки по ширине 8 символов
        }
        cout << endl;                        // переход на новую строку сетки
    }

    //специальный вывод для теста 4 (сравниваем значение в самом центре области с точным аналитическим решением которое равно 1.0)
    if (cfg.TestId == 4) {
        int center_node = NodesY * (NodesX / 2) + NodesY / 2;                 // индекс центрального узла
        cout << endl << "Error in center: " << scientific << setprecision(15) // выводим |полученное - точное|
             << abs(x[center_node]) - 1.0 << endl;
    }

    system("pause"); // чтобы консоль не закрылась сразу
    return 0;
}