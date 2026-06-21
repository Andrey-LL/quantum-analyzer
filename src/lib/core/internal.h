// Внутренние структуры и утилиты для gaussian_parser.cpp, matrix_ops.cpp, quantum_analysis.cpp
#ifndef GAUSSIAN_INTERNAL_H
#define GAUSSIAN_INTERNAL_H

#include <string>
#include <vector>
#include <tuple>
#include <map>
#include <memory>
#include <cmath>
#include <boost/dynamic_bitset.hpp>
#include <Eigen/Dense>

// Предварительные объявления для разрыва циклических зависимостей.
struct GaussianFileImpl;
struct MatrixHandleImpl;
struct GroupHandleImpl;
struct MOData;

// ============================================================================
// Константы и настройки
// ============================================================================

// Точность для проверки симметрии матриц
constexpr double SYMMETRY_TOLERANCE = 1e-10;

// ============================================================================
// Внутренние структуры данных (реализация для непрозрачных хэндлов)
// ============================================================================

/**
 * Данные молекулярных орбиталей
 */
struct MOData {
    Eigen::MatrixXd coefficients;          // nbasis x nmo.
    Eigen::VectorXd eigenvalues;           // nmo.
    std::vector<std::string> symmetries;   // Симметрии МО.
    std::vector<int> ao_to_atom;           // AO -> атом.
    std::vector<std::string> ao_labels;    // Метки AO.
    MOData() {}
};

/**
 * Внутренняя реализация хэндла файла Gaussian
 */
struct GaussianFileImpl {
    std::string filename;
    std::string content;
    int nbasis = 0;
    int natoms = 0;
    int alpha_electrons = 0;
    int beta_electrons = 0;
    std::string basis_set;
    std::string method;
    double nuclear_repulsion = 0.0;
    // Геометрия хранится как (атомный номер, символ, x, y, z).
    std::vector<std::tuple<int, std::string, double, double, double>> geometry;
    bool is_open = false;
    std::map<std::string, Eigen::MatrixXd> matrix_cache;
    // Кэш данных молекулярных орбиталей.
    bool mo_cached = false;
    MOData mo;
};

/**
 * Внутренняя реализация хэндла матрицы
 */
struct MatrixHandleImpl {
    Eigen::MatrixXd data;
    std::string type;
    int rows, cols;
    // Статус симметрии: 0 = неизвестно, 1 = симметрична, 2 = несимметрична.
    int is_symmetric;
    double trace_val;
    double condition_number;

    MatrixHandleImpl() : rows(0), cols(0), is_symmetric(0), trace_val(NAN), condition_number(-1.0) {}

    MatrixHandleImpl(const Eigen::MatrixXd& matrix, const std::string& matrix_type)
            : data(matrix), type(matrix_type), rows(matrix.rows()), cols(matrix.cols()),
              is_symmetric(0), trace_val(NAN), condition_number(-1.0)
    {
        // Только для квадратных матриц проверяем симметрию и вычисляем trace.
        if (rows == cols && rows > 0) {
            trace_val = data.trace();
            refresh_symmetry();

            // Число обусловленности вычисляется только для симметричных overlap матриц.
            if (is_symmetric == 1 && type == "overlap") {
                Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(data);
                auto eigenvalues = solver.eigenvalues();
                if (eigenvalues.size() > 0 && eigenvalues.minCoeff() > 1e-12) {
                    condition_number = eigenvalues.maxCoeff() / eigenvalues.minCoeff();
                } else {
                    condition_number = -1.0;
                }
            } else {
                condition_number = -1.0;
            }
        } else {
            // Неквадратная матрица: trace неприменим, симметрия не определена.
            is_symmetric = 0;
            trace_val = NAN;
            condition_number = -1.0;
        }
    }

    bool refresh_symmetry() {
        if (is_symmetric == 1) return true;
        if (is_symmetric == 2) return false;
        if (rows != cols || rows <= 0) {
            is_symmetric = 0;
            return false;
        }

        bool sym = true;
        for (int i = 0; i < rows && sym; ++i) {
            for (int j = i + 1; j < cols; ++j) {
                double a = data(i, j);
                double b = data(j, i);
                if (!std::isfinite(a) || !std::isfinite(b) || std::fabs(a - b) > SYMMETRY_TOLERANCE) {
                    sym = false;
                    break;
                }
            }
        }
        is_symmetric = sym ? 1 : 2;
        return sym;
    }
};

/**
 * Внутренняя реализация хэндла группы (битовая маска)
 */
struct GroupHandleImpl {
    int nbasis = 0;
    boost::dynamic_bitset<> bits;
    explicit GroupHandleImpl(int n)
            : nbasis(n), bits((n > 0) ? static_cast<size_t>(n) : 0) {}
};

// ============================================================================
// Внутренние утилиты
// ============================================================================

// Настраиваемые ограничения для анализа блоков плотности.
#ifndef QA_MAX_ATOMS
#define QA_MAX_ATOMS 256
#endif
#ifndef QA_MAX_ATOM_PAIRS
#define QA_MAX_ATOM_PAIRS (QA_MAX_ATOMS * (QA_MAX_ATOMS - 1) / 2)
#endif

// Внутренние структуры анализа блоков плотности; соответствуют typedef из api.h.
// Данные матриц выделяются динамически, без фиксированных массивов.
struct DiagonalBlockInfo {
    int atom;              // Индекс атома.
    int size;              // n_A × n_A.
    double trace;          // Tr(P_AA).
    double frob_norm;      // ||P_AA||_F.
    double* data;          // Матрица в row-major формате, длина size*size.
};

struct OffDiagonalBlockInfo {
    int atom_a, atom_b;    // Пара атомов.
    int size_a, size_b;    // n_A, n_B.
    double frob_norm;      // ||P_AB||_F.
    double* data;          // Матрица в row-major формате, длина size_a*size_b.
};

struct DensityBlockAnalysis {
    int natoms;
    int nbasis;
    std::vector<DiagonalBlockInfo> diagonal_blocks;
    std::vector<OffDiagonalBlockInfo> offdiag_blocks;
    int n_diagonal;
    int n_offdiagonal;
    
    DensityBlockAnalysis() : natoms(0), nbasis(0), n_diagonal(0), n_offdiagonal(0) {}
    
    void reserve(int expected_atoms, int expected_pairs) {
        diagonal_blocks.reserve(expected_atoms);
        offdiag_blocks.reserve(expected_pairs);
    }
    
    void resize(int atoms, int pairs) {
        diagonal_blocks.resize(atoms);
        offdiag_blocks.resize(pairs);
    }
};

/**
 * Вычисляет индекс блока столбцов для треугольной матрицы
 * @param n Размер матрицы
 * @param k Номер строки в выводе Gaussian
 * @return Номер блока столбцов (j)
 */
int n_j(int n, int k);

/**
 * Очищает D-нотацию в экспоненциальных числах (D+03 -> E+03)
 * @param text Входная строка
 * @return Очищенная строка
 */
std::string clean_d_notation(const std::string& text);

/**
 * Извлекает числовые значения из строки
 * @param line Строка с числами
 * @return Вектор чисел
 */
std::vector<double> extract_numbers(const std::string& line);

// ============================================================================
// Функции парсинга
// ============================================================================

/**
 * Извлекает метаданные из содержимого файла
 * @param file Указатель на структуру файла
 * @return true при успехе, false при ошибке
 */
bool extract_metadata(GaussianFileImpl* file);

/**
 * Извлекает матрицу заданного типа
 * @param file Указатель на структуру файла
 * @param matrix_type Тип матрицы
 * @return Матрица Eigen
 */
Eigen::MatrixXd extract_matrix(const GaussianFileImpl* file, const std::string& matrix_type);

// ============================================================================
// Константы и настройки
// ============================================================================

// Словарь маркеров для различных типов матриц
extern const std::map<std::string, std::pair<std::string, bool>> MATRIX_MARKERS;

// Словарь атомных номеров -> символы элементов
extern const std::map<int, std::string> ATOMIC_SYMBOLS;

// Максимальное количество итераций в функции n_j
constexpr int MAX_NJ_ITERATIONS = 1000;

// ============================================================================
// Вспомогательные функции преобразования дескрипторов (для внутреннего использования)
// ============================================================================
/*
static inline MatrixHandleImpl* as_matrix(MatrixHandle h) {
    return reinterpret_cast<MatrixHandleImpl*>(h);
}

static inline const MatrixHandleImpl* as_matrix_const(MatrixHandle h) {
    return reinterpret_cast<const MatrixHandleImpl*>(h);
}

static inline GroupHandleImpl* as_group(GroupHandle h) {
    return reinterpret_cast<GroupHandleImpl*>(h);
}

static inline bool is_finite(double x) {
    return std::isfinite(x);
}

// Итерация по установленным битам без сканирования всех битов
template <typename Func>
static inline void for_each_set_bit(const boost::dynamic_bitset<>& b, Func f) {
    for (auto i = b.find_first(); i != boost::dynamic_bitset<>::npos; i = b.find_next(i)) {
        f(static_cast<int>(i));
    }
}
*/
// ============================================================================
// Валидация
// ============================================================================

/**
 * Проверяет корректность извлеченной overlap матрицы
 * @param matrix Overlap матрица
 * @param expected_size Ожидаемый размер
 * @return true если корректная
 */
bool validate_overlap_matrix_internal(const Eigen::MatrixXd& matrix, int expected_size);

/**
 * Проверяет корректность извлеченной density матрицы
 * @param density_matrix Density матрица
 * @param overlap_matrix Overlap матрица
 * @param expected_electrons Ожидаемое количество электронов
 * @return true если корректная
 */
bool validate_density_matrix_internal(const Eigen::MatrixXd& density_matrix,
                                      const Eigen::MatrixXd& overlap_matrix,
                                      int expected_electrons);

// ============================================================================
// Отладка и логирование (только в debug режиме)
// ============================================================================

#ifdef GAUSSIAN_DEBUG
/**
 * Выводит статистику матрицы для отладки
 * @param matrix Матрица
 * @param name Название матрицы
 */
void debug_print_matrix_stats(const Eigen::MatrixXd& matrix, const std::string& name);

/**
 * Выводит первые n×n элементов матрицы
 * @param matrix Матрица
 * @param n Количество строк/столбцов для вывода
 * @param name Название матрицы
 */
void debug_print_matrix_corner(const Eigen::MatrixXd& matrix, int n, const std::string& name);
#endif

#endif
