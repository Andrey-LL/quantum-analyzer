// Математические примитивы для работы с плотными матрицами.
// Реализует линейную алгебру, поэлементные преобразования и блочные редукции.
#include "api.h"
#include "internal.h"
#include <boost/dynamic_bitset/dynamic_bitset.hpp>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <cmath>
#include <limits>
#include <memory>
#include <algorithm>
#include <stdexcept>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::string;

// Предварительное объявление LAPACK dsyev.
extern "C" {
void dsyev_(char* jobz, char* uplo, int* n, double* a, int* lda, double* w,
            double* work, int* lwork, int* info);
}

// Предварительные объявления для внутренних структур.
struct GaussianFileImpl;
struct MatrixHandleImpl;
struct GroupHandleImpl;
struct MOData;

// -----------------------------------------------------------------------------
// Вспомогательные функции преобразования дескрипторов
// -----------------------------------------------------------------------------
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

static inline bool check_same_size(const MatrixHandleImpl* A, const MatrixHandleImpl* B) {
    return A && B && (A->rows == B->rows) && (A->cols == B->cols);
}

static inline bool require_square(const MatrixHandleImpl* A) {
    return A && (A->rows == A->cols);
}

static inline bool require_symmetric(MatrixHandleImpl* A) {
    return A && A->refresh_symmetry();
}

// Итерация по установленным битам без полного сканирования битовой маски.
template <typename Func>
static inline void for_each_set_bit(const boost::dynamic_bitset<>& b, Func f) {
    for (auto i = b.find_first(); i != boost::dynamic_bitset<>::npos; i = b.find_next(i)) {
        f(static_cast<int>(i));
    }
}

// Проверка согласованности размерностей для умножения
static inline bool check_multiplication_size(const MatrixHandleImpl* A, const MatrixHandleImpl* B) {
    return A && B && (A->cols == B->rows);
}

// -----------------------------------------------------------------------------
// C API реализации
// -----------------------------------------------------------------------------
extern "C" {

// ===== Базовые операции с матрицами =====

/**
 * Создание матрицы из row-major массива double.
 */
QUANTUM_ANALYZER_API MatrixHandle matrix_create(int rows, int cols, const double* data) {
    if (rows <= 0 || cols <= 0 || !data) return nullptr;
    try {
        MatrixXd M(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                double v = data[i * cols + j];
                if (!is_finite(v)) return nullptr;
                M(i, j) = v;
            }
        }
        auto impl = std::make_unique<MatrixHandleImpl>(M, "manual");
        return reinterpret_cast<MatrixHandle>(impl.release());
    } catch (...) {
        return nullptr;
    }
}

/**
 * Преобразование матрицы из AO-представления в MO-представление: C^T * M * C
 * Требования:
 *   - Количество строк в матрице коэффициентов MO (C) должно совпадать с размером матрицы M
 *   - Матрица коэффициентов должна иметь корректные размеры (nbasis × nmo)
 * @param matrix_ao Исходная матрица в AO-базисе
 * @param mo_coeff Матрица коэффициентов молекулярных орбиталей (C)
 * @return Новая матрица в MO-базисе или nullptr при ошибке
 */
QUANTUM_ANALYZER_API MatrixHandle matrix_ao_to_mo(MatrixHandle matrix_ao, MatrixHandle mo_coeff) {
    if (!matrix_ao || !mo_coeff) return nullptr;
    try {
        auto impl_ao = as_matrix(matrix_ao);
        auto impl_mo = as_matrix(mo_coeff);
        // Перед умножением проверяем согласованность размерностей.
        const MatrixXd &C = impl_mo->data;
        const MatrixXd &A = impl_ao->data;
        if (C.cols() <= 0 || C.rows() != A.rows()) return nullptr;
        // Преобразование AO -> MO: C^T * A * C.
        MatrixXd Ct = C.transpose();
        MatrixXd temp = Ct * A;
        MatrixXd result = temp * C;
        auto result_impl = std::make_unique<MatrixHandleImpl>(result, impl_ao->type + "_MO");
        return reinterpret_cast<MatrixHandle>(result_impl.release());
    } catch (...) { return nullptr; }
}

/**
 * Получение размеров матрицы
 * @param matrix Хэндл матрицы
 * @param rows Указатель для записи количества строк (может быть nullptr)
 * @param cols Указатель для записи количества столбцов (может быть nullptr)
 * @return 0 при успехе, -1 при ошибке (некорректный хэндл или оба указателя nullptr)
 */
QUANTUM_ANALYZER_API int matrix_get_size(MatrixHandle matrix, int *rows, int *cols) {
    if (!matrix) return -1;
    if (!rows && !cols) return -1; // Требуется хотя бы один выходной указатель.
    auto impl = as_matrix(matrix);
    if (rows) *rows = impl->rows;
    if (cols) *cols = impl->cols;
    return 0;
}

/**
 * Получение элемента матрицы по индексам
 * @param matrix Хэндл матрицы
 * @param row Индекс строки
 * @param col Индекс столбца
 * @return Значение элемента или NAN при ошибке (некорректный хэндл или индексы)
 */
QUANTUM_ANALYZER_API double matrix_get_element(MatrixHandle matrix, int row, int col) {
    if (!matrix) return std::numeric_limits<double>::quiet_NaN();
    auto impl = as_matrix(matrix);
    if (row < 0 || row >= impl->rows || col < 0 || col >= impl->cols)
        return std::numeric_limits<double>::quiet_NaN();
    return impl->data(row, col);
}

/**
 * След матрицы: Tr(A) = sum_i A_ii
 * Требования: квадратная матрица (проверяется при создании хэндла)
 * @param matrix Хэндл матрицы
 * @return Значение следа или NAN при ошибке
 */
QUANTUM_ANALYZER_API double matrix_trace(MatrixHandle matrix) {
    if (!matrix) return std::numeric_limits<double>::quiet_NaN();
    auto impl = as_matrix(matrix);
    return impl->trace_val;
}

/**
 * Проверка симметричности матрицы
 * @param matrix Хэндл матрицы
 * @return 0 = unknown, 1 = yes, 2 = no, -1 при ошибке хэндла
 */
QUANTUM_ANALYZER_API int matrix_is_symmetric(MatrixHandle matrix) {
    if (!matrix) return -1;
    auto impl = as_matrix(matrix);
    impl->refresh_symmetry();
    return impl->is_symmetric;
}

/**
 * Получение числа обусловленности матрицы.
 * Если значение ещё не вычислялось (impl->condition_number < 0), вычисляет и кэширует.
 * @param matrix Хэндл матрицы
 * @return Число обусловленности; NAN при ошибке, INF при вырожденной матрице
 */
QUANTUM_ANALYZER_API double matrix_condition_number(MatrixHandle matrix) {
    if (!matrix) return std::numeric_limits<double>::quiet_NaN();

    auto impl = as_matrix(matrix);

    // Невыраженное отрицательным числом condition_number означает валидное закэшированное значение.
    if (impl->condition_number >= 0.0) {
        return impl->condition_number;
    }

    const auto& A = impl->data;
    if (A.rows() == 0 || A.cols() == 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    double cond;
    const double eps = 1e-15;

    // Для симметричных матриц используем разложение по собственным значениям.
    if (require_symmetric(impl)) {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(A);
        if (solver.info() != Eigen::Success) {
            return std::numeric_limits<double>::quiet_NaN();
        }

        const auto& eig = solver.eigenvalues();
        if (eig.size() == 0) {
            return std::numeric_limits<double>::quiet_NaN();
        }

        const double emax = std::fabs(eig.maxCoeff());
        const double emin = std::fabs(eig.minCoeff());

        if (!std::isfinite(emax) || !std::isfinite(emin)) {
            return std::numeric_limits<double>::quiet_NaN();
        }

        if (emin <= eps) {
            cond = std::numeric_limits<double>::infinity();
        } else {
            cond = emax / emin;
        }
    } else {
        // Для несимметричных/прямоугольных матриц используем SVD
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
        const auto& s = svd.singularValues();

        if (s.size() == 0) {
            return std::numeric_limits<double>::quiet_NaN();
        }

        // Сингулярные значения отсортированы по убыванию: s(0) = максимум, s(end) = минимум.
        const double smax = s(0);
        const double smin = s(s.size() - 1);

        if (!std::isfinite(smax) || !std::isfinite(smin)) {
            return std::numeric_limits<double>::quiet_NaN();
        }

        // Сингулярные значения по определению неотрицательны.
        if (smin <= eps) {
            cond = std::numeric_limits<double>::infinity();
        } else {
            cond = smax / smin;
        }
    }

    impl->condition_number = cond;
    return cond;
}

/**
 * Вычисление собственных значений симметричной матрицы
 * Требования: матрица должна быть симметричной (проверяется флагом)
 * @param matrix Хэндл симметричной матрицы
 * @param eigenvalues Буфер для записи собственных значений (должен быть выделен вызывающей стороной)
 * @return 0 при успехе, -1 при ошибке
 */
QUANTUM_ANALYZER_API int matrix_eigenvalues(MatrixHandle matrix, double *eigenvalues) {
    if (!matrix || !eigenvalues) return -1;
    try {
        auto impl = as_matrix(matrix);
        if (!require_symmetric(impl)) return -1;
        Eigen::SelfAdjointEigenSolver<MatrixXd> solver(impl->data);
        VectorXd eig_vals = solver.eigenvalues();
        for (int i = 0; i < eig_vals.size(); ++i) eigenvalues[i] = eig_vals(i);
        return 0;
    } catch (...) { return -1; }
}

/**
 * Диагонализация симметричной матрицы: возвращает собственные значения И векторы
 * A = U · Λ · U^T, где Λ — диагональная матрица eigenvalues, U — матрица eigenvectors
 * @param matrix матрица
 * @param eigenvalues Буфер для собственных значений [n]
 * @param eigenvectors Выходная матрица [n×n], столбцы = собственные векторы
 * @return 0 при успехе, -1 при ошибке
 */
QUANTUM_ANALYZER_API int matrix_eigensystem(MatrixHandle matrix, double* eigenvalues, MatrixHandle* eigenvectors) {
    if (!matrix || !eigenvalues || !eigenvectors) return -1;
    try {
        auto impl = as_matrix(matrix);
        
        MatrixXd eig_vecs;
        VectorXd eig_vals;

        // Если матрица симметрична, используем быстрый решатель.
        if (require_symmetric(impl)) {
            Eigen::SelfAdjointEigenSolver<MatrixXd> solver(impl->data);
            if (solver.info() != Eigen::Success) return -1;
            eig_vals = solver.eigenvalues();
            eig_vecs = solver.eigenvectors();
        } else {
            // Для несимметричных матриц (например, блоков PS)
            Eigen::EigenSolver<MatrixXd> solver(impl->data);
            if (solver.info() != Eigen::Success) return -1;
            // Берем вещественную часть (мнимая для PS_AA должна быть ~0)
            eig_vals = solver.eigenvalues().real();
            eig_vecs = solver.eigenvectors().real();
        }

        // Записываем собственные значения
        for (int i = 0; i < eig_vals.size(); ++i) {
            eigenvalues[i] = eig_vals(i);
        }

        // Создаём матрицу собственных векторов
        auto evecs_impl = std::make_unique<MatrixHandleImpl>();
        evecs_impl->data = eig_vecs;
        evecs_impl->type = "eigenvectors";
        evecs_impl->rows = eig_vecs.rows();
        evecs_impl->cols = eig_vecs.cols();
        evecs_impl->is_symmetric = 0;  // не симметрична
        evecs_impl->trace_val = eig_vecs.trace();
        evecs_impl->condition_number = -1.0;

        *eigenvectors = reinterpret_cast<MatrixHandle>(evecs_impl.release());
        return 0;
    } catch (...) { return -1; }
}

/**
 * SVD разложение: A = U · Σ · V^T
 * Использует Eigen::JacobiSVD для полной декомпозиции
 * @param matrix Матрица [m×n]
 * @param U Выходная матрица [m×m] — левые сингулярные векторы
 * @param S Выходная диагональная матрица [min(m,n)×min(m,n)] — сингулярные числа
 * @param Vt Выходная матрица [n×n] — транспонированные правые векторы
 * @return 0 при успехе, -1 при ошибке
 */
QUANTUM_ANALYZER_API int matrix_svd(MatrixHandle matrix, MatrixHandle* U, MatrixHandle* S, MatrixHandle* Vt) {
    if (!matrix || !U || !S || !Vt) return -1;
    try {
        auto impl = as_matrix(matrix);
        const MatrixXd& A = impl->data;

        int m = A.rows();
        int n = A.cols();
        int k = std::min(m, n);

        // Полное SVD разложение.
        Eigen::JacobiSVD<MatrixXd> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
        if (svd.info() != Eigen::Success) return -1;

        MatrixXd U_mat = svd.matrixU();
        MatrixXd V_mat = svd.matrixV();
        VectorXd singular_values = svd.singularValues();

        // Создаём U [m×m].
        auto U_impl = std::make_unique<MatrixHandleImpl>();
        U_impl->data = U_mat;
        U_impl->type = "svd_U";
        U_impl->rows = m;
        U_impl->cols = m;
        U_impl->is_symmetric = 0;
        U_impl->trace_val = U_mat.trace();
        U_impl->condition_number = -1.0;
        *U = reinterpret_cast<MatrixHandle>(U_impl.release());

        // Создаём диагональную матрицу Σ [k×k].
        MatrixXd S_mat = MatrixXd::Zero(k, k);
        for (int i = 0; i < k; ++i) {
            S_mat(i, i) = singular_values(i);
        }
        auto S_impl = std::make_unique<MatrixHandleImpl>();
        S_impl->data = S_mat;
        S_impl->type = "svd_S";
        S_impl->rows = k;
        S_impl->cols = k;
        S_impl->is_symmetric = 1;  // диагональная
        S_impl->trace_val = S_mat.trace();
        S_impl->condition_number = (singular_values(k-1) > 1e-15) ? singular_values(0) / singular_values(k-1) : std::numeric_limits<double>::infinity();
        *S = reinterpret_cast<MatrixHandle>(S_impl.release());

        // Создаём V^T [n×n].
        auto Vt_impl = std::make_unique<MatrixHandleImpl>();
        Vt_impl->data = V_mat.transpose();
        Vt_impl->type = "svd_Vt";
        Vt_impl->rows = n;
        Vt_impl->cols = n;
        Vt_impl->is_symmetric = 0;
        Vt_impl->trace_val = Vt_impl->data.trace();
        Vt_impl->condition_number = -1.0;
        *Vt = reinterpret_cast<MatrixHandle>(Vt_impl.release());

        return 0;
    } catch (...) { return -1; }
}

/**
 * Возведение симметричной матрицы в произвольную степень через спектральное разложение
 * Требования: матрица должна быть симметричной
 * @param matrix Хэндл симметричной матрицы
 * @param exponent Степень (может быть дробной или отрицательной)
 * @return Новая матрица A^exponent или nullptr при ошибке
 */
QUANTUM_ANALYZER_API MatrixHandle matrix_power(MatrixHandle matrix, double exponent) {
    if (!matrix) return nullptr;
    try {
        auto impl = as_matrix(matrix);
        if (!require_symmetric(impl)) return nullptr;
        Eigen::SelfAdjointEigenSolver<MatrixXd> solver(impl->data);
        VectorXd eigenvalues = solver.eigenvalues();
        MatrixXd eigenvectors = solver.eigenvectors();
        for (int i = 0; i < eigenvalues.size(); ++i)
            eigenvalues(i) = std::pow(std::max(eigenvalues(i), 1e-12), exponent);
        MatrixXd result = eigenvectors * eigenvalues.asDiagonal() * eigenvectors.transpose();
        auto result_impl = std::make_unique<MatrixHandleImpl>(result, impl->type + "^" + std::to_string(exponent));
        return reinterpret_cast<MatrixHandle>(result_impl.release());
    } catch (...) { return nullptr; }
}

/**
 * Матричное умножение: C = A·B
 * Требования: A.cols == B.rows
 * @param a Хэндл левой матрицы
 * @param b Хэндл правой матрицы
 * @return Новая матрица-произведение или nullptr при ошибке
 */
QUANTUM_ANALYZER_API MatrixHandle matrix_multiply(MatrixHandle a, MatrixHandle b) {
    if (!a || !b) return nullptr;
    try {
        auto impl_a = as_matrix(a);
        auto impl_b = as_matrix(b);
        if (impl_a->cols != impl_b->rows) return nullptr;
        MatrixXd result = impl_a->data * impl_b->data;
        auto result_impl = std::make_unique<MatrixHandleImpl>(result, impl_a->type + "*" + impl_b->type);
        return reinterpret_cast<MatrixHandle>(result_impl.release());
    } catch (...) { return nullptr; }
}

/**
 * Освобождение памяти, занятой матрицей
 * @param matrix Хэндл матрицы для освобождения
 */
QUANTUM_ANALYZER_API void matrix_free(MatrixHandle matrix) {
    if (matrix) delete as_matrix(matrix);
}

// ===== Базовая линейная алгебра =====

/**
 * Линейная комбинация матриц: C = α·A + β·B
 * Требования: A и B одинакового размера
 * Результат симметричен, если обе входные матрицы симметричны
 */
QUANTUM_ANALYZER_API MatrixHandle matrix_add(MatrixHandle a, double alpha, MatrixHandle b, double beta) {
    if (!a || !b || !is_finite(alpha) || !is_finite(beta)) return nullptr;
    try {
        auto A = as_matrix(a);
        auto B = as_matrix(b);
        if (!check_same_size(A, B)) return nullptr;

        MatrixXd C = alpha * A->data + beta * B->data;

        auto out = std::make_unique<MatrixHandleImpl>();
        out->data = std::move(C);
        out->type = "linear_combination";
        out->rows = A->rows;
        out->cols = A->cols;
        // Линейная комбинация симметричных матриц всегда симметрична
        out->is_symmetric = (A->is_symmetric == 1 && B->is_symmetric == 1) ? 1 : 0;
        out->trace_val = (out->rows == out->cols) ? out->data.trace() : std::numeric_limits<double>::quiet_NaN();
        out->condition_number = -1.0;
        return reinterpret_cast<MatrixHandle>(out.release());
    } catch (...) {
        return nullptr;
    }
}

/**
 * Масштабирование матрицы: B = α·A
 */
QUANTUM_ANALYZER_API MatrixHandle matrix_scale(MatrixHandle a, double alpha) {
    if (!a || !is_finite(alpha)) return nullptr;
    try {
        auto A = as_matrix(a);

        MatrixXd B = alpha * A->data;

        auto out = std::make_unique<MatrixHandleImpl>();
        out->data = std::move(B);
        out->type = "scaled";
        out->rows = A->rows;
        out->cols = A->cols;
        out->is_symmetric = A->is_symmetric;  // пробрасываем статус (0/1/2)
        out->trace_val = (out->rows == out->cols) ? out->data.trace() : std::numeric_limits<double>::quiet_NaN();
        out->condition_number = -1.0;
        return reinterpret_cast<MatrixHandle>(out.release());
    } catch (...) {
        return nullptr;
    }
}

/**
 * Транспонирование матрицы: B = A^T
 */
QUANTUM_ANALYZER_API MatrixHandle matrix_transpose(MatrixHandle a) {
    if (!a) return nullptr;
    try {
        auto A = as_matrix(a);

        MatrixXd B = A->data.transpose();

        auto out = std::make_unique<MatrixHandleImpl>();
        out->data = std::move(B);
        out->type = "transpose";
        out->rows = A->cols;
        out->cols = A->rows;
        out->is_symmetric = A->is_symmetric;  // транспонирование сохраняет симметрию
        out->trace_val = (out->rows == out->cols) ? out->data.trace() : std::numeric_limits<double>::quiet_NaN();
        out->condition_number = A->condition_number;
        return reinterpret_cast<MatrixHandle>(out.release());
    } catch (...) {
        return nullptr;
    }
}

/**
 * Тройное произведение: C = A·B·A^T
 * Симметризация применяется если B симметрична (для устранения численного шума).
 * @param a Левая матрица A (размер m×n)
 * @param b Правая матрица B (должна быть квадратной n×n)
 * @return Новая матрица C размером m×m
 */
QUANTUM_ANALYZER_API MatrixHandle matrix_triple_product_symm(MatrixHandle a, MatrixHandle b) {
    if (!a || !b) return nullptr;
    try {
        auto A = as_matrix(a);
        auto B = as_matrix(b);
        if (!require_square(B) || A->cols != B->rows) return nullptr;

        // Вычисляем C = A * B * A^T.
        MatrixXd C = A->data * B->data * A->data.transpose();

        // Симметризация применяется только при симметричной B, чтобы убрать численный шум.
        int result_is_symmetric = (B->is_symmetric == 1) ? 1 : 0;
        if (result_is_symmetric == 1) {
            C = (C + C.transpose()) * 0.5;
        }

        auto out = std::make_unique<MatrixHandleImpl>();
        out->data = std::move(C);
        out->type = "triple_product_symm";
        out->rows = A->rows;
        out->cols = A->rows;
        out->is_symmetric = result_is_symmetric;
        out->trace_val = C.trace();
        out->condition_number = -1.0;
        return reinterpret_cast<MatrixHandle>(out.release());
    } catch (...) {
        return nullptr;
    }
}

// ===== Поэлементные операции =====

/**
 * Поэлементное произведение (Адамар): C_ij = A_ij * B_ij
 * Требования: одинаковые размеры
 */
QUANTUM_ANALYZER_API MatrixHandle matrix_hadamard(MatrixHandle a, MatrixHandle b) {
    if (!a || !b) return nullptr;
    try {
        auto A = as_matrix(a);
        auto B = as_matrix(b);
        if (!check_same_size(A, B)) return nullptr;

        MatrixXd C = A->data.cwiseProduct(B->data);

        auto out = std::make_unique<MatrixHandleImpl>();
        out->data = std::move(C);
        out->type = "hadamard";
        out->rows = A->rows;
        out->cols = A->cols;
        out->is_symmetric = (A->is_symmetric == 1 && B->is_symmetric == 1) ? 1 : 0;
        out->trace_val = (out->rows == out->cols) ? out->data.trace() : std::numeric_limits<double>::quiet_NaN();
        out->condition_number = -1.0;
        return reinterpret_cast<MatrixHandle>(out.release());
    } catch (...) {
        return nullptr;
    }
}

/**
 * Поэлементное возведение в степень: C_ij = A_ij^p
 * Правила домена:
 *   - Для целых степеней: разрешены любые конечные значения (включая отрицательные)
 *   - Для нецелых степеней: требуется A_ij >= 0
 *   - Для отрицательных степеней: требуется A_ij != 0
 */
QUANTUM_ANALYZER_API MatrixHandle matrix_cwise_pow(MatrixHandle a, double exponent) {
    if (!a || !is_finite(exponent)) return nullptr;
    try {
        auto A = as_matrix(a);

        // Проверка целочисленности степени
        double rintp = std::round(exponent);
        bool is_integer = std::fabs(exponent - rintp) < 1e-12;

        // Проверка допустимости элементов
        for (int i = 0; i < A->rows; ++i) {
            for (int j = 0; j < A->cols; ++j) {
                double v = A->data(i, j);
                if (!is_finite(v)) return nullptr;
                if (exponent < 0.0 && std::fabs(v) < 1e-15) return nullptr; // Деление на ноль.
                if (!is_integer && v < 0.0) return nullptr; // Результат был бы комплексным.
            }
        }

        MatrixXd C = A->data.array().pow(exponent).matrix();

        auto out = std::make_unique<MatrixHandleImpl>();
        out->data = std::move(C);
        out->type = "cwise_pow";
        out->rows = A->rows;
        out->cols = A->cols;
        out->is_symmetric = A->is_symmetric;  // Статус симметрии сохраняется.
        out->trace_val = (out->rows == out->cols) ? out->data.trace() : std::numeric_limits<double>::quiet_NaN();
        out->condition_number = -1.0;
        return reinterpret_cast<MatrixHandle>(out.release());
    } catch (...) {
        return nullptr;
    }
}

/**
 * Пороговая фильтрация: C_ij = (|A_ij| >= theta) ? A_ij : 0.0
 */
QUANTUM_ANALYZER_API MatrixHandle matrix_threshold(MatrixHandle a, double theta) {
    if (!a || !is_finite(theta) || theta < 0.0) return nullptr;
    try {
        auto A = as_matrix(a);
        double t = std::fabs(theta);

        MatrixXd C = A->data.unaryExpr([t](double x) {
            return (std::fabs(x) >= t) ? x : 0.0;
        });

        auto out = std::make_unique<MatrixHandleImpl>();
        out->data = std::move(C);
        out->type = "threshold";
        out->rows = A->rows;
        out->cols = A->cols;
        out->is_symmetric = A->is_symmetric;  // пробрасываем статус
        out->trace_val = (out->rows == out->cols) ? out->data.trace() : std::numeric_limits<double>::quiet_NaN();
        out->condition_number = -1.0;
        return reinterpret_cast<MatrixHandle>(out.release());
    } catch (...) {
        return nullptr;
    }
}

/**
 * Ограничение значений (clamp): C_ij = min(max(A_ij, lo), hi)
 */
QUANTUM_ANALYZER_API MatrixHandle matrix_clamp(MatrixHandle a, double lo, double hi) {
    if (!a || !is_finite(lo) || !is_finite(hi) || lo > hi) return nullptr;
    try {
        auto A = as_matrix(a);

        MatrixXd C = A->data.unaryExpr([lo, hi](double x) {
            return std::min(std::max(x, lo), hi);
        });

        auto out = std::make_unique<MatrixHandleImpl>();
        out->data = std::move(C);
        out->type = "clamp";
        out->rows = A->rows;
        out->cols = A->cols;
        out->is_symmetric = A->is_symmetric;  // пробрасываем статус
        out->trace_val = (out->rows == out->cols) ? out->data.trace() : std::numeric_limits<double>::quiet_NaN();
        out->condition_number = -1.0;
        return reinterpret_cast<MatrixHandle>(out.release());
    } catch (...) {
        return nullptr;
    }
}

// ===== Редукции =====

/**
 * Норма Фробениуса: ||A||_F = sqrt(sum_ij A_ij^2)
 */
QUANTUM_ANALYZER_API double matrix_norm_fro(MatrixHandle a) {
    if (!a) return std::numeric_limits<double>::quiet_NaN();
    try {
        auto A = as_matrix(a);
        double norm = A->data.norm();
        return is_finite(norm) ? norm : std::numeric_limits<double>::quiet_NaN();
    } catch (...) {
        return std::numeric_limits<double>::quiet_NaN();
    }
}

/**
 * Максимальный по модулю элемент матрицы
 */
QUANTUM_ANALYZER_API double matrix_max_abs(MatrixHandle a) {
    if (!a) return std::numeric_limits<double>::quiet_NaN();
    try {
        auto A = as_matrix(a);
        double max_val = 0.0;
        for (int i = 0; i < A->rows; ++i) {
            for (int j = 0; j < A->cols; ++j) {
                double v = std::fabs(A->data(i, j));
                if (v > max_val) max_val = v;
            }
        }
        return is_finite(max_val) ? max_val : std::numeric_limits<double>::quiet_NaN();
    } catch (...) {
        return std::numeric_limits<double>::quiet_NaN();
    }
}

/**
 * Минимальный по модулю ненулевой элемент матрицы
 */
QUANTUM_ANALYZER_API double matrix_min_abs_nonzero(MatrixHandle a) {
    if (!a) return std::numeric_limits<double>::quiet_NaN();
    try {
        auto A = as_matrix(a);
        double min_val = std::numeric_limits<double>::max();
        bool found = false;
        for (int i = 0; i < A->rows; ++i) {
            for (int j = 0; j < A->cols; ++j) {
                double v = std::fabs(A->data(i, j));
                if (v > 1e-15 && v < min_val) {
                    min_val = v;
                    found = true;
                }
            }
        }
        return (found && is_finite(min_val)) ? min_val : std::numeric_limits<double>::quiet_NaN();
    } catch (...) {
        return std::numeric_limits<double>::quiet_NaN();
    }
}

/**
 * Извлечение диагонали как вектора
 * Возвращает новый MatrixHandle размером N×1
 */
QUANTUM_ANALYZER_API MatrixHandle matrix_get_diagonal(MatrixHandle a) {
    if (!a) return nullptr;
    try {
        auto A = as_matrix(a);
        if (!require_square(A)) return nullptr;

        VectorXd diag = A->data.diagonal();
        MatrixXd out_mat = diag;

        auto out = std::make_unique<MatrixHandleImpl>();
        out->data = std::move(out_mat);
        out->type = "diagonal";
        out->rows = A->rows;
        out->cols = 1;
        out->is_symmetric = false;
        out->trace_val = diag.sum();
        out->condition_number = -1.0;
        return reinterpret_cast<MatrixHandle>(out.release());
    } catch (...) {
        return nullptr;
    }
}

// ===== Блочные операции =====

/**
 * Материализация плотного блока по маскам строк и столбцов
 * Возвращает новую матрицу размером |rows| × |cols|
 */
QUANTUM_ANALYZER_API MatrixHandle matrix_extract_block(MatrixHandle m, GroupHandle rows, GroupHandle cols) {
    if (!m || !rows || !cols) return nullptr;
    try {
        auto M = as_matrix(m);
        auto R = as_group(rows);
        auto C = as_group(cols);
        if (!M || !R || !C) return nullptr;
        if (R->nbasis != M->rows) return nullptr;
        if (C->nbasis != M->cols) return nullptr;

        int nr = static_cast<int>(R->bits.count());
        int nc = static_cast<int>(C->bits.count());
        if (nr == 0 || nc == 0) return nullptr;

        MatrixXd block(nr, nc);
        int i_idx = 0;
        for_each_set_bit(R->bits, [&](int i) {
            int j_idx = 0;
            for_each_set_bit(C->bits, [&](int j) {
                block(i_idx, j_idx) = M->data(i, j);
                j_idx++;
            });
            i_idx++;
        });

        auto out = std::make_unique<MatrixHandleImpl>();
        out->data = std::move(block);
        out->type = "block";
        out->rows = nr;
        out->cols = nc;
        // Блок симметричен только если он квадратный, использует одинаковые маски строк/столбцов,
        // а исходная матрица симметрична.
        out->is_symmetric = (nr == nc && R->bits == C->bits && M->is_symmetric == 1) ? 1 : 0;
        out->trace_val = (nr == nc) ? out->data.trace() : std::numeric_limits<double>::quiet_NaN();
        out->condition_number = -1.0;
        return reinterpret_cast<MatrixHandle>(out.release());
    } catch (...) {
        return nullptr;
    }
}

/**
 * Сумма диагональных элементов по индексам из группы: sum_{i∈G} A_ii
 * Требования: квадратная матрица, размер битсета группы == размер матрицы
 */
QUANTUM_ANALYZER_API double matrix_block_trace(MatrixHandle m, GroupHandle g) {
    if (!m || !g) return std::numeric_limits<double>::quiet_NaN();
    try {
        auto M = as_matrix(m);
        auto G = as_group(g);
        if (!M || !G) return std::numeric_limits<double>::quiet_NaN();
        if (M->rows != M->cols) return std::numeric_limits<double>::quiet_NaN();
        if (G->nbasis != M->rows) return std::numeric_limits<double>::quiet_NaN();

        double sum = 0.0;
        for_each_set_bit(G->bits, [&](int i) {
            double v = M->data(i, i);
            if (!is_finite(v)) throw std::runtime_error("non-finite value in matrix_block_trace");
            sum += v;
        });
        return sum;
    } catch (...) {
        return std::numeric_limits<double>::quiet_NaN();
    }
}

/**
 * Сумма квадратов элементов по блоку: sum_{i∈rows, j∈cols} A_ij^2
 * Ключевой примитив для порядков связей Виберга
 */
QUANTUM_ANALYZER_API double matrix_block_sum_squares(MatrixHandle m, GroupHandle rows, GroupHandle cols) {
    if (!m || !rows || !cols) return std::numeric_limits<double>::quiet_NaN();
    try {
        auto M = as_matrix(m);
        auto R = as_group(rows);
        auto C = as_group(cols);
        if (!M || !R || !C) return std::numeric_limits<double>::quiet_NaN();
        if (R->nbasis != M->rows) return std::numeric_limits<double>::quiet_NaN();
        if (C->nbasis != M->cols) return std::numeric_limits<double>::quiet_NaN();

        double sum = 0.0;
        for_each_set_bit(R->bits, [&](int i) {
            for_each_set_bit(C->bits, [&](int j) {
                double v = M->data(i, j);
                if (!is_finite(v)) throw std::runtime_error("non-finite value in matrix_block_sum_squares");
                sum += v * v;
            });
        });
        return sum;
    } catch (...) {
        return std::numeric_limits<double>::quiet_NaN();
    }
}

/**
 * Билинейная форма Майера по блоку: sum_{i∈rows, j∈cols} P_ij * S_ji
 * Ключевой примитив для порядков связей Майера
 */
QUANTUM_ANALYZER_API double matrix_block_mayer_pair(MatrixHandle p, MatrixHandle s, GroupHandle rows, GroupHandle cols) {
    if (!p || !s || !rows || !cols) return std::numeric_limits<double>::quiet_NaN();
    try {
        auto P = as_matrix(p);
        auto S = as_matrix(s);
        auto R = as_group(rows);
        auto C = as_group(cols);
        if (!P || !S || !R || !C) return std::numeric_limits<double>::quiet_NaN();
        if (P->rows != S->rows || P->cols != S->cols) return std::numeric_limits<double>::quiet_NaN();
        if (R->nbasis != P->rows) return std::numeric_limits<double>::quiet_NaN();
        if (C->nbasis != P->cols) return std::numeric_limits<double>::quiet_NaN();

        double sum = 0.0;
        for_each_set_bit(R->bits, [&](int i) {
            for_each_set_bit(C->bits, [&](int j) {
                double v1 = P->data(i, j);
                double v2 = S->data(j, i);
                if (!is_finite(v1) || !is_finite(v2)) throw std::runtime_error("non-finite in mayer_pair");
                sum += v1 * v2;
            });
        });
        return sum;
    } catch (...) {
        return std::numeric_limits<double>::quiet_NaN();
    }
}

/**
 * Норма Фробениуса блока: ||A_block||_F
 */
QUANTUM_ANALYZER_API double matrix_block_norm_fro(MatrixHandle m, GroupHandle rows, GroupHandle cols) {
    double sum_sq = matrix_block_sum_squares(m, rows, cols);
    if (!is_finite(sum_sq) || sum_sq < 0.0) return std::numeric_limits<double>::quiet_NaN();
    return std::sqrt(sum_sq);
}

// ===== Групповые операции (битовые маски) =====

QUANTUM_ANALYZER_API GroupHandle group_create(int nbasis) {
    if (nbasis <= 0) return nullptr;
    try {
        auto g = std::make_unique<GroupHandleImpl>(nbasis);
        // все биты 0 по умолчанию
        return reinterpret_cast<GroupHandle>(g.release());
    } catch (...) {
        return nullptr;
    }
}

QUANTUM_ANALYZER_API GroupHandle group_create_full(int nbasis) {
    if (nbasis <= 0) return nullptr;
    try {
        auto g = std::make_unique<GroupHandleImpl>(nbasis);
        g->bits.set(); // все биты 1
        return reinterpret_cast<GroupHandle>(g.release());
    } catch (...) {
        return nullptr;
    }
}

QUANTUM_ANALYZER_API GroupHandle group_from_indices(const int* idx, int count, int nbasis) {
    if (!idx || count < 0 || nbasis <= 0) return nullptr;
    try {
        // Предварительная проверка индексов
        for (int k = 0; k < count; ++k) {
            int i = idx[k];
            if (i < 0 || i >= nbasis) return nullptr;
        }

        auto g = std::make_unique<GroupHandleImpl>(nbasis);
        for (int k = 0; k < count; ++k) {
            g->bits.set(static_cast<size_t>(idx[k]));
        }
        return reinterpret_cast<GroupHandle>(g.release());
    } catch (...) {
        return nullptr;
    }
}

QUANTUM_ANALYZER_API GroupHandle group_from_atom(GaussianFileHandle file, int atom_idx) {
    if (!file) return nullptr;
    try {
        extern int gaussian_get_basis_size(GaussianFileHandle);
        extern int gaussian_get_ao_atom_mapping(GaussianFileHandle, int*, int);
        extern int gaussian_get_num_atoms(GaussianFileHandle);

        int nbasis = gaussian_get_basis_size(file);
        int natoms = gaussian_get_num_atoms(file);
        if (nbasis <= 0 || natoms <= 0) return nullptr;
        if (atom_idx < 0 || atom_idx >= natoms) return nullptr;

        std::vector<int> ao2atom(nbasis);
        if (gaussian_get_ao_atom_mapping(file, ao2atom.data(), nbasis) != 0) return nullptr;

        auto g = std::make_unique<GroupHandleImpl>(nbasis);
        for (int mu = 0; mu < nbasis; ++mu) {
            int a = ao2atom[mu];
            if (a == atom_idx) g->bits.set(static_cast<size_t>(mu));
        }
        return reinterpret_cast<GroupHandle>(g.release());
    } catch (...) {
        return nullptr;
    }
}

QUANTUM_ANALYZER_API void group_set_bit(GroupHandle g, int index, int value) {
    auto gg = as_group(g);
    if (!gg) return;
    if (index < 0 || index >= gg->nbasis) return;
    if (value) gg->bits.set(static_cast<size_t>(index));
    else gg->bits.reset(static_cast<size_t>(index));
}

QUANTUM_ANALYZER_API int group_get_bit(GroupHandle g, int index) {
    auto gg = as_group(g);
    if (!gg) return -1;
    if (index < 0 || index >= gg->nbasis) return -1;
    return gg->bits.test(static_cast<size_t>(index)) ? 1 : 0;
}

QUANTUM_ANALYZER_API int group_count(GroupHandle g) {
    auto gg = as_group(g);
    if (!gg) return -1;
    return static_cast<int>(gg->bits.count());
}

QUANTUM_ANALYZER_API int group_nbasis(GroupHandle g) {
    auto gg = as_group(g);
    if (!gg) return -1;
    return gg->nbasis;
}

QUANTUM_ANALYZER_API GroupHandle group_or(GroupHandle a, GroupHandle b) {
    auto ga = as_group(a);
    auto gb = as_group(b);
    if (!ga || !gb) return nullptr;
    if (ga->nbasis != gb->nbasis) return nullptr;
    try {
        auto out = std::make_unique<GroupHandleImpl>(ga->nbasis);
        out->bits = ga->bits | gb->bits;
        return reinterpret_cast<GroupHandle>(out.release());
    } catch (...) {
        return nullptr;
    }
}

QUANTUM_ANALYZER_API GroupHandle group_and(GroupHandle a, GroupHandle b) {
    auto ga = as_group(a);
    auto gb = as_group(b);
    if (!ga || !gb) return nullptr;
    if (ga->nbasis != gb->nbasis) return nullptr;
    try {
        auto out = std::make_unique<GroupHandleImpl>(ga->nbasis);
        out->bits = ga->bits & gb->bits;
        return reinterpret_cast<GroupHandle>(out.release());
    } catch (...) {
        return nullptr;
    }
}

QUANTUM_ANALYZER_API void group_free(GroupHandle g) {
    if (!g) return;
    delete as_group(g);
}

// ===== Спектральные операции =====

/**
 * Симметричное возведение в степень через диагонализацию: A^p = U·Λ^p·U^T
 * Требования:
 *   - Матрица квадратная и симметричная (проверяется по флагу)
 *   - Для p < 0: строго положительно определённая (все λ > eps)
 *   - Для дробных p и отрицательных λ: возвращает nullptr
 * Параметр eps используется только как порог для "почти нулевых" собственных значений
 */
QUANTUM_ANALYZER_API MatrixHandle matrix_symm_pow(MatrixHandle a, double exponent, double eps) {
    if (!a) return nullptr;
    if (!is_finite(exponent) || !is_finite(eps) || eps < 0.0) return nullptr;
    try {
        auto A = as_matrix(a);
        if (!require_square(A) || !require_symmetric(A)) return nullptr;

        int n = A->rows;
        if (n <= 0) return nullptr;

        // Копируем в column-major формат для LAPACK.
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> Acm = A->data;
        std::vector<double> w(static_cast<size_t>(n));

        char jobz = 'V';
        char uplo = 'U';
        int lda = n;
        int info = 0;

        // Запрос рабочего пространства LAPACK.
        int lwork = -1;
        double wkopt = 0.0;
        dsyev_(&jobz, &uplo, &n, Acm.data(), &lda, w.data(), &wkopt, &lwork, &info);
        if (info != 0) return nullptr;

        lwork = static_cast<int>(wkopt);
        if (lwork < 1) lwork = 1;
        std::vector<double> work(static_cast<size_t>(lwork));

        // Основной вызов LAPACK.
        dsyev_(&jobz, &uplo, &n, Acm.data(), &lda, w.data(), work.data(), &lwork, &info);
        if (info != 0) return nullptr;

        // Проверка собственных значений для отрицательных степеней
        if (exponent < 0.0) {
            for (int i = 0; i < n; ++i) {
                if (!is_finite(w[static_cast<size_t>(i)]) || w[static_cast<size_t>(i)] <= eps) {
                    return nullptr;
                }
            }
        }

        // Дробная степень отрицательных собственных значений дала бы комплексный результат.
        double rintp = std::round(exponent);
        bool p_is_integer = std::fabs(exponent - rintp) < 1e-12;
        for (int i = 0; i < n; ++i) {
            double lam = w[static_cast<size_t>(i)];
            if (!is_finite(lam)) return nullptr;
            if (lam < 0.0 && !p_is_integer) return nullptr;
        }

        // Вычисление Λ^p
        std::vector<double> wpow(static_cast<size_t>(n));
        for (int i = 0; i < n; ++i) {
            wpow[static_cast<size_t>(i)] = std::pow(w[static_cast<size_t>(i)], exponent);
            if (!is_finite(wpow[static_cast<size_t>(i)])) return nullptr;
        }

        // Сборка результата: U * diag(wpow) * U^T.
        MatrixXd U = Acm;
        MatrixXd D = MatrixXd::Zero(n, n);
        for (int i = 0; i < n; ++i) D(i, i) = wpow[static_cast<size_t>(i)];
        MatrixXd R = U * D * U.transpose();

        // Симметризация устраняет численный шум после спектральной сборки.
        R = (R + R.transpose()) * 0.5;

        auto out = std::make_unique<MatrixHandleImpl>();
        out->data = std::move(R);
        out->type = "symm_pow";
        out->rows = n;
        out->cols = n;
        out->is_symmetric = 1;  // Результат симметричен по построению.
        out->trace_val = R.trace();
        out->condition_number = -1.0;
        return reinterpret_cast<MatrixHandle>(out.release());
    } catch (...) {
        return nullptr;
    }
}

}
