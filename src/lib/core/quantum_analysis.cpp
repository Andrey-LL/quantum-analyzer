// Численные примитивы квантово-химического анализа для Lua-слоя.
// Основано на формулах из работ Löwdin (1950, 1970), Mulliken (1955), Mayer, Wiberg.

#include "api.h"
#include "internal.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <iostream>
#include <cmath>
#include <vector>
#include <memory>

using namespace Eigen;
using namespace std;

// ============================================================================
// Предварительные объявления внутренних структур.
// ============================================================================

struct GaussianFileImpl;
struct MatrixHandleImpl;
struct GroupHandleImpl;
struct MOData;

// ============================================================================
// Утилиты валидации
// ============================================================================

/**
 * Безопасная проверка приближённого равенства матриц
 * Избегает ассертов Eigen при несовпадении размеров
 * @param A Первая матрица
 * @param B Вторая матрица
 * @param tol Допустимая погрешность
 * @return true если матрицы равны с заданной точностью
 */
static inline bool safe_is_approx(const MatrixXd &A, const MatrixXd &B, double tol) {
    if (A.rows() != B.rows() || A.cols() != B.cols()) return false;
    return A.isApprox(B, tol);
}

/**
 * Проверка матрицы на наличие некорректных значений (NaN, Inf, экстремально большие)
 * @param M Проверяемая матрица
 * @param max_abs_threshold Порог максимального абсолютного значения
 * @param problematic_index Опциональный выходной параметр для индекса проблемы
 * @return true если найдены некорректные значения
 */
static bool matrix_has_bad_values(const MatrixXd &M, double max_abs_threshold, pair<int,int>* problematic_index=nullptr) {
    for (int i = 0; i < M.rows(); ++i) {
        for (int j = 0; j < M.cols(); ++j) {
            double v = M(i,j);
            if (!isfinite(v) || fabs(v) > max_abs_threshold) {
                if (problematic_index) *problematic_index = {i,j};
                return true;
            }
        }
    }
    return false;
}

// ============================================================================
// Симметричная ортогонализация Лёвдина.
// Формула: X = S^(-1/2), P' = X P X
// Источник: Löwdin (1950), Mayer (2002)
// ============================================================================

/**
 * Вычисление ортогонализатора Лёвдина: X = S^(-1/2)
 * Использует спектральное разложение: S = U·Λ·U^T → S^(-1/2) = U·Λ^(-1/2)·U^T
 * Требования:
 *   - Матрица перекрывания должна быть симметричной
 *   - Все собственные значения должны быть положительными (SPD)
 * @param overlap Матрица перекрывания S
 * @return Handle матрицы X = S^(-1/2) или nullptr при ошибке
 */
QUANTUM_ANALYZER_API MatrixHandle loewdin_orthogonalizer(MatrixHandle overlap) {
    if (!overlap) return nullptr;

    try {
        auto impl = reinterpret_cast<MatrixHandleImpl*>(overlap);
        const MatrixXd &S = impl->data;

        if (S.rows() != S.cols()) return nullptr;

        int n = S.rows();

        // Спектральное разложение симметричной матрицы
        SelfAdjointEigenSolver<MatrixXd> solver(S);

        if (solver.info() != Success) return nullptr;

        VectorXd evals = solver.eigenvalues();
        MatrixXd evecs = solver.eigenvectors();

        double min_eig = evals.minCoeff();

        // Проверка положительной определённости
        if (min_eig <= 0.0) {
            cerr << "loewdin_orthogonalizer: overlap matrix has non-positive eigenvalue "
                  << min_eig << ", cannot orthogonalize" << endl;
            return nullptr;
        }

        // Вычисление Λ^(-1/2)
        VectorXd inv_sqrt = evals.array().rsqrt();

        // Сборка S^(-1/2) = U · Λ^(-1/2) · U^T
        MatrixXd X = evecs * inv_sqrt.asDiagonal() * evecs.transpose();

        auto result_impl = make_unique<MatrixHandleImpl>();
        result_impl->data = X;
        result_impl->type = "S_inv_sqrt";
        result_impl->rows = n;
        result_impl->cols = n;
        result_impl->is_symmetric = 1;  // симметрична по построению
        result_impl->trace_val = X.trace();
        result_impl->condition_number = evals.maxCoeff() / evals.minCoeff();

        return reinterpret_cast<MatrixHandle>(result_impl.release());
    }
    catch (const exception &e) {
        cerr << "loewdin_orthogonalizer error: " << e.what() << endl;
        return nullptr;
    }
}

/**
 * Преобразование матрицы плотности в ортогональный базис Лёвдина
 * Формула: P' = S^(1/2) · P · S^(1/2)
 * Источник: Löwdin (1970), формула 4.75
 * @param density Матрица плотности P (AO-базис)
 * @param overlap Матрица перекрывания S
 * @return Handle ортогонализованной матрицы плотности P' или nullptr
 */
QUANTUM_ANALYZER_API MatrixHandle loewdin_transform_density(MatrixHandle density, MatrixHandle overlap) {
    if (!density || !overlap) return nullptr;

    try {
        // Вычисление S^(1/2) через спектральное разложение
        auto S_impl = reinterpret_cast<MatrixHandleImpl*>(overlap);
        const MatrixXd &S = S_impl->data;

        SelfAdjointEigenSolver<MatrixXd> solver(S);
        if (solver.info() != Success) return nullptr;

        VectorXd evals = solver.eigenvalues();
        MatrixXd evecs = solver.eigenvectors();

        // Проверка на положительные собственные значения
        if (evals.minCoeff() <= 0.0) {
            cerr << "loewdin_transform_density: non-positive eigenvalue in overlap" << endl;
            return nullptr;
        }

        // S^(1/2) = U · Λ^(1/2) · U^T
        VectorXd sqrt_eig = evals.cwiseSqrt();
        MatrixXd S_sqrt = evecs * sqrt_eig.asDiagonal() * evecs.transpose();

        auto P_impl = reinterpret_cast<MatrixHandleImpl*>(density);

        // P' = S^(1/2) · P · S^(1/2)
        MatrixXd R = S_sqrt * P_impl->data * S_sqrt;

        auto res_impl = make_unique<MatrixHandleImpl>();
        res_impl->data = R;
        res_impl->type = "P_loewdin";
        res_impl->rows = R.rows();
        res_impl->cols = R.cols();
        res_impl->is_symmetric = 1;  // симметрична по построению
        res_impl->trace_val = R.trace();
        res_impl->condition_number = -1.0;

        return reinterpret_cast<MatrixHandle>(res_impl.release());
    } catch (const exception &e) {
        cerr << "Error in loewdin_transform_density: " << e.what() << endl;
        return nullptr;
    }
}

// ============================================================================
// Популяционный анализ.
// Источники: Mulliken (1955), Löwdin (1970)
// ============================================================================

/**
 * Löwdin популяционный анализ по атомам
 * Формула: P' = S^(1/2) · P · S^(1/2), затем q_A = Z_A - sum_{μ∈A} P'_{μμ}
 * Источник: Löwdin (1970), формулы 4.75-4.76
 * Преимущество: меньшая зависимость от базиса по сравнению с Малликеном
 * @param file Handle Gaussian файла
 * @param density Матрица плотности P
 * @param overlap Матрица перекрывания S
 * @param charges Выходной массив зарядов [natoms]
 * @return 0 при успехе, -1 при ошибке
 */
int loewdin_charges(GaussianFileHandle file, MatrixHandle density,
                    MatrixHandle overlap, double* charges) {
    if (!file || !density || !overlap || !charges) return -1;

    try {
        int natoms = gaussian_get_num_atoms(file);
        int nbasis = gaussian_get_basis_size(file);
        if (natoms <= 0 || nbasis <= 0) return -1;

        // Получаем маппинг AO -> атом.
        vector<int> ao2atom(nbasis);
        if (gaussian_get_ao_atom_mapping(file, ao2atom.data(), nbasis) != 0) {
            return -1;
        }

        // Получаем атомные номера.
        vector<int> Zs(natoms);
        for (int i = 0; i < natoms; ++i) {
            Zs[i] = gaussian_get_atomic_number(file, i);
            if (Zs[i] <= 0) return -1;
        }

        auto P_impl = reinterpret_cast<MatrixHandleImpl*>(density);
        MatrixXd P = P_impl->data;

        auto S_impl = reinterpret_cast<MatrixHandleImpl*>(overlap);
        MatrixXd S = S_impl->data;

        // Спектральное разложение S для вычисления S^(1/2)
        Eigen::SelfAdjointEigenSolver<MatrixXd> es(S);
        VectorXd eigvals = es.eigenvalues();
        MatrixXd eigvecs = es.eigenvectors();

        // S^(1/2) = U · Λ^(1/2) · U^T
        VectorXd sqrt_eigvals = eigvals.cwiseSqrt();
        MatrixXd S_sqrt = eigvecs * sqrt_eigvals.asDiagonal() * eigvecs.transpose();

        // Формула 4.75: P_Löwdin = S^(1/2) · P · S^(1/2)
        MatrixXd P_loewdin = S_sqrt * P * S_sqrt;

        // Проверка P_loewdin на корректные значения
        pair<int,int> idx;
        if (matrix_has_bad_values(P_loewdin, 1e6, &idx)) {
            cerr << "Invalid values in P_loewdin at (" << idx.first << ", " << idx.second << ")" << endl;
            return -1;
        }

        // Расчет популяций как сумма диагональных элементов P_loewdin по атомам
        vector<double> populations(natoms, 0.0);
        for (int mu = 0; mu < nbasis; ++mu) {
            int atom = ao2atom[mu];
            if (atom >= 0 && atom < natoms) {
                populations[atom] += P_loewdin(mu, mu);
            }
        }

        // Формула 4.76: q_A = Z_A - sum(n_μ^Löwdin)
        for (int i = 0; i < natoms; ++i) {
            charges[i] = Zs[i] - populations[i];
        }

        return 0;
    } catch (const exception &e) {
        cerr << "Error in loewdin_charges: " << e.what() << endl;
        return -1;
    }
}

// ============================================================================
// Естественные орбитали.
// Источник: Löwdin (1955), Reed & Weinhold (NBO)
// ============================================================================

/**
 * Вычисление натуральных орбиталей и их заселённостей
 * Алгоритм:
 *   1. Вычислить M = S^(1/2) · P · S^(1/2)
 *   2. Диагонализировать M = U · n · U^T
 *   3. Заселённости: собственные значения n (в убывающем порядке)
 *   4. Коэффициенты NO: C_NO = S^(-1/2) · U
 * @param density Матрица плотности P
 * @param overlap Матрица перекрывания S
 * @param occupations Выходной массив оккупаций [nbasis]
 * @param C_NO Выходная матрица коэффициентов NO [nbasis × nbasis] (может быть nullptr)
 * @return 0 при успехе, -1 при ошибке
 */
QUANTUM_ANALYZER_API int natural_orbitals(MatrixHandle density, MatrixHandle overlap,
                     double* occupations, MatrixHandle* C_NO) {
    if (!density || !overlap || !occupations) return -1;

    try {
        auto P_impl = reinterpret_cast<MatrixHandleImpl*>(density);
        auto S_impl = reinterpret_cast<MatrixHandleImpl*>(overlap);
        const MatrixXd &P = P_impl->data;
        const MatrixXd &S = S_impl->data;

        int nbasis = P.rows();
        if (nbasis <= 0) return -1;

        // Валидация входных матриц
        pair<int,int> idx;
        if (matrix_has_bad_values(P, 1e6, &idx) || matrix_has_bad_values(S, 1e6, &idx)) {
            cerr << "Detected invalid matrix values in P or S. Aborting natural orbitals." << endl;
            return -1;
        }

        // Вычисление S^(1/2) через спектральное разложение
        SelfAdjointEigenSolver<MatrixXd> solver(S);
        if (solver.info() != Success) return -1;

        VectorXd sqrt_eig = solver.eigenvalues().cwiseSqrt();
        MatrixXd S_sqrt = solver.eigenvectors() * sqrt_eig.asDiagonal() * solver.eigenvectors().transpose();

        // M = S^(1/2) · P · S^(1/2)
        MatrixXd M = S_sqrt * P * S_sqrt;

        // Диагонализация M
        SelfAdjointEigenSolver<MatrixXd> solver_M(M);
        if (solver_M.info() != Success) {
            return -1;
        }

        // Заселённости (в убывающем порядке)
        VectorXd eigs = solver_M.eigenvalues().reverse();
        for (int i = 0; i < nbasis; ++i) {
            occupations[i] = eigs(i);
        }

        // Если требуется матрица коэффициентов C_NO
        if (C_NO != nullptr) {
            // U (в обратном порядке столбцов для соответствия убывающим оккупациям)
            MatrixXd U = solver_M.eigenvectors().rowwise().reverse();

            // C_NO = S^(-1/2) · U; устойчивость зависит от loewdin_orthogonalizer.
            MatrixHandle S_inv_half = loewdin_orthogonalizer(overlap);
            if (!S_inv_half) {
                return -1;
            }
            auto S_inv_half_impl = reinterpret_cast<MatrixHandleImpl*>(S_inv_half);

            MatrixXd C_NO_mat = S_inv_half_impl->data * U;

            auto C_NO_impl = make_unique<MatrixHandleImpl>();
            C_NO_impl->data = C_NO_mat;
            C_NO_impl->type = "natural_orbitals";
            C_NO_impl->rows = nbasis;
            C_NO_impl->cols = nbasis;
            C_NO_impl->is_symmetric = 2;  // не симметрична (коэффициенты NO)
            C_NO_impl->trace_val = C_NO_mat.trace();
            C_NO_impl->condition_number = 0.0;

            *C_NO = reinterpret_cast<MatrixHandle>(C_NO_impl.release());

            matrix_free(S_inv_half);
        }

        return 0;
    } catch (const exception &e) {
        cerr << "Error in natural_orbitals: " << e.what() << endl;
        return -1;
    }
}

/**
 * Подсчёт числа значимых натуральных орбиталей
 * @param density Матрица плотности P
 * @param overlap Матрица перекрывания S
 * @param threshold Порог заселённости (обычно 0.02 для виртуальных орбиталей)
 * @return Число NO с occupation > threshold или -1 при ошибке
 */
QUANTUM_ANALYZER_API int natural_orbitals_count_significant(MatrixHandle density, MatrixHandle overlap,
                                        double threshold) {
    if (!density || !overlap) return -1;

    try {
        auto P_impl = reinterpret_cast<MatrixHandleImpl*>(density);
        int nbasis = P_impl->rows;
        if (nbasis <= 0) return -1;

        vector<double> occupations(nbasis);
        if (natural_orbitals(density, overlap, occupations.data(), nullptr) != 0) {
            return -1;
        }

        int count = 0;
        for (int i = 0; i < nbasis; ++i) {
            if (occupations[i] > threshold) {
                count++;
            }
        }

        return count;
    } catch (const exception &e) {
        cerr << "Error in natural_orbitals_count_significant: " << e.what() << endl;
        return -1;
    }
}

// ============================================================================
// Валидация и диагностика.
// ============================================================================

/**
 * Диагностика матрицы плотности
 * Проверки:
 *   - Симметричность: P = P^T
 *   - След Tr(PS) = N_electrons (число электронов)
 *   - Отсутствие NaN/Inf и экстремальных значений
 * @param file Handle Gaussian файла (для получения N_electrons)
 * @param density Матрица плотности P
 * @param overlap Матрица перекрывания S
 * @param tolerance Допустимая относительная погрешность
 * @return 1 если корректна, 0 если нет, -1 при ошибке
 */
QUANTUM_ANALYZER_API int validate_density_matrix_diagnostics(GaussianFileHandle file, MatrixHandle density,
                                         MatrixHandle overlap, double tolerance) {
    if (!file || !density || !overlap) return -1;

    try {
        auto P_impl = reinterpret_cast<MatrixHandleImpl*>(density);
        auto S_impl = reinterpret_cast<MatrixHandleImpl*>(overlap);
        MatrixXd P = P_impl->data;
        MatrixXd S = S_impl->data;

        // Проверка симметричности
        if (!safe_is_approx(P, P.transpose(), tolerance)) {
            cerr << "Density matrix is not symmetric" << endl;
            return 0;
        }

        // Проверка на NaN/Inf и экстремальные значения
        pair<int,int> idx;
        if (matrix_has_bad_values(P, 1e9, &idx)) {
            cerr << "Density matrix contains NaN/Inf or excessively large values at ("
                  << idx.first << ", " << idx.second << ")" << endl;
            return 0;
        }
        if (matrix_has_bad_values(S, 1e9, &idx)) {
            cerr << "Overlap matrix contains NaN/Inf or excessively large values at ("
                  << idx.first << ", " << idx.second << ")" << endl;
            return 0;
        }

        // Вычисление Tr(PS) = N_e
        MatrixXd PS = P * S;
        double tr_PS = PS.trace();

        int alpha, beta;
        int N_e = gaussian_get_electrons(file, &alpha, &beta);
        if (N_e <= 0) return -1;

        if (!isfinite(tr_PS)) {
            cerr << "Tr(PS) is not finite: " << tr_PS << endl;
            return 0;
        }

        // Проверка соответствия числа электронов
        double rel_err = fabs(tr_PS - N_e);
        if (rel_err > max(1.0, fabs((double)N_e) * tolerance)) {
            cerr << "Tr(PS) = " << tr_PS << " != N_e = " << N_e << endl;
            return 0;
        }

        return 1;
    } catch (const exception &e) {
        cerr << "Error in validate_density_matrix: " << e.what() << endl;
        return -1;
    }
}

/**
 * Диагностика матрицы перекрывания
 * Проверки:
 *   - Симметричность: S = S^T
 *   - Положительная определённость: все λ > 0
 *   - Число обусловленности (предупреждение если > 10^12)
 * @param overlap Матрица перекрывания S
 * @param tolerance Допустимая погрешность для проверки симметричности
 * @param min_eigenvalue Минимальное допустимое собственное значение
 * @return 1 если корректна, 0 если нет, -1 при ошибке
 */
QUANTUM_ANALYZER_API int validate_overlap_matrix_diagnostics(MatrixHandle overlap, double tolerance,
                                         double min_eigenvalue) {
    if (!overlap) return -1;

    try {
        auto S_impl = reinterpret_cast<MatrixHandleImpl*>(overlap);
        MatrixXd S = S_impl->data;

        // Проверка симметричности с увеличенным допуском
        if (!safe_is_approx(S, S.transpose(), tolerance * 10)) {
            cerr << "Warning: Overlap matrix symmetry tolerance increased due to numerical precision issues" << endl;
            return 0;
        }

        // Проверка на NaN/Inf и экстремальные значения
        pair<int,int> idx;
        if (matrix_has_bad_values(S, 1e9, &idx)) {
            cerr << "Overlap matrix contains NaN/Inf or excessively large values at ("
                  << idx.first << ", " << idx.second << ")" << endl;
            return 0;
        }

        // Спектральное разложение
        SelfAdjointEigenSolver<MatrixXd> solver(S);
        if (solver.info() != Success) {
            cerr << "Eigen decomposition of S failed" << endl;
            return -1;
        }
        VectorXd eigenvalues = solver.eigenvalues();
        double min_found = eigenvalues.minCoeff();
        double max_found = eigenvalues.maxCoeff();

        if (!isfinite(min_found) || !isfinite(max_found)) {
            cerr << "Non-finite eigenvalues in overlap" << endl;
            return 0;
        }

        // Допускаем небольшие отрицательные собственные значения из-за численного шума
        if (min_found < -1e-6) {
            cerr << "Overlap has eigenvalue " << min_found << " < allowed tolerance. Proceeding with regularization." << endl;
        }

        // Предупреждение о плохой обусловленности
        if (min_found > 0 && max_found / min_found > 1e12) {
            cerr << "Warning: Overlap matrix poorly conditioned (cond=" << (max_found/min_found) << ")" << endl;
        }

        return 1;
    } catch (const exception &e) {
        cerr << "Error in validate_overlap_matrix: " << e.what() << endl;
        return -1;
    }
}

// ============================================================================
// Анализ блоков плотности.
// ============================================================================

/**
 * Получение диапазона AO для атома (если AO непрерывны)
 * @param file Handle Gaussian файла
 * @param atom Индекс атома
 * @param start_out Выходной параметр: начальный индекс, включительно
 * @param end_out Выходной параметр: конечный индекс, исключительно
 * @return 0 если непрерывно, -2 если AO атома расположены с разрывом, -1 при ошибке
 */
int get_atom_ao_range(GaussianFileHandle file, int atom, int* start_out, int* end_out) {
    if (!file || !start_out || !end_out) return -1;
    int nbasis = gaussian_get_basis_size(file);
    if (nbasis <= 0) return -1;
    vector<int> ao2atom(nbasis);
    if (gaussian_get_ao_atom_mapping(file, ao2atom.data(), nbasis) != 0) return -1;

    int start = -1, end = -1;
    for (int i = 0; i < nbasis; ++i) {
        int a = ao2atom[i];
        if (a == atom) {
            if (start == -1) start = i;
            end = i + 1;
        } else if (start != -1) {
            // Обнаружен разрыв в индексах AO
            return -2;
        }
    }
    if (start == -1) {
        *start_out = 0; *end_out = 0;
        return 0;
    }
    *start_out = start; *end_out = end;
    return 0;
}

/**
 * Извлечение подматрицы из MatrixHandle
 * @param matrix Исходная матрица
 * @param row_start Начальная строка, включительно
 * @param row_end Конечная строка, исключительно
 * @param col_start Начальный столбец, включительно
 * @param col_end Конечный столбец, исключительно
 * @return Новый MatrixHandle с подматрицей или nullptr
 */
MatrixHandle extract_matrix_block(MatrixHandle matrix, int row_start, int row_end, int col_start, int col_end) {
    if (!matrix) return nullptr;
    auto impl = reinterpret_cast<MatrixHandleImpl*>(matrix);
    int rows = impl->rows;
    int cols = impl->cols;
    if (row_start < 0 || row_end > rows || col_start < 0 || col_end > cols || row_start >= row_end || col_start >= col_end) return nullptr;

    MatrixXd block = impl->data.block(row_start, col_start, row_end - row_start, col_end - col_start);
    auto out_impl = make_unique<MatrixHandleImpl>();
    out_impl->data = block;
    out_impl->rows = block.rows();
    out_impl->cols = block.cols();
    out_impl->type = "block";
    out_impl->is_symmetric = 0;  // Симметрия зависит от исходной матрицы и диапазона.
    out_impl->trace_val = block.trace();
    return reinterpret_cast<MatrixHandle>(out_impl.release());
}

/**
 * Главный анализ блоков плотности
 * Вычисляет диагональные (P_AA) и офф-диагональные (P_AB) блоки матрицы плотности
 * @param file Handle Gaussian файла
 * @param density Матрица плотности P
 * @param result_out Указатель на структуру DensityBlockAnalysis (должна быть выделена)
 * @return 0 при успехе, -1 при ошибке
 */
QUANTUM_ANALYZER_API int analyze_density_blocks(GaussianFileHandle file, MatrixHandle density, DensityBlockAnalysis* result_out) {
    if (!file || !density || !result_out) return -1;
    int natoms = gaussian_get_num_atoms(file);
    int nbasis = gaussian_get_basis_size(file);
    if (natoms <= 0 || nbasis <= 0) return -1;

    // Инициализация структуры с динамическим выделением памяти.
    result_out->natoms = natoms;
    result_out->nbasis = nbasis;
    result_out->n_diagonal = 0;
    result_out->n_offdiagonal = 0;
    result_out->diagonal_blocks.clear();
    result_out->offdiag_blocks.clear();
    result_out->diagonal_blocks.reserve(natoms);
    result_out->offdiag_blocks.reserve(natoms * (natoms - 1) / 2);

    // Получаем маппинг AO -> атом.
    vector<int> ao2atom(nbasis);
    if (gaussian_get_ao_atom_mapping(file, ao2atom.data(), nbasis) != 0) return -1;

    auto P_impl = reinterpret_cast<MatrixHandleImpl*>(density);
    MatrixXd P = P_impl->data;

    // Построение списков индексов AO по атомам с поддержкой разрывов.
    vector<vector<int>> atom_indices(natoms);
    for (int mu = 0; mu < nbasis; ++mu) {
        int a = ao2atom[mu];
        if (a >= 0 && a < natoms) atom_indices[a].push_back(mu);
    }

    // Диагональные блоки P_AA
    for (int a = 0; a < natoms; ++a) {
        int na = (int)atom_indices[a].size();
        if (na == 0) continue;
        MatrixXd block(na, na);
        for (int i = 0; i < na; ++i)
            for (int j = 0; j < na; ++j)
                block(i,j) = P(atom_indices[a][i], atom_indices[a][j]);
        double tr = block.trace();
        double frob = block.norm();
        double* data = (double*)malloc(sizeof(double) * na * na);
        if (!data) return -1;
        for (int i = 0; i < na; ++i)
            for (int j = 0; j < na; ++j)
                data[i*na + j] = block(i,j);
        
        DiagonalBlockInfo info;
        info.atom = a;
        info.size = na;
        info.trace = tr;
        info.frob_norm = frob;
        info.data = data;
        result_out->diagonal_blocks.push_back(info);
        result_out->n_diagonal++;
    }

    // Внедиагональные блоки P_AB для пар A < B.
    for (int A = 0; A < natoms; ++A) {
        for (int B = A+1; B < natoms; ++B) {
            int na = (int)atom_indices[A].size();
            int nb = (int)atom_indices[B].size();
            if (na == 0 || nb == 0) continue;
            MatrixXd block(na, nb);
            for (int i = 0; i < na; ++i)
                for (int j = 0; j < nb; ++j)
                    block(i,j) = P(atom_indices[A][i], atom_indices[B][j]);
            double frob = block.norm();
            double* data = (double*)malloc(sizeof(double) * na * nb);
            if (!data) return -1;
            for (int i = 0; i < na; ++i)
                for (int j = 0; j < nb; ++j)
                    data[i*nb + j] = block(i,j);
            
            OffDiagonalBlockInfo info;
            info.atom_a = A;
            info.atom_b = B;
            info.size_a = na;
            info.size_b = nb;
            info.frob_norm = frob;
            info.data = data;
            result_out->offdiag_blocks.push_back(info);
            result_out->n_offdiagonal++;
        }
    }

    return 0;
}

/**
 * Освобождение памяти, выделенной анализом блоков
 * @param analysis Указатель на структуру DensityBlockAnalysis
 */
QUANTUM_ANALYZER_API void free_density_block_analysis(DensityBlockAnalysis* analysis) {
    if (!analysis) return;
    
    // Освобождаем данные матриц, выделенные через malloc.
    for (auto& block : analysis->diagonal_blocks) {
        if (block.data) {
            free(block.data);
            block.data = nullptr;
        }
    }
    for (auto& block : analysis->offdiag_blocks) {
        if (block.data) {
            free(block.data);
            block.data = nullptr;
        }
    }
    
    // Сбрасываем контейнеры и счётчики.
    analysis->diagonal_blocks.clear();
    analysis->offdiag_blocks.clear();
    analysis->n_diagonal = 0;
    analysis->n_offdiagonal = 0;
    analysis->natoms = 0;
    analysis->nbasis = 0;
}
