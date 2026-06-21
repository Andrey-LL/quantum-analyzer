// C API для интерфейса с LuaJIT FFI.

#ifndef QUANTUM_ANALYZER_API_H
#define QUANTUM_ANALYZER_API_H

#ifdef __cplusplus
extern "C" {
#endif

// Макрос экспорта использует настройку QUANTUM_ANALYZER_API_EXPORT из Xmake.
#ifdef _WIN32
    // Windows: __declspec(dllexport/dllimport).
    #ifdef QUANTUM_ANALYZER_API_EXPORT
        #define QUANTUM_ANALYZER_API __declspec(dllexport)
    #else
        #define QUANTUM_ANALYZER_API __declspec(dllimport)
    #endif
#else
    // Linux/macOS: visibility attribute управляет видимостью символов.
    #ifdef __GNUC__
        #define QUANTUM_ANALYZER_API __attribute__((visibility("default")))
    #else
        #define QUANTUM_ANALYZER_API
    #endif
#endif

typedef void* GaussianFileHandle;
typedef void* MatrixHandle;
typedef struct DensityBlockAnalysis DensityBlockAnalysis;

// Операции с файлом Gaussian.
QUANTUM_ANALYZER_API GaussianFileHandle gaussian_open(const char* filename);
QUANTUM_ANALYZER_API void gaussian_close(GaussianFileHandle file);
QUANTUM_ANALYZER_API int gaussian_get_basis_size(GaussianFileHandle file);
QUANTUM_ANALYZER_API const char* gaussian_get_method(GaussianFileHandle file);
QUANTUM_ANALYZER_API int gaussian_get_num_atoms(GaussianFileHandle file);
QUANTUM_ANALYZER_API int gaussian_get_electrons(GaussianFileHandle file, int* alpha, int* beta);

// Доступ к матрицам.
QUANTUM_ANALYZER_API MatrixHandle matrix_create(int rows, int cols, const double* data);
QUANTUM_ANALYZER_API MatrixHandle gaussian_get_matrix(GaussianFileHandle file, const char* matrix_type);
QUANTUM_ANALYZER_API int matrix_get_size(MatrixHandle matrix, int* rows, int* cols);
QUANTUM_ANALYZER_API double matrix_get_element(MatrixHandle matrix, int row, int col);
QUANTUM_ANALYZER_API double matrix_trace(MatrixHandle matrix);
QUANTUM_ANALYZER_API int matrix_is_symmetric(MatrixHandle matrix);
QUANTUM_ANALYZER_API double matrix_condition_number(MatrixHandle matrix);
QUANTUM_ANALYZER_API void matrix_free(MatrixHandle matrix);

// Линейная алгебра на Eigen.
QUANTUM_ANALYZER_API int matrix_eigenvalues(MatrixHandle matrix, double* eigenvalues);
QUANTUM_ANALYZER_API int matrix_eigensystem(MatrixHandle matrix, double* eigenvalues, MatrixHandle* eigenvectors);
QUANTUM_ANALYZER_API int matrix_svd(MatrixHandle matrix, MatrixHandle* U, MatrixHandle* S, MatrixHandle* Vt);

QUANTUM_ANALYZER_API MatrixHandle matrix_power(MatrixHandle matrix, double exponent);
QUANTUM_ANALYZER_API MatrixHandle matrix_multiply(MatrixHandle a, MatrixHandle b);

// API для молекулярных орбиталей.
// Возвращает матрицу коэффициентов МО размером nbasis x nmo или NULL.
QUANTUM_ANALYZER_API MatrixHandle gaussian_get_mo_coefficients(GaussianFileHandle file);
// Заполняет буфер энергий размером size; возвращает полное число доступных энергий или -1 при ошибке.
QUANTUM_ANALYZER_API int gaussian_get_orbital_energies(GaussianFileHandle file, double* energies, int size);
// Заполняет ao2atom маппингом AO -> атом длиной nbasis. Возвращает 0 при успехе.
QUANTUM_ANALYZER_API int gaussian_get_ao_atom_mapping(GaussianFileHandle file, int* ao2atom, int nbasis);
// Возвращает выделенный массив C-строк; вызывающий код освобождает его через gaussian_free_ao_labels.
QUANTUM_ANALYZER_API const char** gaussian_get_ao_labels(GaussianFileHandle file);
QUANTUM_ANALYZER_API void gaussian_free_ao_labels(const char** labels, int count);

// Доступ к атомным номерам и координатам геометрии.
// Возвращает атомный номер Z для указанного атома или -1 при ошибке.
QUANTUM_ANALYZER_API int gaussian_get_atomic_number(GaussianFileHandle file, int atom_idx);
// Заполняет coords длиной 3*max_atoms и atomic_numbers длиной max_atoms. Возвращает natoms или -1 при ошибке.
QUANTUM_ANALYZER_API int gaussian_get_geometry_coordinates(GaussianFileHandle file, double* coords, int* atomic_numbers, int max_atoms);

// Преобразование AO -> MO: C^T * M * C.
QUANTUM_ANALYZER_API MatrixHandle matrix_ao_to_mo(MatrixHandle matrix_ao, MatrixHandle mo_coeff);

// Валидация матриц.
QUANTUM_ANALYZER_API int validate_overlap_matrix(MatrixHandle overlap);
QUANTUM_ANALYZER_API int validate_density_matrix(MatrixHandle density, MatrixHandle overlap, int total_electrons);
// Анализ Mulliken вынесен из парсера в quantum_analysis API.
QUANTUM_ANALYZER_API const char* gaussian_get_basis_name(GaussianFileHandle file);
QUANTUM_ANALYZER_API double gaussian_get_nuclear_repulsion(GaussianFileHandle file);

typedef void* GroupHandle;

// Операции с группами AO на битовых масках.
QUANTUM_ANALYZER_API GroupHandle group_create(int nbasis);                       // Создаёт пустую группу.
QUANTUM_ANALYZER_API GroupHandle group_create_full(int nbasis);                  // Создаёт группу со всеми установленными битами.
QUANTUM_ANALYZER_API GroupHandle group_from_indices(const int* idx, int count, int nbasis); // Создаёт группу из индексов AO.
QUANTUM_ANALYZER_API GroupHandle group_from_atom(GaussianFileHandle file, int atom_idx);    // Создаёт группу для одного атома.
QUANTUM_ANALYZER_API void group_set_bit(GroupHandle g, int index, int value);    // Устанавливает или сбрасывает бит; value: 0 или 1.
QUANTUM_ANALYZER_API int  group_get_bit(GroupHandle g, int index);               // Возвращает значение бита или -1 при ошибке.
QUANTUM_ANALYZER_API int  group_count(GroupHandle g);                            // Возвращает число установленных битов.
QUANTUM_ANALYZER_API int  group_nbasis(GroupHandle g);                           // Возвращает длину битовой маски.
QUANTUM_ANALYZER_API GroupHandle group_or(GroupHandle a, GroupHandle b);         // Объединение групп.
QUANTUM_ANALYZER_API GroupHandle group_and(GroupHandle a, GroupHandle b);        // Пересечение групп.
QUANTUM_ANALYZER_API void group_free(GroupHandle g);                             // Освобождает ресурсы группы.

// Базовая линейная алгебра.
QUANTUM_ANALYZER_API MatrixHandle matrix_add(MatrixHandle a, double alpha, MatrixHandle b, double beta); // C = α·A + β·B.
QUANTUM_ANALYZER_API MatrixHandle matrix_scale(MatrixHandle a, double alpha);    // B = α·A.
QUANTUM_ANALYZER_API MatrixHandle matrix_transpose(MatrixHandle a);              // B = A^T.
QUANTUM_ANALYZER_API MatrixHandle matrix_triple_product_symm(MatrixHandle a, MatrixHandle b); // C = A·B·A^T, результат симметричен при симметричной B.

// Поэлементные операции.
QUANTUM_ANALYZER_API MatrixHandle matrix_hadamard(MatrixHandle a, MatrixHandle b); // C_ij = A_ij * B_ij.
QUANTUM_ANALYZER_API MatrixHandle matrix_cwise_pow(MatrixHandle a, double exponent); // C_ij = A_ij^p.
QUANTUM_ANALYZER_API MatrixHandle matrix_threshold(MatrixHandle a, double theta); // C_ij = (|A_ij| >= θ) ? A_ij : 0.
QUANTUM_ANALYZER_API MatrixHandle matrix_clamp(MatrixHandle a, double lo, double hi); // C_ij = clamp(A_ij, lo, hi).

// Редукции матриц.
QUANTUM_ANALYZER_API double matrix_norm_fro(MatrixHandle a);                     // Норма Фробениуса: ||A||_F.
QUANTUM_ANALYZER_API double matrix_max_abs(MatrixHandle a);                      // Максимальный модуль элемента.
QUANTUM_ANALYZER_API double matrix_min_abs_nonzero(MatrixHandle a);              // Минимальный ненулевой модуль элемента.
QUANTUM_ANALYZER_API MatrixHandle matrix_get_diagonal(MatrixHandle a);           // Возвращает диагональ как матрицу N×1.

// Блочные операции через группы без копирования данных для редукций.
QUANTUM_ANALYZER_API MatrixHandle matrix_extract_block(MatrixHandle m, GroupHandle rows, GroupHandle cols); // Материализует плотный подблок.
QUANTUM_ANALYZER_API double matrix_block_trace(MatrixHandle m, GroupHandle g);   // Сумма диагональных элементов по индексам группы.
QUANTUM_ANALYZER_API double matrix_block_sum_squares(MatrixHandle m, GroupHandle rows, GroupHandle cols); // Σ_{i∈rows,j∈cols} m_ij².
QUANTUM_ANALYZER_API double matrix_block_mayer_pair(MatrixHandle p, MatrixHandle s, GroupHandle rows, GroupHandle cols); // Σ P_ij·S_ji.
QUANTUM_ANALYZER_API double matrix_block_norm_fro(MatrixHandle m, GroupHandle rows, GroupHandle cols); // Норма Фробениуса для блока.

// Степень симметричной матрицы через LAPACK dsyev.
QUANTUM_ANALYZER_API MatrixHandle matrix_symm_pow(MatrixHandle a, double exponent, double eps);

// Ортогонализация Лёвдина.
QUANTUM_ANALYZER_API MatrixHandle loewdin_orthogonalizer(MatrixHandle overlap);
QUANTUM_ANALYZER_API MatrixHandle loewdin_transform_density(MatrixHandle density, MatrixHandle overlap);

// Естественные орбитали.
QUANTUM_ANALYZER_API int natural_orbitals(MatrixHandle density, MatrixHandle overlap,
                     double* occupations, MatrixHandle* C_NO);
QUANTUM_ANALYZER_API int natural_orbitals_count_significant(MatrixHandle density, MatrixHandle overlap,
                                        double threshold);

// Валидация и диагностика.
QUANTUM_ANALYZER_API int validate_density_matrix_diagnostics(GaussianFileHandle file, MatrixHandle density,
                             MatrixHandle overlap, double tolerance);
QUANTUM_ANALYZER_API int validate_overlap_matrix_diagnostics(MatrixHandle overlap, double tolerance,
                             double min_eigenvalue);

// Анализ блоков матрицы плотности.
QUANTUM_ANALYZER_API int analyze_density_blocks(GaussianFileHandle file, MatrixHandle density, DensityBlockAnalysis* result_out);
QUANTUM_ANALYZER_API void free_density_block_analysis(DensityBlockAnalysis* analysis);

#ifdef __cplusplus
}
#endif

#endif
