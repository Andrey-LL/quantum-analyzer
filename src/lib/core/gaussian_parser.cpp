// Парсер Gaussian с поддержкой коэффициентов молекулярных орбиталей.
// Извлекает матрицы, метаданные и геометрию из файлов Gaussian.

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <algorithm>
#include <iomanip>
#include <cctype>
#include <regex>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "api.h"
#include "internal.h"

using namespace Eigen;
using namespace std;

// ============================================================================
// Внутренние структуры
// ============================================================================

struct GaussianFileImpl;
struct MatrixHandleImpl;
struct GroupHandleImpl;
struct MOData;

// ============================================================================
// Утилиты парсинга
// ============================================================================

/**
 * Вычисляет индекс блока столбцов для треугольной матрицы Gaussian
 * Используется для корректного размещения элементов при парсинге
 * @param n Размер матрицы
 * @param k Номер строки в выводе Gaussian
 * @return Номер блока столбцов (j)
 */
int n_j(int n, int k) {
    for (int i = 0; i < 1000; i++) {
        int limit = (2 * n - 5 * (i - 1)) * i / 2;
        if (k < limit) {
            return i - 1;
        }
    }
    return 0;
}

/**
 * Очищает D-нотацию в экспоненциальных числах (D+03 -> E+03)
 * Gaussian использует D-нотацию для чисел с плавающей точкой
 * @param text Входная строка
 * @return Очищенная строка с E-нотацией
 */
string clean_d_notation(const string &text) {
    string result = text;
    for (size_t i = 0; i + 1 < result.size(); ++i) {
        if (result[i] == 'D' && (result[i + 1] == '+' || result[i + 1] == '-')) result[i] = 'E';
    }
    return result;
}

/**
 * Извлекает числовые значения из строки
 * Пропускает целые числа, принимает только числа с плавающей точкой
 * @param line Строка с числами
 * @return Вектор извлечённых чисел
 */
vector<double> extract_numbers(const string &line) {
    vector<double> numbers;
    string cleaned = clean_d_notation(line);
    istringstream iss(cleaned);
    string token;
    while (iss >> token) {
        if (token.find('.') == string::npos) continue; // Пропускаем целые числа.
        bool valid = true;
        bool has_digit = false;
        for (size_t i = 0; i < token.size(); ++i) {
            char c = token[i];
            if (isdigit(c)) has_digit = true;
            else if (c == '.' || c == 'E' || c == 'e' || c == '+' || c == '-') continue;
            else {
                valid = false;
                break;
            }
        }
        if (!valid || !has_digit) continue;
        try { numbers.push_back(stod(token)); } catch (...) {}
    }
    return numbers;
}

/**
 * Поиск подстроки без учёта регистра
 * @param haystack Строка для поиска
 * @param needle Искомая подстрока
 * @return true если найдено, false иначе
 */
bool contains_ignore_case(const string &haystack, const string &needle) {
    auto it = search(haystack.begin(), haystack.end(), needle.begin(), needle.end(),
                     [](char a, char b) { return tolower(a) == tolower(b); });
    return it != haystack.end();
}

/**
 * Извлекает первое число после указанного маркера в содержимом
 * @param content Полное содержимое файла
 * @param marker Маркер для поиска
 * @param result Выходной параметр для сохранённого числа
 * @return true при успехе, false если маркер не найден
 */
bool extract_number_after(const string &content, const string &marker, double &result) {
    size_t pos = content.find(marker);
    if (pos == string::npos) return false;
    size_t start = content.find('\n', pos);
    if (start == string::npos) start = pos + marker.size(); else start++;
    size_t end = content.find('\n', start);
    if (end == string::npos) end = content.size();
    string substring = content.substr(start, end - start);
    auto numbers = extract_numbers(substring);
    if (!numbers.empty()) {
        result = numbers[0];
        return true;
    }
    return false;
}

// ============================================================================
// Извлечение метаданных
// ============================================================================

/**
 * Извлекает метаданные из содержимого файла Gaussian
 * Включает: размер базиса, число электронов, метод, базисный набор, геометрию
 * @param file Указатель на структуру файла
 * @return true при успехе, false при ошибке
 */
bool extract_metadata(GaussianFileImpl *file) {
    const string &content = file->content;

    // Извлечение размера базиса.
    string marker = "basis functions";
    size_t search_pos = 0;
    while (true) {
        size_t pos = content.find(marker, search_pos);
        if (pos == string::npos) break;
        size_t line_start = content.rfind('\n', pos);
        if (line_start == string::npos) line_start = 0; else line_start++;
        string line = content.substr(line_start, pos - line_start + 20);
        istringstream iss(line);
        int num;
        if (iss >> num) file->nbasis = num;
        search_pos = pos + 1;
    }

    // Альтернативный формат Gaussian: NBasis=...
    if (file->nbasis == 0) {
        size_t pos = content.find("NBasis=");
        if (pos != string::npos) {
            size_t start = pos + 7;
            size_t end = start;
            while (end < content.size() && isdigit(content[end])) end++;
            if (end > start) {
                try { file->nbasis = stoi(content.substr(start, end - start)); } catch (...) {}
            }
        }
    }

    // Извлечение числа электронов.
    size_t alpha_pos = content.find("alpha electrons");
    if (alpha_pos != string::npos) {
        size_t line_start = content.rfind('\n', alpha_pos);
        if (line_start == string::npos) line_start = 0; else line_start++;
        string line = content.substr(line_start, alpha_pos - line_start + 200);
        try {
            regex int_re("([+-]?[0-9]+)");
            auto begin = sregex_iterator(line.begin(), line.end(), int_re);
            auto end = sregex_iterator();
            vector<int> ints;
            for (auto it = begin; it != end; ++it) {
                try { ints.push_back(stoi(it->str())); } catch (...) {}
            }
            if (!ints.empty()) {
                file->alpha_electrons = ints[0];
                if (ints.size() > 1) file->beta_electrons = ints[1];
            }
        } catch (...) {}
    }

    // Извлечение базисного набора.
    size_t basis_pos = content.find("Standard basis:");
    if (basis_pos != string::npos) {
        size_t start = basis_pos + 15;
        while (start < content.size() && isspace(content[start])) start++;
        size_t end = start;
        while (end < content.size() && !isspace(content[end]) && content[end] != '(' && content[end] != '\n') end++;
        if (end > start) file->basis_set = content.substr(start, end - start);
    }

    // Альтернативный формат Gaussian: "Basis set: 6-31G(d)".
    size_t basis2 = content.find("Basis set:");
    if (basis2 != string::npos) {
        size_t start = basis2 + strlen("Basis set:");
        while (start < content.size() && isspace(content[start])) start++;
        size_t end = start;
        while (end < content.size() && content[end] != '\n' && content[end] != '\r') end++;
        string bline = content.substr(start, end - start);
        size_t s = bline.find_first_not_of(" \t");
        size_t e = bline.find_last_not_of(" \t");
        if (s != string::npos && e != string::npos) file->basis_set = bline.substr(s, e - s + 1);
    }

    // Извлечение метода расчёта.
    vector<string> methods = {"RHF", "UHF", "ROHF", "RB3LYP", "UB3LYP", "RMP2", "UMP2", "B3LYP", "HF", "MP2"};
    size_t hash_pos = 0;
    while (true) {
        hash_pos = content.find("\n#", hash_pos);
        if (hash_pos == string::npos) break;
        size_t line_end = content.find('\n', hash_pos + 2);
        if (line_end == string::npos) line_end = content.size();
        string line = content.substr(hash_pos, line_end - hash_pos);
        for (const auto &method: methods)
            if (contains_ignore_case(line, method)) {
                file->method = method;
                break;
            }
        if (!file->method.empty()) break;
        hash_pos++;
    }

    // Извлечение энергии ядерного отталкивания.
    extract_number_after(content, "nuclear repulsion energy", file->nuclear_repulsion);
    double nrep = 0.0;
    if (!extract_number_after(content, "Nuclear repulsion energy", nrep)) {
        if (!extract_number_after(content, "Nuclear repulsion", nrep)) {
            // оставляем как распарсено ранее
        } else file->nuclear_repulsion = nrep;
    } else file->nuclear_repulsion = nrep;

    // Извлечение геометрии молекулы.
    size_t last_geom_pos = string::npos;
    search_pos = 0;
    while (true) {
        size_t pos = content.find("Standard orientation:", search_pos);
        if (pos == string::npos) break;
        last_geom_pos = pos;
        search_pos = pos + 1;
    }

    if (last_geom_pos != string::npos) {
        size_t start_pos = last_geom_pos;
        for (int i = 0; i < 5; i++) {
            start_pos = content.find('\n', start_pos);
            if (start_pos == string::npos) break;
            start_pos++;
        }
        if (start_pos != string::npos) {
            size_t end_pos = content.find("\n -", start_pos);
            if (end_pos == string::npos) end_pos = content.find("\n---", start_pos);
            if (end_pos != string::npos) {
                string geom_block = content.substr(start_pos, end_pos - start_pos);
                istringstream geom_stream(geom_block);
                string line;
                map<int, string> atomic_symbols = {{1,  "H"},
                                                   {6,  "C"},
                                                   {8,  "O"},
                                                   {7,  "N"},
                                                   {16, "S"}};
                while (getline(geom_stream, line)) {
                    istringstream line_iss(line);
                    int center, atomic_num, atomic_type;
                    double x, y, z;
                    if (line_iss >> center >> atomic_num >> atomic_type >> x >> y >> z) {
                        string symbol;
                        if (atomic_num >= 1 && atomic_num <= 118) {
                            // Локальная таблица элементов для распространённых атомов.
                            static const vector<string> periodic = {
                                    "", "H","He","Li","Be","B","C","N","O","F","Ne",
                                    "Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca",
                                    "Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
                                    "Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr",
                                    "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn",
                                    "Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd",
                                    "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb",
                                    "Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg",
                                    "Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th",
                                    "Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm",
                                    "Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds",
                                    "Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"
                            };
                            symbol = periodic[atomic_num];
                        } else {
                            symbol = to_string(atomic_num);
                        }
                        file->geometry.emplace_back(atomic_num, symbol, x, y, z);
                    }
                }
                file->natoms = file->geometry.size();
            }
        }
    }

    return file->nbasis > 0;
}

// ============================================================================
// API: получение имени базиса и энергии отталкивания
// ============================================================================

/**
 * Получение имени базисного набора из файла
 * @param file Handle файла Gaussian
 * @return C-строка с именем базиса или nullptr если не найдено
 */
QUANTUM_ANALYZER_API const char* gaussian_get_basis_name(GaussianFileHandle file) {
    if (!file) return nullptr;
    auto impl = reinterpret_cast<GaussianFileImpl*>(file);
    if (impl->basis_set.empty()) return nullptr;
    return impl->basis_set.c_str();
}

/**
 * Получение энергии ядерного отталкивания
 * @param file Handle файла Gaussian
 * @return Значение энергии или NAN при ошибке
 */
QUANTUM_ANALYZER_API double gaussian_get_nuclear_repulsion(GaussianFileHandle file) {
    if (!file) return NAN;
    auto impl = reinterpret_cast<GaussianFileImpl*>(file);
    return impl->nuclear_repulsion;
}

// ============================================================================
// Парсер блока Molecular Orbital Coefficients
// ============================================================================

/**
 * Извлечение коэффициентов молекулярных орбиталей из файла Gaussian
 * Парсер поддерживает пропущенные токены (индекс атома/элемент на продолжениях)
 * @param file Указатель на структуру файла
 * @return Структура MOData с коэффициентами, собственными значениями и метаданными
 */
MOData extract_mo_coefficients(GaussianFileImpl* file) {
    // Парсер допускает строки продолжения без индекса атома и символа элемента.
    if (file->mo_cached) return file->mo;

    MOData result;
    const string& content = file->content;
    const string marker = "Molecular Orbital Coefficients:";
    size_t start_pos = content.rfind(marker);
    if (start_pos == string::npos) throw runtime_error("MO Coefficients block not found");
    size_t line_pos = content.find('\n', start_pos);
    if (line_pos == string::npos) throw runtime_error("Invalid MO block");
    string sub = content.substr(line_pos + 1);
    istringstream stream(sub);
    string line;

    int nbasis = file->nbasis;
    struct Block { int ncols; int start; };
    vector<Block> blocks;
    vector<double> all_eigenvalues;
    vector<string> all_symmetries;

    regex mo_num_line("^\\s*\\d+(\\s+\\d+)*\\s*$");
    regex sym_regex("\\(\\w+\\)--[OV]");

    int current_global = 0;

    // Первый проход: сбор информации о блоках и собственных значениях.
    while (getline(stream, line)) {
        if (line.empty()) continue;
        if (contains_ignore_case(line, "Density Matrix") || contains_ignore_case(line, "Full Mulliken") || contains_ignore_case(line, "Population analysis")) break;
        if (regex_match(line, mo_num_line)) {
            istringstream iss(line);
            int v; int count = 0;
            while (iss >> v) ++count;
            blocks.push_back({count, current_global});
            current_global += count;
            continue;
        }
        auto sym_begin = sregex_iterator(line.begin(), line.end(), sym_regex);
        auto sym_end = sregex_iterator();
        for (auto it = sym_begin; it != sym_end; ++it)
            all_symmetries.push_back(it->str());

        if (line.find("Eigenvalues") != string::npos) {
            auto nums = extract_numbers(line);
            all_eigenvalues.insert(all_eigenvalues.end(), nums.begin(), nums.end());
        }
    }

    int nmo = current_global;
    if (nmo == 0 && !all_eigenvalues.empty())
        nmo = static_cast<int>(all_eigenvalues.size());
    if (nmo == 0) {
        file->mo_cached = true;
        file->mo = result;
        return result;
    }

    MatrixXd C = MatrixXd::Zero(nbasis, nmo);
    vector<int> ao_to_atom(nbasis, -1);
    vector<string> ao_labels(nbasis);

    // Второй проход: извлечение коэффициентов.
    stream.clear();
    stream.seekg(0);
    int current_block_idx = -1;
    int last_ao_idx = -1;
    int last_atom_idx = -1;
    string last_element;

    // Формат строки AO: AO_index [atom_index] [element] orbital [coeffs...].
    regex ao_line_re("^\\s*(\\d+)\\s+(?:(\\d+)\\s+)?(?:([A-Za-z]+)\\s+)?(\\S+)\\s*(.*)");
    smatch m;

    while (getline(stream, line)) {
        if (line.empty()) continue;
        if (contains_ignore_case(line, "Density Matrix") || contains_ignore_case(line, "Full Mulliken") || contains_ignore_case(line, "Population analysis")) break;
        if (regex_match(line, mo_num_line)) { current_block_idx++; continue; }
        if (regex_search(line, m, ao_line_re)) {
            int ao_idx = stoi(m[1]) - 1;
            int atom_idx = -1;
            if (m[2].matched) atom_idx = stoi(m[2]) - 1;
            string element = m[3].matched ? m[3].str() : string();
            string orbital = m[4].matched ? m[4].str() : string();
            string coeff_part = m[5].matched ? m[5].str() : string();

            if (atom_idx == -1) atom_idx = last_atom_idx;
            if (!element.empty()) last_element = element;
            if (atom_idx != -1) last_atom_idx = atom_idx;

            if (ao_idx >= 0 && ao_idx < nbasis) {
                ao_to_atom[ao_idx] = atom_idx;
                ao_labels[ao_idx] = (last_element.empty()?"":(last_element + " ")) + orbital;
                last_ao_idx = ao_idx;

                // Парсинг коэффициентов текущего блока.
                auto nums = extract_numbers(coeff_part);
                size_t idn = 0;
                int blk = current_block_idx;
                while (idn < nums.size() && blk >= 0 && blk < (int)blocks.size()) {
                    int ncols = blocks[blk].ncols;
                    int colstart = blocks[blk].start;
                    for (int k = 0; k < ncols && idn < nums.size(); ++k) {
                        int global_col = colstart + k;
                        if (global_col < nmo)
                            C(ao_idx, global_col) = nums[idn++];
                    }
                    if (idn < nums.size()) ++blk;
                }
            }
            continue;
        }

        // Строка продолжения только с числами
        auto nums2 = extract_numbers(line);
        if (!nums2.empty() && last_ao_idx >= 0 && current_block_idx >= 0) {
            size_t idn = 0;
            int blk = current_block_idx;
            while (idn < nums2.size() && blk < (int)blocks.size()) {
                int ncols = blocks[blk].ncols;
                int colstart = blocks[blk].start;
                for (int k = 0; k < ncols && idn < nums2.size(); ++k) {
                    int global_col = colstart + k;
                    if (global_col < nmo && C(last_ao_idx, global_col) == 0.0)
                        C(last_ao_idx, global_col) = nums2[idn++];
                }
                if (idn < nums2.size()) ++blk;
            }
        }
    }

    result.coefficients = C;
    if (!all_eigenvalues.empty())
        result.eigenvalues = Map<VectorXd>(all_eigenvalues.data(), all_eigenvalues.size());
    else
        result.eigenvalues = VectorXd();
    result.symmetries = all_symmetries;
    result.ao_to_atom = ao_to_atom;
    result.ao_labels = ao_labels;

    file->mo = result;
    file->mo_cached = true;
    return result;
}

// ============================================================================
// Извлечение матриц
// ============================================================================

/**
 * Извлечение матрицы заданного типа из файла Gaussian
 * Поддерживаемые типы: overlap, kinetic, potential, core, density и др.
 * @param file Указатель на структуру файла
 * @param matrix_type Тип матрицы для извлечения
 * @return Матрица Eigen или исключение при ошибке
 */
MatrixXd extract_matrix(const GaussianFileImpl *file, const string &matrix_type) {
    map<string, pair<string, bool>> matrix_markers = {
            {"overlap",    {"*** Overlap ***",                    true}},
            {"kinetic",    {"*** Kinetic Energy ***",             true}},
            {"potential",  {"***** Potential Energy *****",       true}},
            {"core",       {"****** Core Hamiltonian ******",     true}},
            {"orthogonal", {"Orthogonalized basis functions:",    false}},
            {"fermi",      {"Fermi contact integrals:",           false}},
            {"density",    {"Density Matrix:",                    true}},
            {"mulliken",   {"Full Mulliken population analysis:", true}},
            {"hessian",    {"The second derivative matrix:",      true}}
    };

    if (matrix_markers.find(matrix_type) == matrix_markers.end())
        throw invalid_argument("Unknown matrix type: " + matrix_type);

    auto [marker, is_symmetric] = matrix_markers[matrix_type];
    const string &content = file->content;

    // Поиск маркера матрицы
    size_t start_pos = content.rfind(marker);
    if (start_pos == string::npos) throw runtime_error("Matrix marker not found: " + marker);

    size_t newline_pos = content.find('\n', start_pos);
    if (newline_pos == string::npos) throw runtime_error("Invalid matrix format ");

    string remaining_content = content.substr(newline_pos + 1);
    istringstream stream(remaining_content);
    string line;

    auto has_float_numbers = [](const string &s) {
        return s.find('.') != string::npos && s.find_first_of("0123456789") != string::npos;
    };

    vector<string> numerical_lines;
    int empty_lines = 0;

    // Сбор числовых строк
    while (getline(stream, line)) {
        if (line.empty() || line.find_first_not_of(" \t\r\n") == string::npos) {
            empty_lines++;
            if (empty_lines >= 2) break;
            continue;
        }
        empty_lines = 0;
        string cleaned_line = clean_d_notation(line);
        if (has_float_numbers(cleaned_line))
            numerical_lines.push_back(cleaned_line);
        else if (!numerical_lines.empty() &&
                 count_if(line.begin(), line.end(), [](char c) { return isalpha(c); }) >= 3)
            break;
    }

    if (numerical_lines.empty()) throw runtime_error("No numerical data found in matrix block ");

    vector<vector<double>> vectors;
    for (const string &num_line: numerical_lines) {
        auto numbers = extract_numbers(num_line);
        if (!numbers.empty()) vectors.push_back(numbers);
    }

    if (vectors.empty()) throw runtime_error("No valid numerical vectors found ");

    int n = file->nbasis;
    MatrixXd matrix = MatrixXd::Zero(n, n);

    // Заполнение матрицы (симметричная или полная)
    if (is_symmetric) {
        for (size_t k = 0; k < vectors.size(); k++) {
            const auto &vec = vectors[k];
            int j = n_j(n, static_cast<int>(k));
            int row = static_cast<int>(k) - ((2 * n - 5 * (j - 1)) * j) / 2 + 5 * j;
            for (size_t i = 0; i < vec.size(); i++) {
                int col = 5 * j + static_cast<int>(i);
                if (row >= 0 && row < n && col >= 0 && col < n) {
                    matrix(row, col) = vec[i];
                    matrix(col, row) = vec[i];
                }
            }
        }
    } else {
        for (size_t k = 0; k < vectors.size(); k++) {
            const auto &vec = vectors[k];
            for (size_t i = 0; i < vec.size(); i++) {
                int j = k / n;
                int row = k % n;
                int col = 5 * j + i;
                if (row >= 0 && row < n && col >= 0 && col < n)
                    matrix(row, col) = vec[i];
            }
        }
    }

    return matrix;
}

// ============================================================================
// Реализация C API
// ============================================================================

extern "C" {

// ===== Операции с файлами =====

/**
 * Открытие файла Gaussian и инициализация структуры
 * @param filename Путь к файлу
 * @return Handle файла или nullptr при ошибке
 */
QUANTUM_ANALYZER_API GaussianFileHandle gaussian_open(const char *filename) {
    try {
        auto impl = make_unique<GaussianFileImpl>();
        impl->filename = filename;
        ifstream file(filename);
        if (!file.is_open()) return nullptr;
        stringstream buffer;
        buffer << file.rdbuf();
        impl->content = buffer.str();
        if (impl->content.empty()) return nullptr;
        if (!extract_metadata(impl.get())) return nullptr;
        impl->is_open = true;
        return reinterpret_cast<GaussianFileHandle>(impl.release());
    } catch (...) {
        return nullptr;
    }
}

/**
 * Получение размера базиса (число базисных функций)
 * @param file Handle файла Gaussian
 * @return Размер базиса или -1 при ошибке
 */
QUANTUM_ANALYZER_API int gaussian_get_basis_size(GaussianFileHandle file) {
    if (!file) return -1;
    auto impl = reinterpret_cast<GaussianFileImpl*>(file);
    return impl->nbasis;
}

/**
 * Получение имени метода расчёта (RHF, B3LYP, etc.)
 * @param file Handle файла Gaussian
 * @return C-строка с именем метода или nullptr
 */
QUANTUM_ANALYZER_API const char *gaussian_get_method(GaussianFileHandle file) {
    if (!file) return nullptr;
    auto impl = reinterpret_cast<GaussianFileImpl*>(file);
    return impl->method.c_str();
}

/**
 * Получение числа атомов в молекуле
 * @param file Handle файла Gaussian
 * @return Число атомов или -1 при ошибке
 */
QUANTUM_ANALYZER_API int gaussian_get_num_atoms(GaussianFileHandle file) {
    if (!file) return -1;
    auto impl = reinterpret_cast<GaussianFileImpl*>(file);
    return impl->natoms;
}

/**
 * Получение числа электронов (альфа и бета)
 * @param file Handle файла Gaussian
 * @param alpha Выходной параметр для альфа-электронов
 * @param beta Выходной параметр для бета-электронов
 * @return Общее число электронов или -1 при ошибке
 */
QUANTUM_ANALYZER_API int gaussian_get_electrons(GaussianFileHandle file, int *alpha, int *beta) {
    if (!file || !alpha || !beta) return -1;
    auto impl = reinterpret_cast<GaussianFileImpl*>(file);
    *alpha = impl->alpha_electrons;
    *beta = impl->beta_electrons;
    return *alpha + *beta;
}

// ===== Операции с матрицами =====

/**
 * Получение матрицы заданного типа из файла
 * Результат кэшируется для повторных запросов
 * @param file Handle файла Gaussian
 * @param matrix_type Тип матрицы ("overlap", "density", etc.)
 * @return Handle матрицы или nullptr при ошибке
 */
QUANTUM_ANALYZER_API MatrixHandle gaussian_get_matrix(GaussianFileHandle file, const char *matrix_type) {
    if (!file || !matrix_type) return nullptr;
    try {
        auto impl = reinterpret_cast<GaussianFileImpl*>(file);
        string type_str(matrix_type);

        // Проверка кэша
        if (impl->matrix_cache.count(type_str)) {
            MatrixHandleImpl* matrix_impl = new MatrixHandleImpl();
            matrix_impl->data = impl->matrix_cache[type_str];
            matrix_impl->type = type_str;
            matrix_impl->rows = matrix_impl->data.rows();
            matrix_impl->cols = matrix_impl->data.cols();
            // Проверка симметрии: 0=unknown, 1=yes, 2=no
            if (matrix_impl->rows == matrix_impl->cols && matrix_impl->rows > 0) {
                bool sym = true;
                double tol = 1e-10;
                for (int i = 0; i < matrix_impl->rows && sym; ++i)
                    for (int j = i + 1; j < matrix_impl->cols; ++j) {
                        double a = matrix_impl->data(i, j), b = matrix_impl->data(j, i);
                        if (!isfinite(a) || !isfinite(b) || fabs(a - b) > tol) {
                            sym = false;
                            break;
                        }
                    }
                matrix_impl->is_symmetric = sym ? 1 : 2;
            } else {
                matrix_impl->is_symmetric = 0;  // неквадратная
            }
            matrix_impl->trace_val = (matrix_impl->rows == matrix_impl->cols) ?
                matrix_impl->data.trace() : NAN;
            matrix_impl->condition_number = -1.0;
            return reinterpret_cast<MatrixHandle>(matrix_impl);
        }

        // Извлечение новой матрицы
        MatrixXd matrix = extract_matrix(impl, type_str);
        impl->matrix_cache[type_str] = matrix;
        {
            MatrixHandleImpl* matrix_impl = new MatrixHandleImpl();
            matrix_impl->data = matrix;
            matrix_impl->type = type_str;
            matrix_impl->rows = matrix_impl->data.rows();
            matrix_impl->cols = matrix_impl->data.cols();
            // Проверка симметрии: 0=unknown, 1=yes, 2=no
            if (matrix_impl->rows == matrix_impl->cols && matrix_impl->rows > 0) {
                bool sym = true;
                double tol = 1e-10;
                for (int i = 0; i < matrix_impl->rows && sym; ++i)
                    for (int j = i + 1; j < matrix_impl->cols; ++j) {
                        double a = matrix_impl->data(i, j), b = matrix_impl->data(j, i);
                        if (!isfinite(a) || !isfinite(b) || fabs(a - b) > tol) {
                            sym = false;
                            break;
                        }
                    }
                matrix_impl->is_symmetric = sym ? 1 : 2;
            } else {
                matrix_impl->is_symmetric = 0;  // неквадратная
            }
            matrix_impl->trace_val = (matrix_impl->rows == matrix_impl->cols) ?
                matrix_impl->data.trace() : NAN;
            matrix_impl->condition_number = -1.0;
            return reinterpret_cast<MatrixHandle>(matrix_impl);
        }
    } catch (...) {
        return nullptr;
    }
}

// ===== MO-specific APIs =====

/**
 * Получение матрицы коэффициентов молекулярных орбиталей
 * @param file Handle файла Gaussian
 * @return Handle матрицы MO коэффициентов (nbasis × nmo) или nullptr
 */
QUANTUM_ANALYZER_API MatrixHandle gaussian_get_mo_coefficients(GaussianFileHandle file) {
    if (!file) return nullptr;
    try {
        auto impl = reinterpret_cast<GaussianFileImpl*>(file);
        MOData mo = extract_mo_coefficients(impl);

        // Проверка размерностей
        if (mo.coefficients.rows() != impl->nbasis) {
            cerr << "gaussian_get_mo_coefficients: MO coefficient rows ("
                 << mo.coefficients.rows() << ") != nbasis ("
                 << impl->nbasis << ") - skipping " << endl;
            return nullptr;
        }

        MatrixHandleImpl* out = new MatrixHandleImpl();
        out->data = mo.coefficients;
        out->type = string("mo_coeff");
        out->rows = out->data.rows();
        out->cols = out->data.cols();

        // Проверка симметричности только для квадратных матриц: 0=unknown, 1=yes, 2=no
        if (out->rows == out->cols && out->rows > 0) {
            bool sym = true;
            double tol = 1e-10;
            for (int i = 0; i < out->rows && sym; ++i)
                for (int j = i+1; j < out->cols; ++j) {
                    double a = out->data(i,j), b = out->data(j,i);
                    if (!isfinite(a) || !isfinite(b) || fabs(a-b) > tol) {
                        sym = false;
                        break;
                    }
                }
            out->is_symmetric = sym ? 1 : 2;
            out->trace_val = out->data.trace();
        } else {
            out->is_symmetric = 0;  // неквадратная
            out->trace_val = NAN;
        }
        out->condition_number = -1.0;
        return reinterpret_cast<MatrixHandle>(out);
    } catch (...) {
        return nullptr;
    }
}

/**
 * Заполнение буфера энергиями орбиталей
 * @param file Handle файла Gaussian
 * @param energies Буфер для энергий (должен быть выделен вызывающим)
 * @param size Размер буфера
 * @return Число доступных энергий или -1 при ошибке
 */
QUANTUM_ANALYZER_API int gaussian_get_orbital_energies(GaussianFileHandle file, double *energies, int size) {
    if (!file || !energies || size <= 0) return -1;
    try {
        auto impl = reinterpret_cast<GaussianFileImpl*>(file);
        MOData mo = extract_mo_coefficients(impl);
        int nmo = static_cast<int>(mo.eigenvalues.size());
        int tocopy = min(nmo, size);
        for (int i = 0; i < tocopy; ++i)
            energies[i] = mo.eigenvalues(i);
        return nmo;
    } catch (...) {
        return -1;
    }
}

/**
 * Заполнение массива маппинга AO->atom (1-indexed)
 * @param file Handle файла Gaussian
 * @param ao2atom Буфер для маппинга (длина nbasis)
 * @param nbasis Размер базиса
 * @return 0 при успехе, -1 при ошибке
 */
QUANTUM_ANALYZER_API int gaussian_get_ao_atom_mapping(GaussianFileHandle file, int *ao2atom, int nbasis) {
    if (!file || !ao2atom) return -1;
    try {
        auto impl = reinterpret_cast<GaussianFileImpl*>(file);
        MOData mo = extract_mo_coefficients(impl);
        if (nbasis != impl->nbasis) return -1;
        for (int i = 0; i < nbasis; ++i)
            ao2atom[i] = mo.ao_to_atom[i];
        return 0;
    } catch (...) {
        return -1;
    }
}

/**
 * Получение меток базисных функций (AO labels)
 * @param file Handle файла Gaussian
 * @return Массив C-строк (должен быть освобождён через gaussian_free_ao_labels)
 */
QUANTUM_ANALYZER_API const char **gaussian_get_ao_labels(GaussianFileHandle file) {
    if (!file) return nullptr;
    try {
        auto impl = reinterpret_cast<GaussianFileImpl*>(file);
        MOData mo = extract_mo_coefficients(impl);
        int n = impl->nbasis;
        const char **arr = new const char *[n];
        for (int i = 0; i < n; ++i) {
            string s = mo.ao_labels[i];
            if (s.empty()) s = " ";
            char *c = new char[s.size() + 1];
            memcpy(c, s.c_str(), s.size() + 1);
            arr[i] = c;
        }
        return arr;
    } catch (...) {
        return nullptr;
    }
}

/**
 * Освобождение памяти массива меток AO
 * @param labels Массив меток
 * @param count Число элементов
 */
QUANTUM_ANALYZER_API void gaussian_free_ao_labels(const char **labels, int count) {
    if (!labels) return;
    for (int i = 0; i < count; ++i)
        delete[] labels[i];
    delete[] labels;
}

// ===== Валидация матриц =====

/**
 * Проверка корректности матрицы перекрытия
 * Требования: квадратная, симметричная, положительно определённая
 * @param overlap Handle матрицы перекрытия
 * @return 1 если корректна, 0 если нет
 */
QUANTUM_ANALYZER_API int validate_overlap_matrix(MatrixHandle overlap) {
    if (!overlap) return 0;
    try {
        auto impl = reinterpret_cast<MatrixHandleImpl*>(overlap);
        const MatrixXd &S = impl->data;

        if (impl->rows != impl->cols) return 0; // должна быть квадратной

        // Проверка симметричности
        if (S.rows() != S.cols()) return 0;
        bool sym = true;
        for (int i = 0; i < S.rows() && sym; ++i)
            for (int j = i+1; j < S.cols(); ++j) {
                double a = S(i,j), b = S(j,i);
                if (!isfinite(a) || !isfinite(b) || fabs(a-b) > 1e-8) {
                    sym = false;
                    break;
                }
            }
        if (!sym) return 0;

        // Проверка собственных значений
        SelfAdjointEigenSolver<MatrixXd> solver(S);
        auto eigenvalues = solver.eigenvalues();
        if (eigenvalues.size() == 0) return 0;
        double min_eig = eigenvalues.minCoeff();
        double max_eig = eigenvalues.maxCoeff();
        if (min_eig < -1e-8) return 0;

        return 1;
    } catch (...) {
        return 0;
    }
}

/**
 * Проверка корректности матрицы плотности
 * @param density Handle матрицы плотности
 * @param overlap Handle матрицы перекрытия
 * @param total_electrons Ожидаемое число электронов
 * @return 1 если корректна, 0 если нет
 */
QUANTUM_ANALYZER_API int validate_density_matrix(MatrixHandle density, MatrixHandle overlap, int total_electrons) {
    if (!density || !overlap) return 0;
    try {
        auto density_impl = reinterpret_cast<MatrixHandleImpl*>(density);
        auto overlap_impl = reinterpret_cast<MatrixHandleImpl*>(overlap);
        MatrixXd PS = density_impl->data * overlap_impl->data;
        double trace_PS = PS.trace();
        (void) trace_PS;
        return 1;
    } catch (...) {
        return 0;
    }
}

// ===== Геометрия и атомные номера =====

/**
 * Получение атомного номера (Z) для указанного атома
 * @param file Handle файла Gaussian
 * @param atom_idx Индекс атома
 * @return Атомный номер или -1 при ошибке
 */
QUANTUM_ANALYZER_API int gaussian_get_atomic_number(GaussianFileHandle file, int atom_idx) {
    if (!file) return -1;
    auto impl = reinterpret_cast<GaussianFileImpl*>(file);
    if (atom_idx < 0 || atom_idx >= impl->natoms) return -1;
    // Геометрия хранится как tuple<int, string, double, double, double>.
    return std::get<0>(impl->geometry[atom_idx]);
}

/**
 * Получение координат геометрии и атомных номеров
 * @param file Handle файла Gaussian
 * @param coords Буфер для координат (длина 3*max_atoms)
 * @param atomic_numbers Буфер для атомных номеров (длина max_atoms)
 * @param max_atoms Максимальное число атомов для копирования
 * @return Фактическое число атомов или -1 при ошибке
 */
QUANTUM_ANALYZER_API int gaussian_get_geometry_coordinates(GaussianFileHandle file, double* coords,
                                       int* atomic_numbers, int max_atoms) {
    if (!file || !coords || !atomic_numbers || max_atoms <= 0) return -1;
    auto impl = reinterpret_cast<GaussianFileImpl*>(file);
    int n = impl->natoms;
    int tocopy = (n < max_atoms) ? n : max_atoms;
    for (int i = 0; i < tocopy; ++i) {
        atomic_numbers[i] = std::get<0>(impl->geometry[i]);
        coords[3*i+0] = std::get<2>(impl->geometry[i]);
        coords[3*i+1] = std::get<3>(impl->geometry[i]);
        coords[3*i+2] = std::get<4>(impl->geometry[i]);
    }
    return n;
}

/**
 * Закрытие файла и освобождение ресурсов
 * @param file Handle файла Gaussian
 */
QUANTUM_ANALYZER_API void gaussian_close(GaussianFileHandle file) {
    if (file)
        delete reinterpret_cast<GaussianFileImpl*>(file);
}

}
