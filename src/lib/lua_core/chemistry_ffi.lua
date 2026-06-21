-- LuaJIT FFI-обёртка над C API.
-- Molecule является корневым владельцем C-ресурсов и ведёт реестр refcount/LRU.
-- Matrix и Group являются лёгкими дескрипторами и не владеют памятью напрямую.

local ffi = require("ffi")

-- ============================================================================
-- Объявления LuaJIT FFI.
-- ============================================================================

local is_windows = package.config:sub(1,1) == "\\"
local libname = is_windows and "quantum_analyzer_core" or "libquantum_analyzer_core"

local gaussian = QUANTUM_ANALYZER_EMBEDDED and ffi.C or nil
if not gaussian then
    local lib_path = package.searchpath(libname, package.cpath)
    if not lib_path then
        error("Library not found: " .. libname)
    end
    gaussian = ffi.load(lib_path)
end

ffi.cdef[[
    typedef void* GaussianFileHandle;
    typedef void* MatrixHandle;
    typedef void* GroupHandle;

    GaussianFileHandle gaussian_open(const char* filename);
    void gaussian_close(GaussianFileHandle file);
    int gaussian_get_basis_size(GaussianFileHandle file);
    int gaussian_get_num_atoms(GaussianFileHandle file);
    int gaussian_get_atomic_number(GaussianFileHandle file, int atom_idx);
    const char* gaussian_get_basis_name(GaussianFileHandle file);
    double gaussian_get_nuclear_repulsion(GaussianFileHandle file);
    int gaussian_get_ao_atom_mapping(GaussianFileHandle file, int* ao2atom, int nbasis);
    MatrixHandle gaussian_get_mo_coefficients(GaussianFileHandle file);
    int gaussian_get_orbital_energies(GaussianFileHandle file, double* energies, int size);

    MatrixHandle matrix_create(int rows, int cols, const double* data);
    MatrixHandle gaussian_get_matrix(GaussianFileHandle file, const char* matrix_type);
    int matrix_get_size(MatrixHandle matrix, int* rows, int* cols);
    double matrix_get_element(MatrixHandle matrix, int row, int col);
    double matrix_trace(MatrixHandle matrix);
    int matrix_is_symmetric(MatrixHandle matrix);
    double matrix_condition_number(MatrixHandle matrix);
    void matrix_free(MatrixHandle matrix);

    int matrix_eigenvalues(MatrixHandle matrix, double* eigenvalues);
    int matrix_eigensystem(MatrixHandle matrix, double* eigenvalues, MatrixHandle* eigenvectors);
    int matrix_svd(MatrixHandle matrix, MatrixHandle* U, MatrixHandle* S, MatrixHandle* Vt);
    MatrixHandle matrix_power(MatrixHandle matrix, double exponent);
    MatrixHandle matrix_multiply(MatrixHandle a, MatrixHandle b);
    MatrixHandle matrix_ao_to_mo(MatrixHandle matrix_ao, MatrixHandle mo_coeff);
    MatrixHandle matrix_add(MatrixHandle a, double alpha, MatrixHandle b, double beta);
    MatrixHandle matrix_scale(MatrixHandle a, double alpha);
    MatrixHandle matrix_transpose(MatrixHandle a);
    MatrixHandle matrix_triple_product_symm(MatrixHandle a, MatrixHandle b);
    MatrixHandle matrix_hadamard(MatrixHandle a, MatrixHandle b);
    MatrixHandle matrix_cwise_pow(MatrixHandle a, double exponent);
    MatrixHandle matrix_threshold(MatrixHandle a, double theta);
    MatrixHandle matrix_clamp(MatrixHandle a, double lo, double hi);
    double matrix_norm_fro(MatrixHandle a);
    double matrix_max_abs(MatrixHandle a);
    double matrix_min_abs_nonzero(MatrixHandle a);
    MatrixHandle matrix_get_diagonal(MatrixHandle a);
    MatrixHandle matrix_extract_block(MatrixHandle m, GroupHandle rows, GroupHandle cols);
    double matrix_block_trace(MatrixHandle m, GroupHandle g);
    double matrix_block_sum_squares(MatrixHandle m, GroupHandle rows, GroupHandle cols);
    double matrix_block_mayer_pair(MatrixHandle p, MatrixHandle s, GroupHandle rows, GroupHandle cols);
    double matrix_block_norm_fro(MatrixHandle m, GroupHandle rows, GroupHandle cols);
    MatrixHandle matrix_symm_pow(MatrixHandle a, double exponent, double eps);

    GroupHandle group_create(int nbasis);
    GroupHandle group_create_full(int nbasis);
    GroupHandle group_from_indices(const int* idx, int count, int nbasis);
    GroupHandle group_from_atom(GaussianFileHandle file, int atom_idx);
    void group_set_bit(GroupHandle g, int index, int value);
    int group_get_bit(GroupHandle g, int index);
    int group_count(GroupHandle g);
    int group_nbasis(GroupHandle g);
    GroupHandle group_or(GroupHandle a, GroupHandle b);
    GroupHandle group_and(GroupHandle a, GroupHandle b);
    void group_free(GroupHandle g);
]]

-- ============================================================================
-- Molecule: корневой владелец C-ресурсов и реестр дескрипторов.
-- ============================================================================

local Molecule = {}
Molecule.__index = Molecule

local DEFAULT_LIMITS = {
    matrix = 16,
    group  = 100
}

function Molecule.load(filename)
    local handle = gaussian.gaussian_open(filename)
    if handle == nil then
        error("Cannot open Gaussian file: " .. filename)
    end

    local self = setmetatable({
        handle = handle,
        filename = filename,
        _closed = false,
        _registry = {
            matrix = {},
            group  = {},
            next_id  = 1,
            limits   = {
                matrix = DEFAULT_LIMITS.matrix,
                group  = DEFAULT_LIMITS.group
            },
            -- Кэш хранит только ID persistent-ресурсов, а не Lua-wrapper'ы.
            persistent_ids = {
                overlap = nil,
                density = nil
            }
        }
    }, Molecule)

    return self
end

-- Внутренняя регистрация C-ресурса.
function Molecule:_register(kind, handle, opts)
    opts = opts or {}
    local persistent = opts.persistent or false
    local m_type = opts.type or "unknown"

    local reg = self._registry
    local pool = reg[kind]
    local limit = reg.limits[kind]

    -- Эвикция перед добавлением применяется только к temporary-ресурсам.
    if not persistent then
        self:_evict_if_needed(kind)
    end

    local id = reg.next_id
    reg.next_id = id + 1

    pool[id] = {
        handle     = handle,
        refcount   = 0,
        last_used  = os.clock(),
        persistent = persistent,
        type       = m_type,
        freed      = false
    }

    return id
end

-- Эвикция по LRU выполняется до попадания в заданный лимит.
function Molecule:_evict_if_needed(kind)
    local reg = self._registry
    local pool = reg[kind]
    local limit = reg.limits[kind]

    while true do
        local temp_count = 0
        local oldest_id = nil
        local oldest_time = math.huge

        for id, entry in pairs(pool) do
            if not entry.persistent then
                temp_count = temp_count + 1
                if entry.refcount <= 0 and not entry.freed then
                    if entry.last_used < oldest_time then
                        oldest_time = entry.last_used
                        oldest_id = id
                    end
                end
            end
        end

        -- Эвикция вызывается перед добавлением нового ресурса, поэтому при temp_count == limit
        -- нужно освободить один уже отпущенный дескриптор.
        if temp_count < limit or not oldest_id then
            break
        end

        -- Удаляем самый старый ресурс без активных ссылок.
        local entry = pool[oldest_id]
        if kind == "matrix" then
            gaussian.matrix_free(entry.handle)
        elseif kind == "group" then
            gaussian.group_free(entry.handle)
        end
        entry.freed = true
        entry.handle = nil
        pool[oldest_id] = nil
    end
end

-- Возвращает запись реестра с проверкой жизненного цикла.
function Molecule:_get_entry(kind, id)
    if self._closed then
        error("Molecule is closed")
    end
    local pool = self._registry[kind]
    local entry = pool and pool[id]
    if not entry or entry.freed then
        error(kind:sub(1,1):upper() .. kind:sub(2) .. " resource is invalid or freed")
    end
    return entry
end

-- Обновляет время доступа и возвращает C-handle.
function Molecule:_touch(kind, id)
    local entry = self:_get_entry(kind, id)
    entry.last_used = os.clock()
    return entry.handle
end

-- Уменьшает refcount Lua-дескриптора.
function Molecule:_release(kind, id)
    if self._closed then return end
    local pool = self._registry[kind]
    local entry = pool and pool[id]
    if not entry then return end

    if entry.refcount > 0 then
        entry.refcount = entry.refcount - 1
    end
    -- Эвикция выполняется при регистрации новых ресурсов.
end

-- Закрывает молекулу и освобождает все C-ресурсы, которыми она владеет.
function Molecule:close()
    if self._closed then return end

    local reg = self._registry

    for id, entry in pairs(reg.matrix) do
        if not entry.freed and entry.handle then
            gaussian.matrix_free(entry.handle)
            entry.freed = true
            entry.handle = nil
        end
    end

    for id, entry in pairs(reg.group) do
        if not entry.freed and entry.handle then
            gaussian.group_free(entry.handle)
            entry.freed = true
            entry.handle = nil
        end
    end

    gaussian.gaussian_close(self.handle)
    self.handle = nil
    self._closed = true
    reg.matrix = {}
    reg.group = {}
    reg.persistent_ids = {}
end

function Molecule:is_closed()
    return self._closed
end

-- ============================================================================
-- Matrix: Lua-дескриптор матрицы без прямого владения памятью.
-- ============================================================================

local Matrix = {}
Matrix.__index = Matrix
local Group

-- Создаёт Matrix-дескриптор и увеличивает refcount ресурса.
function Matrix._from_id(molecule, id, m_type)
    local entry = molecule:_get_entry("matrix", id)
    entry.refcount = entry.refcount + 1

    local self = setmetatable({
        molecule = molecule,
        _id = id,
        _type = m_type,
        _released = false
    }, Matrix)
    return self
end

local function check_matrix(self, method)
    if self._released then
        error("Matrix:" .. method .. " called on released matrix")
    end
    self.molecule:_touch("matrix", self._id)
end

local function get_handle(self)
    return self.molecule:_touch("matrix", self._id)
end

function Matrix.from_table(molecule, rows, cols, values)
    if molecule:is_closed() then error("Molecule is closed") end
    if type(rows) ~= "number" or type(cols) ~= "number" or rows <= 0 or cols <= 0 then
        error("Matrix.from_table expects positive rows and cols")
    end
    if type(values) ~= "table" or #values ~= rows * cols then
        error("Matrix.from_table expects rows*cols values")
    end

    local arr = ffi.new("double[?]", rows * cols)
    for i = 1, rows * cols do
        local v = tonumber(values[i])
        if not v then error("Matrix.from_table expects numeric values") end
        arr[i - 1] = v
    end

    local h = gaussian.matrix_create(rows, cols, arr)
    if h == nil then error("matrix_create failed") end
    local id = molecule:_register("matrix", h, {type = "manual"})
    return Matrix._from_id(molecule, id, "manual")
end

function Matrix:size()
    check_matrix(self, "size")
    local h = get_handle(self)
    local rows = ffi.new("int[1]")
    local cols = ffi.new("int[1]")
    gaussian.matrix_get_size(h, rows, cols)
    return rows[0], cols[0]
end

function Matrix:get(i, j)
    check_matrix(self, "get")
    return gaussian.matrix_get_element(get_handle(self), i - 1, j - 1)
end

function Matrix:trace()
    check_matrix(self, "trace")
    return gaussian.matrix_trace(get_handle(self))
end

function Matrix:is_symmetric()
    check_matrix(self, "is_symmetric")
    return gaussian.matrix_is_symmetric(get_handle(self)) == 1
end

function Matrix:condition_number()
    check_matrix(self, "condition_number")
    return gaussian.matrix_condition_number(get_handle(self))
end

function Matrix:pow(exponent)
    check_matrix(self, "pow")
    local h_res = gaussian.matrix_power(get_handle(self), exponent)
    if h_res == nil then error("matrix_power failed") end
    local id = self.molecule:_register("matrix", h_res, {type = self._type .. "^" .. exponent})
    return Matrix._from_id(self.molecule, id, self._type .. "^" .. exponent)
end

function Matrix:symm_pow(exponent, eps)
    check_matrix(self, "symm_pow")
    eps = eps or 1e-12
    local h_res = gaussian.matrix_symm_pow(get_handle(self), exponent, eps)
    if h_res == nil then error("matrix_symm_pow failed") end
    local id = self.molecule:_register("matrix", h_res, {type = self._type .. "^" .. exponent})
    return Matrix._from_id(self.molecule, id, self._type .. "^" .. exponent)
end

function Matrix:mul(other)
    check_matrix(self, "mul")
    if getmetatable(other) ~= Matrix then error("Matrix:mul expects Matrix") end
    check_matrix(other, "mul")
    local h_res = gaussian.matrix_multiply(get_handle(self), get_handle(other))
    if h_res == nil then error("matrix_multiply failed") end
    local id = self.molecule:_register("matrix", h_res, {type = self._type .. "*" .. other._type})
    return Matrix._from_id(self.molecule, id, self._type .. "*" .. other._type)
end

function Matrix:ao_to_mo(mo_coeff)
    check_matrix(self, "ao_to_mo")
    if getmetatable(mo_coeff) ~= Matrix then error("Matrix:ao_to_mo expects Matrix") end
    check_matrix(mo_coeff, "ao_to_mo")
    local h_res = gaussian.matrix_ao_to_mo(get_handle(self), get_handle(mo_coeff))
    if h_res == nil then error("matrix_ao_to_mo failed") end
    local id = self.molecule:_register("matrix", h_res, {type = self._type .. "_MO"})
    return Matrix._from_id(self.molecule, id, self._type .. "_MO")
end

function Matrix:add(other, alpha, beta)
    check_matrix(self, "add")
    if getmetatable(other) ~= Matrix then error("Matrix:add expects Matrix") end
    check_matrix(other, "add")
    alpha = alpha or 1.0
    beta = beta or 1.0
    local h_res = gaussian.matrix_add(get_handle(self), alpha, get_handle(other), beta)
    if h_res == nil then error("matrix_add failed") end
    local id = self.molecule:_register("matrix", h_res, {type = self._type .. "+" .. other._type})
    return Matrix._from_id(self.molecule, id, self._type .. "+" .. other._type)
end

function Matrix:scale(alpha)
    check_matrix(self, "scale")
    alpha = alpha or 1.0
    local h_res = gaussian.matrix_scale(get_handle(self), alpha)
    if h_res == nil then error("matrix_scale failed") end
    local id = self.molecule:_register("matrix", h_res, {type = self._type .. "*s"})
    return Matrix._from_id(self.molecule, id, self._type .. "*s")
end

function Matrix:transpose()
    check_matrix(self, "transpose")
    local h_res = gaussian.matrix_transpose(get_handle(self))
    if h_res == nil then error("matrix_transpose failed") end
    local id = self.molecule:_register("matrix", h_res, {type = self._type .. "T"})
    return Matrix._from_id(self.molecule, id, self._type .. "T")
end

function Matrix:eigenvalues()
    check_matrix(self, "eigenvalues")
    local rows, cols = self:size()
    if rows ~= cols then error("eigenvalues requires square matrix") end
    local values = ffi.new("double[?]", rows)
    local ret = gaussian.matrix_eigenvalues(get_handle(self), values)
    if ret ~= 0 then error("matrix_eigenvalues failed") end
    local out = {}
    for i = 0, rows - 1 do out[i + 1] = values[i] end
    return out
end

function Matrix:eigensystem()
    check_matrix(self, "eigensystem")
    local rows, cols = self:size()
    if rows ~= cols then error("eigensystem requires square matrix") end
    local n = rows
    local eigenvalues = ffi.new("double[?]", n)
    local eigenvectors_handle = ffi.new("MatrixHandle[1]")
    local ret = gaussian.matrix_eigensystem(get_handle(self), eigenvalues, eigenvectors_handle)
    if ret ~= 0 then error("matrix_eigensystem failed") end
    local evals = {}
    for i = 0, n - 1 do evals[i + 1] = eigenvalues[i] end
    local id = self.molecule:_register("matrix", eigenvectors_handle[0], {type = "eigenvectors"})
    return evals, Matrix._from_id(self.molecule, id, "eigenvectors")
end

function Matrix:svd()
    check_matrix(self, "svd")
    local U_h = ffi.new("MatrixHandle[1]")
    local S_h = ffi.new("MatrixHandle[1]")
    local Vt_h = ffi.new("MatrixHandle[1]")
    local ret = gaussian.matrix_svd(get_handle(self), U_h, S_h, Vt_h)
    if ret ~= 0 then error("matrix_svd failed") end
    local id_U = self.molecule:_register("matrix", U_h[0], {type = "svd_U"})
    local id_S = self.molecule:_register("matrix", S_h[0], {type = "svd_S"})
    local id_Vt = self.molecule:_register("matrix", Vt_h[0], {type = "svd_Vt"})
    return Matrix._from_id(self.molecule, id_U, "svd_U"),
           Matrix._from_id(self.molecule, id_S, "svd_S"),
           Matrix._from_id(self.molecule, id_Vt, "svd_Vt")
end

function Matrix:triple_product_symm(other)
    check_matrix(self, "triple_product_symm")
    if getmetatable(other) ~= Matrix then error("Matrix:triple_product_symm expects Matrix") end
    check_matrix(other, "triple_product_symm")
    local h_res = gaussian.matrix_triple_product_symm(get_handle(self), get_handle(other))
    if h_res == nil then error("matrix_triple_product_symm failed") end
    local id = self.molecule:_register("matrix", h_res, {type = "triple_product_symm"})
    return Matrix._from_id(self.molecule, id, "triple_product_symm")
end

function Matrix:hadamard(other)
    check_matrix(self, "hadamard")
    if getmetatable(other) ~= Matrix then error("Matrix:hadamard expects Matrix") end
    check_matrix(other, "hadamard")
    local h_res = gaussian.matrix_hadamard(get_handle(self), get_handle(other))
    if h_res == nil then error("matrix_hadamard failed") end
    local id = self.molecule:_register("matrix", h_res, {type = "hadamard"})
    return Matrix._from_id(self.molecule, id, "hadamard")
end

function Matrix:cwise_pow(exponent)
    check_matrix(self, "cwise_pow")
    local h_res = gaussian.matrix_cwise_pow(get_handle(self), exponent)
    if h_res == nil then error("matrix_cwise_pow failed") end
    local id = self.molecule:_register("matrix", h_res, {type = "cwise_pow"})
    return Matrix._from_id(self.molecule, id, "cwise_pow")
end

function Matrix:threshold(theta)
    check_matrix(self, "threshold")
    local h_res = gaussian.matrix_threshold(get_handle(self), theta)
    if h_res == nil then error("matrix_threshold failed") end
    local id = self.molecule:_register("matrix", h_res, {type = "threshold"})
    return Matrix._from_id(self.molecule, id, "threshold")
end

function Matrix:clamp(lo, hi)
    check_matrix(self, "clamp")
    local h_res = gaussian.matrix_clamp(get_handle(self), lo, hi)
    if h_res == nil then error("matrix_clamp failed") end
    local id = self.molecule:_register("matrix", h_res, {type = "clamp"})
    return Matrix._from_id(self.molecule, id, "clamp")
end

function Matrix:norm_fro()
    check_matrix(self, "norm_fro")
    return gaussian.matrix_norm_fro(get_handle(self))
end

function Matrix:max_abs()
    check_matrix(self, "max_abs")
    return gaussian.matrix_max_abs(get_handle(self))
end

function Matrix:min_abs_nonzero()
    check_matrix(self, "min_abs_nonzero")
    return gaussian.matrix_min_abs_nonzero(get_handle(self))
end

function Matrix:diagonal()
    check_matrix(self, "diagonal")
    local h_res = gaussian.matrix_get_diagonal(get_handle(self))
    if h_res == nil then error("matrix_get_diagonal failed") end
    local id = self.molecule:_register("matrix", h_res, {type = "diagonal"})
    return Matrix._from_id(self.molecule, id, "diagonal")
end

function Matrix:extract_block(rows_group, cols_group)
    check_matrix(self, "extract_block")
    if getmetatable(rows_group) ~= Group or getmetatable(cols_group) ~= Group then
        error("Matrix:extract_block expects Group arguments")
    end
    rows_group:_check("extract_block")
    cols_group:_check("extract_block")
    local h_res = gaussian.matrix_extract_block(get_handle(self), rows_group:_get_handle(), cols_group:_get_handle())
    if h_res == nil then error("matrix_extract_block failed") end
    local id = self.molecule:_register("matrix", h_res, {type = "block"})
    return Matrix._from_id(self.molecule, id, "block")
end

function Matrix:block_trace(group)
    check_matrix(self, "block_trace")
    if getmetatable(group) ~= Group then error("Matrix:block_trace expects Group") end
    group:_check("block_trace")
    return gaussian.matrix_block_trace(get_handle(self), group:_get_handle())
end

function Matrix:block_sum_squares(rows_group, cols_group)
    check_matrix(self, "block_sum_squares")
    if getmetatable(rows_group) ~= Group or getmetatable(cols_group) ~= Group then
        error("Matrix:block_sum_squares expects Group arguments")
    end
    rows_group:_check("block_sum_squares")
    cols_group:_check("block_sum_squares")
    return gaussian.matrix_block_sum_squares(get_handle(self), rows_group:_get_handle(), cols_group:_get_handle())
end

function Matrix:block_mayer_pair(other, rows_group, cols_group)
    check_matrix(self, "block_mayer_pair")
    if getmetatable(other) ~= Matrix then error("Matrix:block_mayer_pair expects Matrix") end
    if getmetatable(rows_group) ~= Group or getmetatable(cols_group) ~= Group then
        error("Matrix:block_mayer_pair expects Group arguments")
    end
    check_matrix(other, "block_mayer_pair")
    rows_group:_check("block_mayer_pair")
    cols_group:_check("block_mayer_pair")
    return gaussian.matrix_block_mayer_pair(
        get_handle(self),
        get_handle(other),
        rows_group:_get_handle(),
        cols_group:_get_handle()
    )
end

function Matrix:block_norm_fro(rows_group, cols_group)
    check_matrix(self, "block_norm_fro")
    if getmetatable(rows_group) ~= Group or getmetatable(cols_group) ~= Group then
        error("Matrix:block_norm_fro expects Group arguments")
    end
    rows_group:_check("block_norm_fro")
    cols_group:_check("block_norm_fro")
    return gaussian.matrix_block_norm_fro(get_handle(self), rows_group:_get_handle(), cols_group:_get_handle())
end

function Matrix:free()
    if self._released then return end
    self.molecule:_release("matrix", self._id)
    self._released = true
    self._id = nil
end

function Matrix:is_freed()
    return self._released or self._id == nil
end

Matrix.__gc = function(self)
    if not self._released then
        pcall(self.free, self)
    end
end

Matrix.__tostring = function(self)
    if self._released then return "Matrix(released)" end
    local ok, rows, cols = pcall(self.size, self)
    if not ok then return "Matrix(invalid)" end
    return string.format("Matrix(%dx%d, type='%s', id=%d)", rows, cols, self._type or "?", self._id)
end

Matrix.__mul = function(A, B)
    if getmetatable(A) == Matrix and getmetatable(B) == Matrix then
        return A:mul(B)
    elseif getmetatable(A) == Matrix and type(B) == "number" then
        return A:scale(B)
    elseif type(A) == "number" and getmetatable(B) == Matrix then
        return B:scale(A)
    else
        error("Invalid multiplication")
    end
end

Matrix.__add = function(A, B)
    return A:add(B, 1.0, 1.0)
end

Matrix.__sub = function(A, B)
    return A:add(B, 1.0, -1.0)
end

Matrix.__unm = function(A)
    return A:scale(-1.0)
end

-- ============================================================================
-- Group: Lua-дескриптор битовой маски AO.
-- ============================================================================

Group = {}
Group.__index = Group

function Group._from_id(molecule, id)
    local entry = molecule:_get_entry("group", id)
    entry.refcount = entry.refcount + 1

    local self = setmetatable({
        molecule = molecule,
        _id = id,
        _released = false
    }, Group)
    return self
end

function Group:_check(method)
    if self._released then
        error("Group:" .. method .. " called on released group")
    end
    self.molecule:_touch("group", self._id)
end

function Group:_get_handle()
    return self.molecule:_touch("group", self._id)
end

function Group.create(nbasis, molecule)
    if molecule:is_closed() then error("Molecule is closed") end
    local h = gaussian.group_create(nbasis)
    if h == nil then error("group_create failed") end
    local id = molecule:_register("group", h)
    return Group._from_id(molecule, id)
end

function Group.create_full(nbasis, molecule)
    if molecule:is_closed() then error("Molecule is closed") end
    local h = gaussian.group_create_full(nbasis)
    if h == nil then error("group_create_full failed") end
    local id = molecule:_register("group", h)
    return Group._from_id(molecule, id)
end

function Group.from_atom(molecule, atom_idx_1based)
    if molecule:is_closed() then error("Molecule is closed") end
    local h = gaussian.group_from_atom(molecule.handle, atom_idx_1based - 1)
    if h == nil then error("group_from_atom failed") end
    local id = molecule:_register("group", h)
    return Group._from_id(molecule, id)
end

function Group.from_indices(indices, molecule)
    if molecule:is_closed() then error("Molecule is closed") end
    local nbasis = molecule:get_basis_size()
    local n = #indices
    local arr = ffi.new("int[?]", n)
    for i = 0, n - 1 do
        arr[i] = indices[i + 1] - 1
    end
    local h = gaussian.group_from_indices(arr, n, nbasis)
    if h == nil then error("group_from_indices failed") end
    local id = molecule:_register("group", h)
    return Group._from_id(molecule, id)
end

function Group:count()
    self:_check("count")
    return gaussian.group_count(self:_get_handle())
end

function Group:nbasis()
    self:_check("nbasis")
    return gaussian.group_nbasis(self:_get_handle())
end

function Group:set(index, value)
    self:_check("set")
    gaussian.group_set_bit(self:_get_handle(), index - 1, value and 1 or 0)
end

function Group:get(index)
    self:_check("get")
    local v = gaussian.group_get_bit(self:_get_handle(), index - 1)
    if v < 0 then error("group_get_bit failed") end
    return v == 1
end

function Group:union(other)
    self:_check("union")
    other:_check("union")
    local h = gaussian.group_or(self:_get_handle(), other:_get_handle())
    if h == nil then error("group_or failed") end
    local id = self.molecule:_register("group", h)
    return Group._from_id(self.molecule, id)
end

function Group:intersection(other)
    self:_check("intersection")
    other:_check("intersection")
    local h = gaussian.group_and(self:_get_handle(), other:_get_handle())
    if h == nil then error("group_and failed") end
    local id = self.molecule:_register("group", h)
    return Group._from_id(self.molecule, id)
end

function Group:free()
    if self._released then return end
    self.molecule:_release("group", self._id)
    self._released = true
    self._id = nil
end

function Group:is_freed()
    return self._released or self._id == nil
end

Group.__gc = function(self)
    if not self._released then
        pcall(self.free, self)
    end
end

Group.__tostring = function(self)
    if self._released then return "Group(released)" end
    local ok, c, n = pcall(function() return self:count(), self:nbasis() end)
    if not ok then return "Group(invalid)" end
    return string.format("Group(count=%d, nbasis=%d, id=%d)", c, n, self._id)
end

-- ============================================================================
-- Методы доступа к данным молекулы.
-- ============================================================================

function Molecule:get_basis_size()
    if self._closed then error("Molecule is closed") end
    return gaussian.gaussian_get_basis_size(self.handle)
end

function Molecule:get_num_atoms()
    if self._closed then error("Molecule is closed") end
    return gaussian.gaussian_get_num_atoms(self.handle)
end

function Molecule:get_atomic_number(atom_idx)
    if self._closed then error("Molecule is closed") end
    return gaussian.gaussian_get_atomic_number(self.handle, atom_idx - 1)
end

function Molecule:get_basis_name()
    if self._closed then error("Molecule is closed") end
    local cstr = gaussian.gaussian_get_basis_name(self.handle)
    if cstr == nil then return nil end
    local s = ffi.string(cstr)
    if s == "" then return nil end
    return s
end

function Molecule:get_nuclear_repulsion()
    if self._closed then error("Molecule is closed") end
    return gaussian.gaussian_get_nuclear_repulsion(self.handle)
end

function Molecule:ao_atom_mapping()
    if self._closed then error("Molecule is closed") end
    local nbasis = self:get_basis_size()
    local arr = ffi.new("int[?]", nbasis)
    local rc = gaussian.gaussian_get_ao_atom_mapping(self.handle, arr, nbasis)
    if rc ~= 0 then error("ao_atom_mapping failed") end
    local t = {}
    for i = 0, nbasis - 1 do
        t[i+1] = tonumber(arr[i]) + 1
    end
    return t
end

function Molecule:mo_coefficients()
    if self._closed then error("Molecule is closed") end
    local h = gaussian.gaussian_get_mo_coefficients(self.handle)
    if h == nil then return nil end
    local id = self:_register("matrix", h, {type = "mo_coeff"})
    return Matrix._from_id(self, id, "mo_coeff")
end

function Molecule:orbital_energies()
    if self._closed then error("Molecule is closed") end
    local nbasis = self:get_basis_size()
    local arr = ffi.new("double[?]", nbasis)
    local n = gaussian.gaussian_get_orbital_energies(self.handle, arr, nbasis)
    if n < 0 then return nil end
    local out = {}
    for i = 0, math.min(n, nbasis) - 1 do
        out[i + 1] = arr[i]
    end
    return out, n
end

-- Persistent-матрицы кэшируются по ID; каждый вызов возвращает новый Lua-дескриптор.
function Molecule:overlap()
    if self._closed then error("Molecule is closed") end
    local pid = self._registry.persistent_ids.overlap
    if not pid then
        local h = gaussian.gaussian_get_matrix(self.handle, "overlap")
        if h == nil then error("Cannot get overlap matrix") end
        pid = self:_register("matrix", h, {persistent = true, type = "overlap"})
        self._registry.persistent_ids.overlap = pid
    end
    -- Новый дескриптор увеличивает refcount в _from_id.
    return Matrix._from_id(self, pid, "overlap")
end

function Molecule:density()
    if self._closed then error("Molecule is closed") end
    local pid = self._registry.persistent_ids.density
    if not pid then
        local h = gaussian.gaussian_get_matrix(self.handle, "density")
        if h == nil then error("Cannot get density matrix") end
        pid = self:_register("matrix", h, {persistent = true, type = "density"})
        self._registry.persistent_ids.density = pid
    end
    return Matrix._from_id(self, pid, "density")
end

-- ============================================================================
-- Экспорт модуля.
-- ============================================================================

return {
    Molecule = Molecule,
    Matrix = Matrix,
    Group = Group,
    gaussian = gaussian
}
