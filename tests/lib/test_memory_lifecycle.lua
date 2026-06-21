local t = dofile((os.getenv("PROJECT_ROOT") or ".") .. "/tests/support.lua")

local chemistry = require("chemistry_ffi")

local function is_finite(value)
    return type(value) == "number" and value == value and value ~= math.huge and value ~= -math.huge
end

t.run_suite("memory_lifecycle", {
    {"temporary Matrix and Group stress", function()
        local iterations = tonumber(os.getenv("QA_MEMORY_STRESS_ITERS") or "500")

        for i = 1, iterations do
            local mol = assert(chemistry.Molecule.load(t.fixture))
            local P = mol:density()
            local S = mol:overlap()
            local PS = P * S
            local T = PS:transpose()
            local H = PS:hadamard(T)
            local root_pow = S:pow(0.5)
            local root_symm = S:symm_pow(0.5)
            local carbon = chemistry.Group.from_atom(mol, 1)
            local all = chemistry.Group.create_full(mol:get_basis_size(), mol)

            local trace = PS:trace()
            local norm = H:norm_fro()
            local block_norm = H:block_norm_fro(all, all)
            local carbon_norm = H:block_norm_fro(carbon, carbon)

            t.assert_true(is_finite(trace), "trace is finite")
            t.assert_true(is_finite(norm), "norm is finite")
            t.assert_true(is_finite(block_norm), "block norm is finite")
            t.assert_true(is_finite(carbon_norm), "atom block norm is finite")
            t.assert_near(trace, 10.0, 1e-3, "electron count remains stable")

            if i % 2 == 0 then
                all:free()
                carbon:free()
                root_symm:free()
                root_pow:free()
                H:free()
                T:free()
                PS:free()
                S:free()
                P:free()
            else
                all:free()
                H:free()
                PS:free()
                P = nil
                S = nil
                T = nil
                root_pow = nil
                root_symm = nil
                carbon = nil
            end

            mol:close()
            mol = nil

            if i % 50 == 0 then
                collectgarbage("collect")
            end
        end

        collectgarbage("collect")
    end},
})
