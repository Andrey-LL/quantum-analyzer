local t = dofile((os.getenv("PROJECT_ROOT") or ".") .. "/tests/support.lua")

local viz = require("viz_gnuplot")

local function assert_svg(path)
    local f = io.open(path, "rb")
    t.assert_true(f ~= nil, "svg file exists")
    if f then
        local head = f:read(120) or ""
        f:close()
        t.assert_true(head:find("<svg", 1, true) ~= nil, "svg header")
    end
end

t.run_suite("viz_gnuplot", {
    {"bar returns custom link and writes output", function()
        local out = t.project_root .. "/.cache/test_outputs/viz/bar.svg"
        os.remove(out)

        local md = viz.bar({
            {label = "C", value = -0.5781},
            {label = "H1", value = 0.1445},
            {label = "H2", value = 0.1445},
        }, {
            caption = "Mulliken charges",
            output = out,
            link = "bar.svg",
            xlabel = "Atom",
            ylabel = "Charge",
        })

        t.assert_eq(md, "![Mulliken charges](bar.svg)", "markdown link")
        assert_svg(out)
    end},

    {"heatmap returns custom link and writes output", function()
        local out = t.project_root .. "/.cache/test_outputs/viz/heatmap.svg"
        os.remove(out)

        local matrix = {
            size = function() return 2, 2 end,
            get = function(_, i, j)
                if i == j then return 1.0 end
                return 0.25
            end,
        }

        local md = viz.heatmap(matrix, {
            caption = "Small heatmap",
            output = out,
            link = "heatmap.svg",
        })

        t.assert_eq(md, "![Small heatmap](heatmap.svg)", "markdown link")
        assert_svg(out)
    end},
})
