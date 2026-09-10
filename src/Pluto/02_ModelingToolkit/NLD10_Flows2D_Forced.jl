### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ ac7f2521-b8d4-5325-8900-741d42e907a4
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots, ModelingToolkit
end

# ╔═╡ ea6459cf-14bf-5d87-b25a-2f840c9aa640
md"""
# Periodically forced Duffing oscillator

Periodic forcing sustains motion across the double well; sampling once per forcing period gives a Poincaré section.
"""

# ╔═╡ 84259192-162a-5b41-bef1-26037d81f20b
# Variables
begin
    @independent_variables t
    @variables x(t) v(t)
end

# ╔═╡ e6ffb7b8-c5a0-581b-9dd7-5827e20d0711
# Parameters
@parameters γ β F ω

# ╔═╡ 24deb66d-6967-51ce-9ed5-7da31e47451b
# Time derivative
D = Differential(t)

# ╔═╡ 9ac7ff3d-dda4-56b9-a76b-13b7d09954b8
equations = [
    D(x) ~ v,
    D(v) ~ -γ*v + β*x - x^3 + F*cos(ω*t)
]

# ╔═╡ 40049141-27c1-514b-b19f-1ca6a7fd6087
@named system = ODESystem(equations, t)

# ╔═╡ dd04c6e8-d9cc-5da1-a61c-09b6cb97f4f7
simplified = structural_simplify(system)

# ╔═╡ c470703c-a743-5755-a3f0-029b73169ada
u0 = [x => 0.1, v => 0.0]

# ╔═╡ 33dc34dc-dcd7-57ac-8f71-044b4811f08c
tspan = (0.0, 200.0)

# ╔═╡ 7efc732c-7351-58b0-8389-b55b3159b1c8
p = [γ => 0.2, β => 1.0, F => 0.3, ω => 1.2]

# ╔═╡ db4319a3-0a7f-577f-9b97-0d1d3167f4d8
prob = ODEProblem(simplified, u0, tspan, p)

# ╔═╡ c619ccf6-d5e9-59ea-8cfc-1c4db616481a
sol = solve(prob, Tsit5(); abstol=1e-9, reltol=1e-7, saveat=0.02);

# ╔═╡ 991b960d-e54f-554e-8a2c-1ec9f1fc3eb5
plot(sol; idxs=[x, v], xlabel="t", ylabel="state")

# ╔═╡ c49d47b6-09aa-5e15-a82d-5017222774bb
plot(sol; idxs=(x, v), legend=false)

# ╔═╡ 2963638a-58fd-5a55-b4a2-f026cc0affa8
let period = 2pi/(only(last(pair) for pair in p if isequal(first(pair), ω)))
    times = collect(20period:period:tspan[2])
    points = sol(times; idxs=[x, v])
    scatter(Array(points)[1,:], Array(points)[2,:];
        xlabel="x", ylabel="v", label="once per forcing period")
end

# ╔═╡ Cell order:
# ╠═ac7f2521-b8d4-5325-8900-741d42e907a4
# ╟─ea6459cf-14bf-5d87-b25a-2f840c9aa640
# ╠═84259192-162a-5b41-bef1-26037d81f20b
# ╠═e6ffb7b8-c5a0-581b-9dd7-5827e20d0711
# ╠═24deb66d-6967-51ce-9ed5-7da31e47451b
# ╠═9ac7ff3d-dda4-56b9-a76b-13b7d09954b8
# ╠═40049141-27c1-514b-b19f-1ca6a7fd6087
# ╠═dd04c6e8-d9cc-5da1-a61c-09b6cb97f4f7
# ╠═c470703c-a743-5755-a3f0-029b73169ada
# ╠═33dc34dc-dcd7-57ac-8f71-044b4811f08c
# ╠═7efc732c-7351-58b0-8389-b55b3159b1c8
# ╠═db4319a3-0a7f-577f-9b97-0d1d3167f4d8
# ╠═c619ccf6-d5e9-59ea-8cfc-1c4db616481a
# ╠═991b960d-e54f-554e-8a2c-1ec9f1fc3eb5
# ╠═c49d47b6-09aa-5e15-a82d-5017222774bb
# ╠═2963638a-58fd-5a55-b4a2-f026cc0affa8
