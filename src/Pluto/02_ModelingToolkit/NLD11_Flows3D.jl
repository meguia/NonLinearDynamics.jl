### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 05c4fe8b-54ae-56ba-9f21-22ddd5de4508
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots, ModelingToolkit
end

# ╔═╡ eeb1e51a-4bc4-57b6-9eee-4cd4a56d439b
md"""
# Flows in 3D: Lorenz system

Three coupled equations produce the Lorenz attractor; short-time agreement is meaningful even when long chaotic trajectories separate.
"""

# ╔═╡ ae71c4bf-bbfb-5405-9e4b-3f1308237d5a
# Variables
begin
    @independent_variables t
    @variables x(t) y(t) z(t)
end

# ╔═╡ 14866df1-b0ab-5e95-ae65-9f78e349a049
# Parameters
@parameters σ ρ β

# ╔═╡ 6d626cef-b3eb-5372-9743-c46bdc8d7375
# Time derivative
D = Differential(t)

# ╔═╡ f0f0f12d-810e-537a-9366-b7c3174c7bbd
equations = [
    D(x) ~ σ*(y-x),
    D(y) ~ x*(ρ-z)-y,
    D(z) ~ x*y-β*z
]

# ╔═╡ bd5d2ef3-d34c-5311-8d70-182c9e9ee634
@named system = ODESystem(equations, t)

# ╔═╡ 370db2c7-062a-51b1-853d-f48df2234260
simplified = structural_simplify(system)

# ╔═╡ 98d51537-9f1c-522c-baaf-1cfa4ae21336
u0 = [x => 1.0, y => 0.0, z => 0.0]

# ╔═╡ 24c02a8c-5434-5277-bb9f-18090289521a
tspan = (0.0, 40.0)

# ╔═╡ 938b738d-611b-5422-bc7d-dfae96d4269f
p = [σ => 10.0, ρ => 28.0, β => 8/3]

# ╔═╡ bf52676b-1e97-5f43-930a-c96e2c7567e3
prob = ODEProblem(simplified, u0, tspan, p)

# ╔═╡ 5bfede41-9413-5413-9ae9-ae2a95043203
sol = solve(prob, Tsit5(); abstol=1e-9, reltol=1e-7, saveat=0.02);

# ╔═╡ 9957df2f-64df-5a96-a228-a607727cd274
plot(sol; idxs=[x, y, z], xlabel="t", ylabel="state")

# ╔═╡ eccd5980-29b0-537b-a491-aa729fce8357
plot(sol; idxs=(x, y, z), legend=false)

# ╔═╡ Cell order:
# ╠═05c4fe8b-54ae-56ba-9f21-22ddd5de4508
# ╟─eeb1e51a-4bc4-57b6-9eee-4cd4a56d439b
# ╠═ae71c4bf-bbfb-5405-9e4b-3f1308237d5a
# ╠═14866df1-b0ab-5e95-ae65-9f78e349a049
# ╠═6d626cef-b3eb-5372-9743-c46bdc8d7375
# ╠═f0f0f12d-810e-537a-9366-b7c3174c7bbd
# ╠═bd5d2ef3-d34c-5311-8d70-182c9e9ee634
# ╠═370db2c7-062a-51b1-853d-f48df2234260
# ╠═98d51537-9f1c-522c-baaf-1cfa4ae21336
# ╠═24c02a8c-5434-5277-bb9f-18090289521a
# ╠═938b738d-611b-5422-bc7d-dfae96d4269f
# ╠═bf52676b-1e97-5f43-930a-c96e2c7567e3
# ╠═5bfede41-9413-5413-9ae9-ae2a95043203
# ╠═9957df2f-64df-5a96-a228-a607727cd274
# ╠═eccd5980-29b0-537b-a491-aa729fce8357
