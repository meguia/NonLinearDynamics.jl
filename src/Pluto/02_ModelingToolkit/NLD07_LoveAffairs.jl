### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ c7efb5d3-9d22-5a51-bea0-21430ae5afd6
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots, ModelingToolkit
end

# ╔═╡ 35ab07b6-7d1b-5d81-8e23-eb2d0bf313db
md"""
# Love affairs: a linear system

The coefficients a and d describe each person’s own response, while b and c describe their response to the partner.
"""

# ╔═╡ d2311d25-06f9-5bf1-ad81-f2d9588c40a4
# Variables
begin
    @independent_variables t
    @variables R(t) J(t)
end

# ╔═╡ 8e118e71-cef6-51be-9955-b627fc83cd14
# Parameters
@parameters a b c d

# ╔═╡ 6ded1864-fa6d-5dc5-ae41-a805f32b8d2b
# Time derivative
D = Differential(t)

# ╔═╡ cd162748-d5db-5e4a-905d-a04b97d6383d
equations = [
    D(R) ~ a*R + b*J,
    D(J) ~ c*R + d*J
]

# ╔═╡ 4a6a4bc4-e7df-52d5-9a57-afde5ee44817
@named system = ODESystem(equations, t)

# ╔═╡ 588b5359-d5d4-5ee3-a7a7-9fe0077bf29a
simplified = structural_simplify(system)

# ╔═╡ 542563d9-a65d-518c-a2b0-8357ed5c2452
u0 = [R => 1.0, J => 0.0]

# ╔═╡ bbc192af-cb07-56e1-b2a1-738f5b47f1dc
tspan = (0.0, 30.0)

# ╔═╡ 95689f62-e9f0-54fe-9110-693f4d622c8e
p = [a => -0.1, b => 1.0, c => -1.0, d => -0.1]

# ╔═╡ 0d2d1757-07fd-51b7-b28d-7b181f242e74
prob = ODEProblem(simplified, u0, tspan, p)

# ╔═╡ 296a65e4-c11e-5d95-bd3d-0feb6b82d02b
sol = solve(prob, Tsit5(); abstol=1e-9, reltol=1e-7, saveat=0.02);

# ╔═╡ a4299bc0-865c-513f-bc2c-b4c9a9cd93a5
plot(sol; idxs=[R, J], xlabel="t", ylabel="state")

# ╔═╡ 29ae7e4d-f7a3-5ccb-8d2c-12119950fb86
plot(sol; idxs=(R, J), legend=false)

# ╔═╡ Cell order:
# ╠═c7efb5d3-9d22-5a51-bea0-21430ae5afd6
# ╟─35ab07b6-7d1b-5d81-8e23-eb2d0bf313db
# ╠═d2311d25-06f9-5bf1-ad81-f2d9588c40a4
# ╠═8e118e71-cef6-51be-9955-b627fc83cd14
# ╠═6ded1864-fa6d-5dc5-ae41-a805f32b8d2b
# ╠═cd162748-d5db-5e4a-905d-a04b97d6383d
# ╠═4a6a4bc4-e7df-52d5-9a57-afde5ee44817
# ╠═588b5359-d5d4-5ee3-a7a7-9fe0077bf29a
# ╠═542563d9-a65d-518c-a2b0-8357ed5c2452
# ╠═bbc192af-cb07-56e1-b2a1-738f5b47f1dc
# ╠═95689f62-e9f0-54fe-9110-693f4d622c8e
# ╠═0d2d1757-07fd-51b7-b28d-7b181f242e74
# ╠═296a65e4-c11e-5d95-bd3d-0feb6b82d02b
# ╠═a4299bc0-865c-513f-bc2c-b4c9a9cd93a5
# ╠═29ae7e4d-f7a3-5ccb-8d2c-12119950fb86
