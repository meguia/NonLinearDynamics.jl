### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ b3663b29-e52c-5f40-b346-34f42587f63c
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots, ModelingToolkit
end

# ╔═╡ 1ed575ee-62ff-568c-921a-daa6eabb5528
md"""
# One-dimensional bifurcations: pitchfork

Change μ through zero: the stable origin gives way to two stable equilibria in this supercritical pitchfork.
"""

# ╔═╡ 254cb1bc-2bda-5ac9-8c87-74663ef288d4
# Variables
begin
    @independent_variables t
    @variables x(t)
end

# ╔═╡ 1d8c8838-a425-515b-807b-bddb551f7b9a
# Parameters
@parameters μ

# ╔═╡ c7dca9bd-aa8a-5739-b813-e7f976fd562e
# Time derivative
D = Differential(t)

# ╔═╡ 1f025e78-c985-5533-9186-3f27eda526ab
equations = [
    D(x) ~ μ*x - x^3
]

# ╔═╡ 6b1d4177-0366-53bd-83ba-0dec43446bca
@named system = ODESystem(equations, t)

# ╔═╡ 358a5300-550a-5a4a-a2e2-0026f7f65ab0
simplified = structural_simplify(system)

# ╔═╡ e22fd2fd-3cdd-51a6-b50d-a7a73b5ce7e7
u0 = [x => 0.1]

# ╔═╡ 1c581090-de0d-5f88-8f5e-84ba344f7f1b
tspan = (0.0, 30.0)

# ╔═╡ 36d8a749-f199-5485-acd1-27aab9ca3be9
p = [μ => 0.5]

# ╔═╡ 8439ad75-035e-5581-9cf5-ab29e4ee5e7f
prob = ODEProblem(simplified, u0, tspan, p)

# ╔═╡ 211e59ed-cb23-58d5-b422-bab23ab8a2a6
sol = solve(prob, Tsit5(); abstol=1e-9, reltol=1e-7, saveat=0.02);

# ╔═╡ 187a46ea-bdaf-556a-b73c-09d77cb3da89
plot(sol; idxs=[x], xlabel="t", ylabel="state")

# ╔═╡ 7313ab16-9814-5296-8f37-b1abfc44763f
let μs = range(-1, 1; length=401)
    figure = plot(μs, zero.(μs); label="origin", xlabel="μ", ylabel="equilibrium x")
    positive = filter(>=(0), μs)
    plot!(figure, positive, sqrt.(positive); label="positive branch")
    plot!(figure, positive, -sqrt.(positive); label="negative branch")
end

# ╔═╡ Cell order:
# ╠═b3663b29-e52c-5f40-b346-34f42587f63c
# ╟─1ed575ee-62ff-568c-921a-daa6eabb5528
# ╠═254cb1bc-2bda-5ac9-8c87-74663ef288d4
# ╠═1d8c8838-a425-515b-807b-bddb551f7b9a
# ╠═c7dca9bd-aa8a-5739-b813-e7f976fd562e
# ╠═1f025e78-c985-5533-9186-3f27eda526ab
# ╠═6b1d4177-0366-53bd-83ba-0dec43446bca
# ╠═358a5300-550a-5a4a-a2e2-0026f7f65ab0
# ╠═e22fd2fd-3cdd-51a6-b50d-a7a73b5ce7e7
# ╠═1c581090-de0d-5f88-8f5e-84ba344f7f1b
# ╠═36d8a749-f199-5485-acd1-27aab9ca3be9
# ╠═8439ad75-035e-5581-9cf5-ab29e4ee5e7f
# ╠═211e59ed-cb23-58d5-b422-bab23ab8a2a6
# ╠═187a46ea-bdaf-556a-b73c-09d77cb3da89
# ╠═7313ab16-9814-5296-8f37-b1abfc44763f
