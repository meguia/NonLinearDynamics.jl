### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 3097502a-2441-5747-bce3-8dbf282381cf
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots, ModelingToolkit
end

# ╔═╡ ff2e4910-ec17-5060-b1b9-e42a3e65bc64
md"""
# Bifurcations: the Bogdanov–Takens model

This is the model used in interactive NLD12: change μ1 to explore its trajectories, then use that lesson for numerical continuation.
"""

# ╔═╡ 66df619b-28cb-53fc-8875-364291ef0a19
# Variables
begin
    @independent_variables t
    @variables x(t) v(t)
end

# ╔═╡ ee79c88a-ef43-5ffb-9790-fe2166f8fea4
# Parameters
@parameters μ1 μ2

# ╔═╡ 1b05853c-7f5b-575a-bdb2-bcedba9afe45
# Time derivative
D = Differential(t)

# ╔═╡ 72a8e160-78b4-5eed-8c29-b8aa806642ed
equations = [
    D(x) ~ v,
    D(v) ~ μ1 + x*(μ2-v+x*(1-x-v))
]

# ╔═╡ 5902957b-9c5c-5f75-9f51-fbe64f163e17
@named system = ODESystem(equations, t)

# ╔═╡ 2f587b3b-0815-5a7d-b251-b03662b9f7b5
simplified = structural_simplify(system)

# ╔═╡ b0da7319-8d46-5092-9afb-42012a34cbbe
u0 = [x => 0.1, v => 0.0]

# ╔═╡ 15a968bd-d4e2-58d6-b68e-04786ab84569
tspan = (0.0, 100.0)

# ╔═╡ 1cbec4bb-bc36-58cc-b25e-ee6a2dad1594
p = [μ1 => -0.02, μ2 => -0.16]

# ╔═╡ 6a7491d9-d8a8-5bd2-bf4f-47d9e244aa98
prob = ODEProblem(simplified, u0, tspan, p)

# ╔═╡ 414ee0fa-612f-5ab6-8a01-9ae45ce64c61
sol = solve(prob, Tsit5(); abstol=1e-9, reltol=1e-7, saveat=0.02);

# ╔═╡ 01ab060d-04d4-545a-bbe6-b924fbfd0010
plot(sol; idxs=[x, v], xlabel="t", ylabel="state")

# ╔═╡ bb65e6e1-29be-58d3-a68d-d791284564eb
plot(sol; idxs=(x, v), legend=false)

# ╔═╡ Cell order:
# ╠═3097502a-2441-5747-bce3-8dbf282381cf
# ╟─ff2e4910-ec17-5060-b1b9-e42a3e65bc64
# ╠═66df619b-28cb-53fc-8875-364291ef0a19
# ╠═ee79c88a-ef43-5ffb-9790-fe2166f8fea4
# ╠═1b05853c-7f5b-575a-bdb2-bcedba9afe45
# ╠═72a8e160-78b4-5eed-8c29-b8aa806642ed
# ╠═5902957b-9c5c-5f75-9f51-fbe64f163e17
# ╠═2f587b3b-0815-5a7d-b251-b03662b9f7b5
# ╠═b0da7319-8d46-5092-9afb-42012a34cbbe
# ╠═15a968bd-d4e2-58d6-b68e-04786ab84569
# ╠═1cbec4bb-bc36-58cc-b25e-ee6a2dad1594
# ╠═6a7491d9-d8a8-5bd2-bf4f-47d9e244aa98
# ╠═414ee0fa-612f-5ab6-8a01-9ae45ce64c61
# ╠═01ab060d-04d4-545a-bbe6-b924fbfd0010
# ╠═bb65e6e1-29be-58d3-a68d-d791284564eb
