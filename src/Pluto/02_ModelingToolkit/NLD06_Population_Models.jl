### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 68f13362-c8f3-57dd-9663-7eeadb9037e2
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots, ModelingToolkit
end

# ╔═╡ 056b5053-c22f-5e1c-b448-da2e617681c2
md"""
# Population growth with harvesting

A constant harvest lowers the logistic equilibrium; this example stops at extinction instead of continuing to negative population.
"""

# ╔═╡ fd5e701d-03ff-50e5-9143-b55d0be67d99
md"""
```math
\begin{aligned}
\dot{N} &= rN\left(1-\frac{N}{K}\right)-H.
\end{aligned}
```
"""

# ╔═╡ 52518659-834c-5fac-8c55-f61fadd80186
# Variables
begin
    @independent_variables t
    @variables N(t)
end

# ╔═╡ b407cfaf-d818-5082-afa7-a7848f3e75ed
# Parameters
@parameters r K H

# ╔═╡ fd822a3c-e1d9-567e-bb4e-46222640acab
# Time derivative
D = Differential(t)

# ╔═╡ cac348fa-3aa9-5a82-9afa-64a3f21ea857
equations = [
    D(N) ~ r*N*(1-N/K) - H
]

# ╔═╡ 221ef5e2-d45a-5564-86c6-f639ae8204be
@named system = ODESystem(equations, t)

# ╔═╡ bd024617-e133-58a0-a3d1-1d9470ac9272
simplified = structural_simplify(system)

# ╔═╡ 30cd1aa5-1402-5502-bb0d-438f8b5636a7
u0 = [N => 4.0]

# ╔═╡ 4873e92d-d82c-534c-9a11-0055e023448c
tspan = (0.0, 20.0)

# ╔═╡ 55a09355-3c16-5a8e-8593-c90c5877821c
p = [r => 1.0, K => 10.0, H => 1.5]

# ╔═╡ 34431db3-926c-5d2b-b730-b0289ef8deea
prob = ODEProblem(simplified, u0, tspan, p)

# ╔═╡ e7a57361-84f3-576b-8d0d-4751f817812f
population_index = findfirst(isequal(N), unknowns(simplified))

# ╔═╡ 412967f5-4733-5d15-b604-bc9070d18f5e
extinction = ContinuousCallback((u,t,integrator) -> u[population_index], nothing, terminate!)

# ╔═╡ bbe4ed49-5e4f-5c8d-9d56-9f3a58516d5b
sol = solve(prob, Tsit5(); abstol=1e-9, reltol=1e-7, callback=extinction, saveat=0.02);

# ╔═╡ 6433b736-275b-50df-b267-4d9aa0a39892
plot(sol; idxs=[N], xlabel="t", ylabel="state")

# ╔═╡ Cell order:
# ╠═68f13362-c8f3-57dd-9663-7eeadb9037e2
# ╟─056b5053-c22f-5e1c-b448-da2e617681c2
# ╟─fd5e701d-03ff-50e5-9143-b55d0be67d99
# ╠═52518659-834c-5fac-8c55-f61fadd80186
# ╠═b407cfaf-d818-5082-afa7-a7848f3e75ed
# ╠═fd822a3c-e1d9-567e-bb4e-46222640acab
# ╠═cac348fa-3aa9-5a82-9afa-64a3f21ea857
# ╠═221ef5e2-d45a-5564-86c6-f639ae8204be
# ╠═bd024617-e133-58a0-a3d1-1d9470ac9272
# ╠═30cd1aa5-1402-5502-bb0d-438f8b5636a7
# ╠═4873e92d-d82c-534c-9a11-0055e023448c
# ╠═55a09355-3c16-5a8e-8593-c90c5877821c
# ╠═34431db3-926c-5d2b-b730-b0289ef8deea
# ╠═e7a57361-84f3-576b-8d0d-4751f817812f
# ╠═412967f5-4733-5d15-b604-bc9070d18f5e
# ╠═bbe4ed49-5e4f-5c8d-9d56-9f3a58516d5b
# ╠═6433b736-275b-50df-b267-4d9aa0a39892
