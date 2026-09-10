### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 01e969ca-39f5-5a3b-894b-8e01b17f5923
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots, ModelingToolkit
end

# ╔═╡ 3fe09a16-5545-5a66-aed2-8bed4f998959
md"""
# The Duffing oscillator

A double-well potential gives two stable resting positions; change the initial state to explore their basins.
"""

# ╔═╡ ed54c963-a36a-57bd-9c5c-c6c2ca6ae354
md"""
```math
\begin{aligned}
\dot{x} &= v, & \dot{v} &= -\gamma v+\beta x-x^3.
\end{aligned}
```
"""

# ╔═╡ fc191cf0-2cf4-566a-938a-11e78ec4fd2c
# Variables
begin
    @independent_variables t
    @variables x(t) v(t)
end

# ╔═╡ 2018fff6-3db6-57cf-acc4-85e7a89a1eb5
# Parameters
@parameters γ β

# ╔═╡ 8b415f2f-b81c-54b7-9391-8fa80ddaca4e
# Time derivative
D = Differential(t)

# ╔═╡ 0b1a69f0-a1b3-51c9-8267-5ae75a3aa623
equations = [
    D(x) ~ v,
    D(v) ~ -γ*v + β*x - x^3
]

# ╔═╡ 9af91a9f-ef55-5c33-94c9-378dcbe5ad98
@named system = ODESystem(equations, t)

# ╔═╡ 59ff61ca-d831-5073-8b40-562c6ced1b43
simplified = structural_simplify(system)

# ╔═╡ cf6dfe5e-68a7-552c-b8fd-1ebf166dc5c4
u0 = [x => 0.1, v => 0.7]

# ╔═╡ 4957ef8f-d839-5b05-b2a4-173b717ef999
tspan = (0.0, 60.0)

# ╔═╡ 2ba12462-b6dd-54ee-91a6-ee9de422d517
p = [γ => 0.15, β => 1.0]

# ╔═╡ b0862ac0-4a73-5e0f-90f6-2f0d1c5c4f3f
prob = ODEProblem(simplified, u0, tspan, p)

# ╔═╡ 4811718e-aeb7-5b98-9f14-d627a8e4bd61
sol = solve(prob, Tsit5(); abstol=1e-9, reltol=1e-7, saveat=0.02);

# ╔═╡ 1477f42b-7aa3-574b-9211-3f19b0528aef
plot(sol; idxs=[x, v], xlabel="t", ylabel="state")

# ╔═╡ 2cfcd68b-7d7f-5c70-9727-d33ed406904d
plot(sol; idxs=(x, v), legend=false)

# ╔═╡ Cell order:
# ╠═01e969ca-39f5-5a3b-894b-8e01b17f5923
# ╟─3fe09a16-5545-5a66-aed2-8bed4f998959
# ╟─ed54c963-a36a-57bd-9c5c-c6c2ca6ae354
# ╠═fc191cf0-2cf4-566a-938a-11e78ec4fd2c
# ╠═2018fff6-3db6-57cf-acc4-85e7a89a1eb5
# ╠═8b415f2f-b81c-54b7-9391-8fa80ddaca4e
# ╠═0b1a69f0-a1b3-51c9-8267-5ae75a3aa623
# ╠═9af91a9f-ef55-5c33-94c9-378dcbe5ad98
# ╠═59ff61ca-d831-5073-8b40-562c6ced1b43
# ╠═cf6dfe5e-68a7-552c-b8fd-1ebf166dc5c4
# ╠═4957ef8f-d839-5b05-b2a4-173b717ef999
# ╠═2ba12462-b6dd-54ee-91a6-ee9de422d517
# ╠═b0862ac0-4a73-5e0f-90f6-2f0d1c5c4f3f
# ╠═4811718e-aeb7-5b98-9f14-d627a8e4bd61
# ╠═1477f42b-7aa3-574b-9211-3f19b0528aef
# ╠═2cfcd68b-7d7f-5c70-9727-d33ed406904d
