### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 5d02de77-e876-5180-92da-f00f47b29eeb
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots, ModelingToolkit
end

# ╔═╡ 4db81724-5e37-5de6-a381-f0d012c46aab
md"""
# Flows in 2D: damped oscillator

Damping turns the harmonic oscillator’s closed orbits into trajectories approaching the origin.
"""

# ╔═╡ 02dec78d-ae5b-5774-9142-92f95e9319b5
# Variables
begin
    @independent_variables t
    @variables x(t) v(t)
end

# ╔═╡ 7f2fdb00-8147-5de5-9c01-7874ddde01c4
# Parameters
@parameters k γ

# ╔═╡ 9b914812-8797-510d-914a-a04bbc1973c8
# Time derivative
D = Differential(t)

# ╔═╡ c4ba2c7e-6ccc-590f-84d2-763ca7f99c64
equations = [
    D(x) ~ v,
    D(v) ~ -k*x - γ*v
]

# ╔═╡ 6e724170-faa9-5aac-bb03-f4c5f460e2e7
@named system = ODESystem(equations, t)

# ╔═╡ 2c8d49e9-80cd-529a-9d3d-3477dc1b4264
simplified = structural_simplify(system)

# ╔═╡ c4713a29-629c-53f1-8ecd-8764fc2a59df
u0 = [x => 1.0, v => 0.0]

# ╔═╡ de90ba47-ecf0-5a53-9c75-a691726c2dda
tspan = (0.0, 30.0)

# ╔═╡ f37cdef0-da7e-539a-a2a3-fd1994167c77
p = [k => 1.0, γ => 0.15]

# ╔═╡ ee9e847a-fe92-5278-9536-59490771675e
prob = ODEProblem(simplified, u0, tspan, p)

# ╔═╡ bd5176e5-f824-5194-ac7a-95a9f892ef1b
sol = solve(prob, Tsit5(); abstol=1e-9, reltol=1e-7, saveat=0.02);

# ╔═╡ 28230d3e-4495-5541-b036-a55fa6ce1659
plot(sol; idxs=[x, v], xlabel="t", ylabel="state")

# ╔═╡ 5193f780-84fb-537a-bc6f-9638249499a2
plot(sol; idxs=(x, v), legend=false)

# ╔═╡ Cell order:
# ╠═5d02de77-e876-5180-92da-f00f47b29eeb
# ╟─4db81724-5e37-5de6-a381-f0d012c46aab
# ╠═02dec78d-ae5b-5774-9142-92f95e9319b5
# ╠═7f2fdb00-8147-5de5-9c01-7874ddde01c4
# ╠═9b914812-8797-510d-914a-a04bbc1973c8
# ╠═c4ba2c7e-6ccc-590f-84d2-763ca7f99c64
# ╠═6e724170-faa9-5aac-bb03-f4c5f460e2e7
# ╠═2c8d49e9-80cd-529a-9d3d-3477dc1b4264
# ╠═c4713a29-629c-53f1-8ecd-8764fc2a59df
# ╠═de90ba47-ecf0-5a53-9c75-a691726c2dda
# ╠═f37cdef0-da7e-539a-a2a3-fd1994167c77
# ╠═ee9e847a-fe92-5278-9536-59490771675e
# ╠═bd5176e5-f824-5194-ac7a-95a9f892ef1b
# ╠═28230d3e-4495-5541-b036-a55fa6ce1659
# ╠═5193f780-84fb-537a-bc6f-9638249499a2
