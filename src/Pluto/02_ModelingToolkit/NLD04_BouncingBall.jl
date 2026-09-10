### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 8f189a87-a105-54f1-a778-ca7e8e6b644c
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots, ModelingToolkit
end

# ╔═╡ e435737f-6596-5cb2-b2df-7015c42d1c77
md"""
# Bouncing ball: compliant contact

A spring acts only below the floor; its finite stiffness allows a small penetration during each bounce.
"""

# ╔═╡ 6754a089-f1a6-5c4c-836d-9db2fd008b4d
# Variables
begin
    @independent_variables t
    @variables h(t) v(t)
end

# ╔═╡ 19ea0a6e-55f6-5b7c-a05e-16b311897f55
# Parameters
@parameters g K γ

# ╔═╡ 555fb835-11e7-55ad-ad92-1063e4f27bce
# Time derivative
D = Differential(t)

# ╔═╡ 02dbdccf-6635-54b0-bd3b-3ff85596a65f
equations = [
    D(h) ~ v,
    D(v) ~ -g - K*min(h, 0) - γ*v
]

# ╔═╡ be68a442-4c83-532f-b270-c0f930cc2e25
@named system = ODESystem(equations, t)

# ╔═╡ a858c6d8-2d22-5122-b682-75c505658912
simplified = structural_simplify(system)

# ╔═╡ 736e53f4-bcd9-57ae-bc5d-9a5c7aee2aa9
u0 = [h => 1.0, v => 0.0]

# ╔═╡ ff5da50e-1af2-5669-92e1-82258c82136d
tspan = (0.0, 5.0)

# ╔═╡ 0e15c98b-8f81-590f-8909-efd5cb9396ed
p = [g => 9.8, K => 1000.0, γ => 0.2]

# ╔═╡ 5722d552-4cb3-528c-883b-32b985e71ce8
prob = ODEProblem(simplified, u0, tspan, p)

# ╔═╡ a85f6095-fef3-5471-83e5-e61dafff4c29
sol = solve(prob, Rosenbrock23(); abstol=1e-9, reltol=1e-7, saveat=0.02, dtmax=0.002);

# ╔═╡ a4987814-a342-508c-a746-3fbc23c817a3
plot(sol; idxs=[h, v], xlabel="t", ylabel="state")

# ╔═╡ 91ac87be-7db6-5e2d-b157-2b9f531ac18a
plot(sol; idxs=(h, v), legend=false)

# ╔═╡ Cell order:
# ╠═8f189a87-a105-54f1-a778-ca7e8e6b644c
# ╟─e435737f-6596-5cb2-b2df-7015c42d1c77
# ╠═6754a089-f1a6-5c4c-836d-9db2fd008b4d
# ╠═19ea0a6e-55f6-5b7c-a05e-16b311897f55
# ╠═555fb835-11e7-55ad-ad92-1063e4f27bce
# ╠═02dbdccf-6635-54b0-bd3b-3ff85596a65f
# ╠═be68a442-4c83-532f-b270-c0f930cc2e25
# ╠═a858c6d8-2d22-5122-b682-75c505658912
# ╠═736e53f4-bcd9-57ae-bc5d-9a5c7aee2aa9
# ╠═ff5da50e-1af2-5669-92e1-82258c82136d
# ╠═0e15c98b-8f81-590f-8909-efd5cb9396ed
# ╠═5722d552-4cb3-528c-883b-32b985e71ce8
# ╠═a85f6095-fef3-5471-83e5-e61dafff4c29
# ╠═a4987814-a342-508c-a746-3fbc23c817a3
# ╠═91ac87be-7db6-5e2d-b157-2b9f531ac18a
