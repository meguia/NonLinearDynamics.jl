### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 1913f953-d75c-5dcf-85c3-fd4e3f96a91d
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots, ModelingToolkit
end

# ╔═╡ 500a6ff3-d550-56b1-a366-86b910b60fe5
md"""
# Self-oscillation: Rayleigh reed

Negative damping feeds small motions and cubic damping limits their amplitude, producing self-oscillation.
"""

# ╔═╡ c6fad613-67de-570c-80a2-e60896a15a15
# Variables
begin
    @independent_variables t
    @variables x(t) v(t)
end

# ╔═╡ 2717f8f0-ab12-5adb-9758-d558db295cea
# Parameters
@parameters μ K sc

# ╔═╡ bcc09890-28b2-5c8c-bae3-1f705929790b
# Time derivative
D = Differential(t)

# ╔═╡ f308b885-e463-5956-a8e6-abc84edd979a
equations = [
    D(x) ~ v,
    D(v) ~ -K*x + (μ-sc^2*v^2)*v
]

# ╔═╡ 3b899fe6-1614-5459-9cf5-bd893264ae87
@named system = ODESystem(equations, t)

# ╔═╡ 3281d720-6313-5ae9-9d16-6999f1bd5762
simplified = structural_simplify(system)

# ╔═╡ 5948c899-0d4b-518d-8772-de200f7cfb9c
u0 = [x => 0.1, v => 0.0]

# ╔═╡ 0724d600-26b5-51d9-b122-a79317e1e83c
tspan = (0.0, 100.0)

# ╔═╡ 2b612fca-6af4-59af-9c53-64891ef48ce1
p = [μ => 0.3, K => 1.0, sc => 1.0]

# ╔═╡ 1bd9ccd5-eaeb-54ac-a7b1-26eb606e95ca
prob = ODEProblem(simplified, u0, tspan, p)

# ╔═╡ 960624e0-18c1-54f9-9875-a3462698cb9b
sol = solve(prob, Tsit5(); abstol=1e-9, reltol=1e-7, saveat=0.02);

# ╔═╡ cc503c76-522c-59f9-a594-85433db4c159
plot(sol; idxs=[x, v], xlabel="t", ylabel="state")

# ╔═╡ e035a626-5463-55b4-857c-24b07aad4519
plot(sol; idxs=(x, v), legend=false)

# ╔═╡ Cell order:
# ╠═1913f953-d75c-5dcf-85c3-fd4e3f96a91d
# ╟─500a6ff3-d550-56b1-a366-86b910b60fe5
# ╠═c6fad613-67de-570c-80a2-e60896a15a15
# ╠═2717f8f0-ab12-5adb-9758-d558db295cea
# ╠═bcc09890-28b2-5c8c-bae3-1f705929790b
# ╠═f308b885-e463-5956-a8e6-abc84edd979a
# ╠═3b899fe6-1614-5459-9cf5-bd893264ae87
# ╠═3281d720-6313-5ae9-9d16-6999f1bd5762
# ╠═5948c899-0d4b-518d-8772-de200f7cfb9c
# ╠═0724d600-26b5-51d9-b122-a79317e1e83c
# ╠═2b612fca-6af4-59af-9c53-64891ef48ce1
# ╠═1bd9ccd5-eaeb-54ac-a7b1-26eb606e95ca
# ╠═960624e0-18c1-54f9-9875-a3462698cb9b
# ╠═cc503c76-522c-59f9-a594-85433db4c159
# ╠═e035a626-5463-55b4-857c-24b07aad4519
