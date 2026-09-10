### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ ff555515-e183-5212-836c-842b5fb4e043
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots
end

# ╔═╡ 4b6c1283-065c-5492-8815-05773296f288
md"""
# Bifurcations: the Bogdanov–Takens unfolding

This is the model used in interactive NLD12: change μ1 to explore its trajectories, then use that lesson for numerical continuation.
"""

# ╔═╡ da860a78-a2d3-5d01-b0e2-3611026f6c7f
function model!(du, u, p, t)
    x, v = u
    μ1, μ2 = p
    du[1] = v
    du[2] = μ1 + x*(μ2-v+x*(1-x-v))
    nothing
end

# ╔═╡ 1fbfcdab-37de-5762-97b8-b287998ca1eb
u0 = [0.1, 0.0]

# ╔═╡ 87f60eb2-2360-5c27-8256-3831e3f095eb
tspan = (0.0, 1000.0)

# ╔═╡ 4a92bfa0-2d35-50b3-999f-d0a08d838804
p = [-0.02, -0.16]  # μ1, μ2

# ╔═╡ 930b0512-c8d5-57ce-9091-8c9ade3935dc
prob = ODEProblem(model!, u0, tspan, p)

# ╔═╡ 54c1cff6-bd44-50a5-a6ce-bca3b4e62e0e
sol = solve(prob, Tsit5(); abstol=1e-9, reltol=1e-7, saveat=0.02);

# ╔═╡ 023ada46-bf58-5ad0-b33a-fc6a8401a79d
plot(sol; idxs=[1, 2], xlabel="t", ylabel="state")

# ╔═╡ 2ad41245-d960-5c10-9ff3-cf349f6dbc88
plot(sol; idxs=(1, 2), legend=false)

# ╔═╡ Cell order:
# ╠═ff555515-e183-5212-836c-842b5fb4e043
# ╟─4b6c1283-065c-5492-8815-05773296f288
# ╠═da860a78-a2d3-5d01-b0e2-3611026f6c7f
# ╠═1fbfcdab-37de-5762-97b8-b287998ca1eb
# ╠═87f60eb2-2360-5c27-8256-3831e3f095eb
# ╠═4a92bfa0-2d35-50b3-999f-d0a08d838804
# ╠═930b0512-c8d5-57ce-9091-8c9ade3935dc
# ╠═54c1cff6-bd44-50a5-a6ce-bca3b4e62e0e
# ╠═023ada46-bf58-5ad0-b33a-fc6a8401a79d
# ╠═2ad41245-d960-5c10-9ff3-cf349f6dbc88
