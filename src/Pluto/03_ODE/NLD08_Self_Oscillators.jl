### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 04c345c5-fd34-5830-97c4-edadbb05055c
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots
end

# ╔═╡ 79124cae-eeb5-5104-9c14-f317202308dc
md"""
# Self-oscillation: Rayleigh reed

Negative damping feeds small motions and cubic damping limits their amplitude, producing self-oscillation.
"""

# ╔═╡ 0a72128d-f117-5d09-8fac-0050ff3d6dea
function model!(du, u, p, t)
    x, v = u
    μ, K, sc = p
    du[1] = v
    du[2] = -K*x + (μ-sc^2*v^2)*v
    nothing
end

# ╔═╡ 30daf2b7-9d99-596e-a0fe-822e9a7fcf59
u0 = [0.1, 0.0]

# ╔═╡ bf02afe1-302b-5060-abd2-0bf1eac1a3a3
tspan = (0.0, 100.0)

# ╔═╡ dfa42413-a9f5-5a28-b603-3e33e473504f
p = [0.3, 1.0, 1.0]  # μ, K, sc

# ╔═╡ 23966858-70a2-5835-bd94-373f629d0b02
prob = ODEProblem(model!, u0, tspan, p)

# ╔═╡ c4de1520-29c1-573a-a809-1d1f9c0b155b
sol = solve(prob, Tsit5(); abstol=1e-9, reltol=1e-7, saveat=0.02);

# ╔═╡ af23b3c0-fd6f-5eac-976b-ee99eb9ef52a
plot(sol; idxs=[1, 2], xlabel="t", ylabel="state")

# ╔═╡ 9bad4688-3b4c-55d1-8e8e-eedcc90151a5
plot(sol; idxs=(1, 2), legend=false)

# ╔═╡ Cell order:
# ╠═04c345c5-fd34-5830-97c4-edadbb05055c
# ╟─79124cae-eeb5-5104-9c14-f317202308dc
# ╠═0a72128d-f117-5d09-8fac-0050ff3d6dea
# ╠═30daf2b7-9d99-596e-a0fe-822e9a7fcf59
# ╠═bf02afe1-302b-5060-abd2-0bf1eac1a3a3
# ╠═dfa42413-a9f5-5a28-b603-3e33e473504f
# ╠═23966858-70a2-5835-bd94-373f629d0b02
# ╠═c4de1520-29c1-573a-a809-1d1f9c0b155b
# ╠═af23b3c0-fd6f-5eac-976b-ee99eb9ef52a
# ╠═9bad4688-3b4c-55d1-8e8e-eedcc90151a5
