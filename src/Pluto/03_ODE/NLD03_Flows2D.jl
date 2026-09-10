### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ c54bb4b0-bec2-51bb-85b2-8a3590b080a4
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots
end

# ╔═╡ e5bb4bb8-33da-51a4-9b15-4682aaddb356
md"""
# Flows in 2D: damped oscillator

Damping turns the harmonic oscillator’s closed orbits into trajectories approaching the origin.
"""

# ╔═╡ 36f31570-0571-5bc7-bbd6-72723edab481
md"""
```math
\begin{aligned}
\dot{x} &= v, & \dot{v} &= -kx-\gamma v.
\end{aligned}
```
"""

# ╔═╡ f2f50938-2660-5469-9566-d5f8cd81d6bc
function model!(du, u, p, t)
    x, v = u
    k, γ = p
    du[1] = v
    du[2] = -k*x - γ*v
    nothing
end

# ╔═╡ 89f5b149-ac4b-5568-9fcd-ecb7afb95d75
u0 = [1.0, 0.0]

# ╔═╡ 625b7408-2575-5935-aac1-ab893a7defc0
tspan = (0.0, 30.0)

# ╔═╡ 288a546b-d6d0-53fa-aefb-60dbbc4046c1
p = [1.0, 0.15]  # k, γ

# ╔═╡ 6c495b7d-ca5b-50e7-a22c-a1bf0be9a8ab
prob = ODEProblem(model!, u0, tspan, p)

# ╔═╡ 1c17895b-de8e-5de0-9dc2-4c7d37075a67
sol = solve(prob, Tsit5(); abstol=1e-9, reltol=1e-7, saveat=0.02);

# ╔═╡ 045859c2-40a0-5c51-a906-3d066a1be3a4
plot(sol; idxs=[1, 2], xlabel="t", ylabel="state")

# ╔═╡ b1fcd466-82e5-5f70-814d-5f1bd10c22ec
plot(sol; idxs=(1, 2), legend=false)

# ╔═╡ Cell order:
# ╠═c54bb4b0-bec2-51bb-85b2-8a3590b080a4
# ╟─e5bb4bb8-33da-51a4-9b15-4682aaddb356
# ╟─36f31570-0571-5bc7-bbd6-72723edab481
# ╠═f2f50938-2660-5469-9566-d5f8cd81d6bc
# ╠═89f5b149-ac4b-5568-9fcd-ecb7afb95d75
# ╠═625b7408-2575-5935-aac1-ab893a7defc0
# ╠═288a546b-d6d0-53fa-aefb-60dbbc4046c1
# ╠═6c495b7d-ca5b-50e7-a22c-a1bf0be9a8ab
# ╠═1c17895b-de8e-5de0-9dc2-4c7d37075a67
# ╠═045859c2-40a0-5c51-a906-3d066a1be3a4
# ╠═b1fcd466-82e5-5f70-814d-5f1bd10c22ec
