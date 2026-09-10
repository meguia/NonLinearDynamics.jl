### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 19c8f521-0f53-5605-8d98-168c4bc62a19
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots
end

# ╔═╡ 026e2d1b-6002-5fab-86c5-9f1ab58ddeb4
md"""
# Periodically forced Duffing oscillator

Periodic forcing sustains motion across the double well; sampling once per forcing period gives a Poincaré section.
"""

# ╔═╡ c46b85c7-8b5d-501d-b9cc-96ada7f08959
md"""
```math
\begin{aligned}
\dot{x} &= v, & \dot{v} &= -\gamma v+\beta x-x^3+F\cos(\omega t).
\end{aligned}
```
"""

# ╔═╡ 8094c580-0d07-55c8-8c80-613c8eb24bfc
function model!(du, u, p, t)
    x, v = u
    γ, β, F, ω = p
    du[1] = v
    du[2] = -γ*v + β*x - x^3 + F*cos(ω*t)
    nothing
end

# ╔═╡ 3e5e2194-e93e-53bb-b20c-c7d5f9a51871
u0 = [0.1, 0.0]

# ╔═╡ afb6eb92-3047-507d-aa31-dcab7ef550df
tspan = (0.0, 2000.0)

# ╔═╡ bed30766-4b12-5bbf-9e87-9159401a26a4
p = [0.15, 1.0, 0.3, 1.2]  # γ, β, F, ω

# ╔═╡ 666bd2cb-7448-5c85-bb6e-9512bb9dcded
prob = ODEProblem(model!, u0, tspan, p)

# ╔═╡ 59b12c63-d720-5986-a77f-46200f2a23b7
sol = solve(prob, Tsit5(); abstol=1e-9, reltol=1e-7, saveat=0.02);

# ╔═╡ 8a7f0137-2f1a-539c-a462-599444407b58
plot(sol; idxs=[1, 2], xlabel="t", ylabel="state")

# ╔═╡ 6958bc97-85ec-53a2-aaf5-6775e861f647
plot(sol; idxs=(1, 2), legend=false)

# ╔═╡ f45654b1-bb8e-5156-8c44-95b104195ef5
let period = 2pi/(p[4])
    times = collect(20period:period:tspan[2])
    points = sol(times; idxs=[1, 2])
    scatter(Array(points)[1,:], Array(points)[2,:];
        xlabel="x", ylabel="v", ms=2, msw=0, label="Poincare section")
end

# ╔═╡ Cell order:
# ╠═19c8f521-0f53-5605-8d98-168c4bc62a19
# ╟─026e2d1b-6002-5fab-86c5-9f1ab58ddeb4
# ╟─c46b85c7-8b5d-501d-b9cc-96ada7f08959
# ╠═8094c580-0d07-55c8-8c80-613c8eb24bfc
# ╠═3e5e2194-e93e-53bb-b20c-c7d5f9a51871
# ╠═afb6eb92-3047-507d-aa31-dcab7ef550df
# ╠═bed30766-4b12-5bbf-9e87-9159401a26a4
# ╠═666bd2cb-7448-5c85-bb6e-9512bb9dcded
# ╠═59b12c63-d720-5986-a77f-46200f2a23b7
# ╠═8a7f0137-2f1a-539c-a462-599444407b58
# ╠═6958bc97-85ec-53a2-aaf5-6775e861f647
# ╠═f45654b1-bb8e-5156-8c44-95b104195ef5
