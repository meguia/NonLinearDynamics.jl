### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ db4f14e7-32f9-5c87-93ee-f988cc5aacca
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots
end

# ╔═╡ 62518ca0-5495-5907-93cd-6535924c39f1
md"""
# One-dimensional bifurcations: pitchfork

Change μ through zero: the stable origin gives way to two stable equilibria in this supercritical pitchfork.
"""

# ╔═╡ 800d9c5c-bae0-5ed5-91bb-8fdc7fcade41
md"""
```math
\begin{aligned}
\dot{x} &= \mu x-x^3.
\end{aligned}
```
"""

# ╔═╡ 9bd171c2-0dbd-5768-91f4-a5023265fbb0
function model!(du, u, p, t)
    x = only(u)
    μ = only(p)
    du[1] = μ*x - x^3
    nothing
end

# ╔═╡ 7bb14eaf-4b73-5036-a893-afbbcbcc932e
u0 = [0.1]

# ╔═╡ d13525f4-c075-5f9b-a75d-e1a55ca1e303
tspan = (0.0, 30.0)

# ╔═╡ f00db1e4-6298-5fb1-b3fe-df8bc73d8365
p = [0.5]  # μ

# ╔═╡ ba0270a1-ef90-586c-8682-fd104d5b00d1
prob = ODEProblem(model!, u0, tspan, p)

# ╔═╡ 23ede884-9339-547d-8047-550719420738
sol = solve(prob, Tsit5(); abstol=1e-9, reltol=1e-7, saveat=0.02);

# ╔═╡ 50fca1fa-bcab-5504-af22-0302e55ea7df
plot(sol; idxs=[1], xlabel="t", ylabel="state")

# ╔═╡ 92cd25ba-520d-5fba-936b-5b37aa16bcc3
let μs = range(-1, 1; length=401)
    figure = plot(μs, zero.(μs); label="origin", xlabel="μ", ylabel="equilibrium x")
    positive = filter(>=(0), μs)
    plot!(figure, positive, sqrt.(positive); label="positive branch")
    plot!(figure, positive, -sqrt.(positive); label="negative branch")
end

# ╔═╡ Cell order:
# ╠═db4f14e7-32f9-5c87-93ee-f988cc5aacca
# ╟─62518ca0-5495-5907-93cd-6535924c39f1
# ╟─800d9c5c-bae0-5ed5-91bb-8fdc7fcade41
# ╠═9bd171c2-0dbd-5768-91f4-a5023265fbb0
# ╠═7bb14eaf-4b73-5036-a893-afbbcbcc932e
# ╠═d13525f4-c075-5f9b-a75d-e1a55ca1e303
# ╠═f00db1e4-6298-5fb1-b3fe-df8bc73d8365
# ╠═ba0270a1-ef90-586c-8682-fd104d5b00d1
# ╠═23ede884-9339-547d-8047-550719420738
# ╠═50fca1fa-bcab-5504-af22-0302e55ea7df
# ╠═92cd25ba-520d-5fba-936b-5b37aa16bcc3
