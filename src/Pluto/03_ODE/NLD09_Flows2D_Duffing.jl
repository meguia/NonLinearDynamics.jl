### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ ff4ccfc6-c478-5eb8-ad55-6a263f056929
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots
end

# ╔═╡ cadba734-5b4b-510b-abd1-bddc2792cbec
md"""
# The Duffing oscillator

A double-well potential gives two stable resting positions; change the initial state to explore their basins.
"""

# ╔═╡ 68460e98-3998-537c-8a5f-be8c7b6b1331
md"""
```math
\begin{aligned}
\dot{x} &= v, & \dot{v} &= -\gamma v+\beta x-x^3.
\end{aligned}
```
"""

# ╔═╡ 4bc05396-18f8-59fb-be29-1f988fa376a9
function model!(du, u, p, t)
    x, v = u
    γ, β = p
    du[1] = v
    du[2] = -γ*v + β*x - x^3
    nothing
end

# ╔═╡ b9eb2d07-4791-520a-a3c8-6c2b1fbf55a4
u0 = [0.1, 0.7]

# ╔═╡ a1df769d-1131-5b35-86cc-2c969cb2d19c
tspan = (0.0, 60.0)

# ╔═╡ d1ba8bea-1d62-51f2-8044-4cf10095fae1
p = [0.15, 1.0]  # γ, β

# ╔═╡ c67fa65d-32d0-5cef-a854-33de4d66d906
prob = ODEProblem(model!, u0, tspan, p)

# ╔═╡ 9a436f19-7c6a-5fe4-a884-4d0e268b0f8b
sol = solve(prob, Tsit5(); abstol=1e-9, reltol=1e-7, saveat=0.02);

# ╔═╡ 3125de08-bb5b-5dee-b608-cda49b83d138
plot(sol; idxs=[1, 2], xlabel="t", ylabel="state")

# ╔═╡ a8bb9eb9-8562-5c7f-8656-ed6ab53c0dd3
plot(sol; idxs=(1, 2), legend=false)

# ╔═╡ Cell order:
# ╠═ff4ccfc6-c478-5eb8-ad55-6a263f056929
# ╟─cadba734-5b4b-510b-abd1-bddc2792cbec
# ╟─68460e98-3998-537c-8a5f-be8c7b6b1331
# ╠═4bc05396-18f8-59fb-be29-1f988fa376a9
# ╠═b9eb2d07-4791-520a-a3c8-6c2b1fbf55a4
# ╠═a1df769d-1131-5b35-86cc-2c969cb2d19c
# ╠═d1ba8bea-1d62-51f2-8044-4cf10095fae1
# ╠═c67fa65d-32d0-5cef-a854-33de4d66d906
# ╠═9a436f19-7c6a-5fe4-a884-4d0e268b0f8b
# ╠═3125de08-bb5b-5dee-b608-cda49b83d138
# ╠═a8bb9eb9-8562-5c7f-8656-ed6ab53c0dd3
