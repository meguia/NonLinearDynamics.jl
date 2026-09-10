### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 05aebc47-b68a-528f-91fa-606bcd7b4d7f
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots
end

# ╔═╡ a0b70844-c2c8-55c1-a948-624f24633d00
md"""
# Bouncing ball: compliant contact

A spring acts only below the floor; its finite stiffness allows a small penetration during each bounce.
"""

# ╔═╡ 2b8a2280-2dd8-5d4b-924c-f9b1c09a5e16
md"""
```math
\begin{aligned}
\dot{h} &= v, & \dot{v} &= -g-K\min(h,0)-\gamma v.
\end{aligned}
```
"""

# ╔═╡ 89a292c1-f23c-584d-bfc7-2eea5984c3da
function model!(du, u, p, t)
    h, v = u
    g, K, γ = p
    du[1] = v
    du[2] = -g - K*min(h, 0) - γ*v
    nothing
end

# ╔═╡ 5be6f146-d797-54ac-893b-71e830d53112
u0 = [1.0, 0.0]

# ╔═╡ f13bd655-2c2a-5f7b-a5fc-a1bdac9501ce
tspan = (0.0, 5.0)

# ╔═╡ 4b8258d4-f5e6-5887-90f7-a82a2a1ec3b7
p = [9.8, 1000.0, 0.2]  # g, K, γ

# ╔═╡ 435a180e-09cb-5a88-a3f1-001095d89e08
prob = ODEProblem(model!, u0, tspan, p)

# ╔═╡ ec99df3d-37aa-5bd2-b79f-874011816030
sol = solve(prob, Rosenbrock23(); abstol=1e-9, reltol=1e-7, saveat=0.02, dtmax=0.002);

# ╔═╡ b59acc55-95ad-5ea3-90fc-c9a167d656d7
plot(sol; idxs=[1, 2], xlabel="t", ylabel="state")

# ╔═╡ 6cb85a20-27f4-597f-a002-54a220d1ed11
plot(sol; idxs=(1, 2), legend=false)

# ╔═╡ Cell order:
# ╠═05aebc47-b68a-528f-91fa-606bcd7b4d7f
# ╟─a0b70844-c2c8-55c1-a948-624f24633d00
# ╟─2b8a2280-2dd8-5d4b-924c-f9b1c09a5e16
# ╠═89a292c1-f23c-584d-bfc7-2eea5984c3da
# ╠═5be6f146-d797-54ac-893b-71e830d53112
# ╠═f13bd655-2c2a-5f7b-a5fc-a1bdac9501ce
# ╠═4b8258d4-f5e6-5887-90f7-a82a2a1ec3b7
# ╠═435a180e-09cb-5a88-a3f1-001095d89e08
# ╠═ec99df3d-37aa-5bd2-b79f-874011816030
# ╠═b59acc55-95ad-5ea3-90fc-c9a167d656d7
# ╠═6cb85a20-27f4-597f-a002-54a220d1ed11
