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
sol = solve(prob, Rosenbrock23(); abstol=1e-9, reltol=1e-7, saveat=0.01, dtmax=0.002);

# ╔═╡ b59acc55-95ad-5ea3-90fc-c9a167d656d7
plot(sol; idxs=[1, 2], xlabel="t", ylabel="state")

# ╔═╡ 6cb85a20-27f4-597f-a002-54a220d1ed11
plot(sol; idxs=(1, 2), legend=false)

# ╔═╡ ade92bc5-1af1-5d51-85c2-71448074ad98
md"""
## A different idealization: instantaneous impacts

The compliant model resolves contact with a spring. A rigid-floor model
instead integrates free flight and changes velocity at each downward crossing:

```math
\dot h=v,\quad \dot v=-g\quad(h>0),\qquad
h^+=0,\quad v^+=-e v^-\quad(h=0,\ v^-<0).
```

Here $0<e\leq1$ is restitution. The kinetic energy immediately after impact
is $e^2$ times that before impact. This is a **hybrid system**: a continuous
ODE plus a discrete reset. A root-locating callback finds the impact between
solver steps; it is not an extra contact force.
"""

# ╔═╡ d79bce35-20c2-5649-8eba-dec461c9f078
function course_flight!(du,u,p,t)
    g,e=p
    du[1]=u[2]
    du[2]=-g
    nothing
end

# ╔═╡ fa26fe16-877e-5013-b04b-d5af9c199526
course_restitution = 0.9  # restitution e

# ╔═╡ 55cd11db-ce2f-580a-a537-cca7fd509139
course_impact = ContinuousCallback(
    (u,t,integrator)->u[1], nothing,
    integrator->begin
        integrator.u[1]=0
        integrator.u[2]=-integrator.p[2]*integrator.u[2]
        abs(integrator.u[2])<1e-6 && terminate!(integrator)
    end;save_positions=(true,true))

# ╔═╡ 891df09d-fec3-5f95-a222-d6780dc6c07d
course_flight_u0 = [1.0,0.0]

# ╔═╡ c5ebe5f4-942f-5570-9001-24f6ddb8810a
course_flight_tspan = (0.0,5.0)

# ╔═╡ c92f0f9d-2c92-57fc-9b7b-a9f7f97c818a
course_flight_p = [9.8,course_restitution]

# ╔═╡ 5ac8105b-d8fa-596a-a98b-77c6f13e96d3
course_flight_prob = ODEProblem(course_flight!, course_flight_u0, course_flight_tspan, course_flight_p)

# ╔═╡ f1186a7c-6d8d-5e84-ba1d-7a62091d5301
course_flight_sol = solve(course_flight_prob, Tsit5(); abstol=1e-9, reltol=1e-8, callback=course_impact, dtmax=0.02);

# ╔═╡ e9216006-3993-501c-9912-54a015c42765
let heights=plot(course_flight_sol;idxs=(0,1),xlabel="t",ylabel="h",legend=false),
    phase=plot(course_flight_sol;idxs=(1,2),xlabel="h",ylabel="v",legend=false)
    plot(heights,phase;layout=(1,2),size=(900,370),margin=5*Plots.mm)
end

# ╔═╡ efe02735-3d53-5a9b-a0cb-190d09d573fd
md"""
**Try:** set $e=1$ and then $e=0.9$. Check the sequence of maximum heights:
$h_{n+1}/h_n=e^2$. The vertical jumps in the phase plane are instantaneous
resets, not ODE trajectories. Compare these with the short, continuous contact
arcs of the compliant model. A tiny-speed termination avoids accumulating
infinitely many impacts if you later use smaller restitution and longer times.
"""

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
# ╟─ade92bc5-1af1-5d51-85c2-71448074ad98
# ╠═d79bce35-20c2-5649-8eba-dec461c9f078
# ╠═fa26fe16-877e-5013-b04b-d5af9c199526
# ╠═55cd11db-ce2f-580a-a537-cca7fd509139
# ╠═891df09d-fec3-5f95-a222-d6780dc6c07d
# ╠═c5ebe5f4-942f-5570-9001-24f6ddb8810a
# ╠═c92f0f9d-2c92-57fc-9b7b-a9f7f97c818a
# ╠═5ac8105b-d8fa-596a-a98b-77c6f13e96d3
# ╠═f1186a7c-6d8d-5e84-ba1d-7a62091d5301
# ╠═e9216006-3993-501c-9912-54a015c42765
# ╟─efe02735-3d53-5a9b-a0cb-190d09d573fd
