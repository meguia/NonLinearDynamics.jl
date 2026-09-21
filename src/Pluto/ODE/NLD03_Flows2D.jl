### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ c54bb4b0-bec2-51bb-85b2-8a3590b080a4
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using NonLinearDynamics: flow2d_nullclines
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

# ╔═╡ 15c9ea85-4382-5ee4-b938-4cda1928369d
md"""
## Nullclines

For the damped oscillator, $\dot x=0$ on $v=0$ (red) and $\dot v=0$
on $v=-kx/\gamma$ (green, when $\gamma>0$). On the red line the arrows are
vertical; on the green line they are horizontal. Only their intersection is
an equilibrium. A trajectory generally **crosses** a nullcline.

The phase portrait records motion in state space; the time trace also records
how quickly it moves. Predict the turning points of $x(t)$ from the crossings
of $v=0$ before inspecting the trace.
"""

# ╔═╡ 9814d360-a81f-5d3a-83a0-bd5db9cb8e7c
let parameters=p
    f=flow2d_nullclines(model!,parameters;xlims=[-1.5,1.5],ylims=[-1.5,1.5],
        regions=false,vectorfield=true,title="Nullclines and a trajectory")
    plot!(f,sol;idxs=(1,2),color=:black,linewidth=2)
    scatter!(f,[0],[0];color=:black,markersize=4)
end

# ╔═╡ 3a0ab20e-b1e2-54d7-83ab-c8c258e7c2ef
md"""
## Why damping removes the closed orbits

For $k>0$, $E=(kx^2+v^2)/2$ measures the oscillator's energy and
$\dot E=-\gamma v^2$. At zero damping, different initial conditions lie
on different constant-energy ellipses. At positive damping, energy decreases
and the origin attracts. A closed orbit of the undamped system is a **center
orbit**, not an attracting limit cycle.

"""

# ╔═╡ 6ebff00f-d4a0-526e-a8ab-c5967eb585a1
let parameters=p, s=sol
    energy=(parameters[1].*s[1,:].^2 .+s[2,:].^2)./2
    plot(s.t,energy;xlabel="t",ylabel="E",label="oscillator energy")
end

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
# ╟─15c9ea85-4382-5ee4-b938-4cda1928369d
# ╠═9814d360-a81f-5d3a-83a0-bd5db9cb8e7c
# ╟─3a0ab20e-b1e2-54d7-83ab-c8c258e7c2ef
# ╠═6ebff00f-d4a0-526e-a8ab-c5967eb585a1
