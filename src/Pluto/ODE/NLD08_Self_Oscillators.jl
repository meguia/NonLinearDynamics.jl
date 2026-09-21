### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 04c345c5-fd34-5830-97c4-edadbb05055c
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots, PlutoUI
end

# ╔═╡ a57d36fd-758f-4de7-ae18-29214e5c3565
TableOfContents()

# ╔═╡ f7505e1a-89cd-5006-b8f1-f42026ef3eda
md"""
# Self-oscillation and excitability

Negative damping can supply energy at small amplitude while nonlinear
dissipation limits large motion. Follow the circuit oscillator into its
slow–fast description, then compare reed, friction, and excitable systems.
"""

# ╔═╡ 9d1e4712-71cc-54c1-bbaf-99e1520e5274
md"""
## Van der Pol: from negative damping to relaxation

A nonlinear resistor in an RLC circuit motivates displacement-dependent
damping. In dimensionless variables the van der Pol oscillator is

```math
\dot x=v,\qquad \dot v=\mu(1-x^2)v-x,\qquad
\dot E=\mu(1-x^2)v^2,\quad E=(x^2+v^2)/2.
```

For $\mu>0$, energy enters when $|x|<1$ and leaves when $|x|>1$.
Increasing $\mu$ changes a nearly sinusoidal oscillation into slow drifts
separated by fast jumps. Compare this with Rayleigh's **velocity-dependent**
damping below.
"""

# ╔═╡ e6eff0eb-9700-59ec-b0c1-9d71a5e7eaa9
function course_vdp!(du,u,p,t)
    μ=only(p)
    x,v=u
    du[1]=v
    du[2]=μ*(1-x^2)*v-x
    nothing
end

# ╔═╡ e34cbc3b-4ceb-5a58-a2f8-42f9a87c6154
course_vdp_mu = 3.0  # μ

# ╔═╡ f14251bd-0578-5e4a-81bf-752d089a6620
course_vdp_u0 = [0.1,0.0]

# ╔═╡ f8e42bcf-fd2a-5b3d-9629-8f510a549e89
course_vdp_tspan = (0.0,80.0)

# ╔═╡ 818ef04d-2ff4-507b-8690-ed8331c23b57
course_vdp_p = [course_vdp_mu]

# ╔═╡ bc511430-4a94-5514-8a23-56ad852237fc
course_vdp_prob = ODEProblem(course_vdp!, course_vdp_u0, course_vdp_tspan, course_vdp_p)

# ╔═╡ 1ad45939-6c73-5879-bb5e-184115d89996
course_vdp_sol = solve(course_vdp_prob, Tsit5(); abstol=1e-9, reltol=1e-8);

# ╔═╡ 5b8f7715-dd61-54ab-9c85-1720ef9f9ee7
let a=plot(course_vdp_sol;idxs=(0,1),xlabel="t",ylabel="x",legend=false),
    b=plot(course_vdp_sol;idxs=(1,2),xlabel="x",ylabel="v",legend=false)
    plot(a,b;layout=(1,2),size=(900,370),margin=5*Plots.mm)
end

# ╔═╡ a885d36c-998e-5655-bb34-050a5aa40db7
md"""
### The same oscillator in Liénard coordinates

For $\mu>0$, set $z=x-x^3/3-v/\mu$. Differentiation gives

```math
\dot x=\mu(x-x^3/3-z),\qquad \dot z=x/\mu.
```

The cubic $z=x-x^3/3$ is the $\dot x=0$ nullcline; $x=0$ is the
$\dot z=0$ nullcline. Away from the cubic, $x$ changes quickly when $\mu$
is large. The motion slows near its outer branches and jumps near the folds.
The coordinate transformation changes the picture, not the underlying trajectory.
"""

# ╔═╡ feb3747f-f9bd-5f1c-90d3-ef67c9a520d1
function course_lienard!(du,u,p,t)
    μ=only(p)
    x,z=u
    du[1]=μ*(x-x^3/3-z)
    du[2]=x/μ
    nothing
end

# ╔═╡ 82fefbf8-76a0-5110-bfd4-71ca3c54a9e9
course_lienard_u0 = [course_vdp_u0[1],course_vdp_u0[1]-course_vdp_u0[1]^3/3-course_vdp_u0[2]/course_vdp_mu]

# ╔═╡ a84faf77-7103-5030-9fa3-01ae932b2f3b
course_lienard_tspan = course_vdp_tspan

# ╔═╡ 3e5f1e47-16a1-5b7e-b529-3a74c941bde6
course_lienard_p = course_vdp_p

# ╔═╡ 2ac8a1bb-710d-50e1-a77e-291dceaed81d
course_lienard_prob = ODEProblem(course_lienard!, course_lienard_u0, course_lienard_tspan, course_lienard_p)

# ╔═╡ 6d50fecd-cee6-5bf8-8a78-c2e44ec027b9
course_lienard_sol = solve(course_lienard_prob, Tsit5(); abstol=1e-9, reltol=1e-8);

# ╔═╡ 7246f732-7143-549b-8f4c-84a3ca01cf84
let ts=range(course_vdp_tspan...;length=2001), xs=range(-2.3,2.3;length=401)
    original=Array(course_vdp_sol(ts))
    transformed=original[1,:] .-original[1,:].^3 ./3 .-original[2,:]./course_vdp_mu
    a=plot(original[1,:],transformed;label="transformed (x,v)",xlabel="x",ylabel="z")
    plot!(a,course_lienard_sol;idxs=(1,2),linestyle=:dash,label="solve (x,z)")
    plot!(a,xs,xs .-xs.^3 ./3;label="dx/dt = 0",color=:red)
    vline!(a,[0];label="dz/dt = 0",color=:green)
    b=plot(ts,abs.(original[1,:] .-course_lienard_sol(ts;idxs=1).u);
        xlabel="t",ylabel="|x_original − x_Liénard|",legend=false)
    plot(a,b;layout=(1,2),size=(950,390),margin=5*Plots.mm)
end

# ╔═╡ ebabfb23-459f-5d5a-87f8-e9628d5d59cf
md"""
**Try:** compare $\mu=0.2$, $3$, and $8$. Where does the motion spend most of
its time? In the standard van der Pol family, $\mu=0$ removes the nonlinear
damping entirely and leaves a center. It is a degenerate onset, not the
generic small-amplitude Hopf bifurcation studied in NLD12; the limit-cycle
amplitude tends to about $2$ as $\mu\to0^+$.
"""

# ╔═╡ 79124cae-eeb5-5104-9c14-f317202308dc
md"""
## Self-oscillator: simple reed model (Rayleigh)

Rayleigh's simple model of a blown reed replaces the constant damping of a
harmonic oscillator with a velocity-dependent dissipation coefficient:

```math
\gamma(v)=s_c^2v^2-\mu.
```

For $\mu>0$, this coefficient is negative at small velocities, so the damping
force feeds energy into the motion. At large velocities it becomes positive
and removes energy, limiting the amplitude.

Here $x$ is displacement, $v$ is velocity, $K>0$ sets the restoring force,
$\mu$ controls the excitation, and $s_c>0$ sets the nonlinear damping scale.
The choice $s_c=1$ gives the Rayleigh oscillator introduced in the 2D flows notebook.
"""

# ╔═╡ f7e71652-3503-5e50-9997-66693bd59334
md"""
The first-order equations are

```math
\begin{aligned}
\dot{x} &= v,\\
\dot{v} &= -Kx-\gamma(v)v=-Kx+(\mu-s_c^2v^2)v.
\end{aligned}
```

The cubic term $-s_c^2v^3$ makes the system nonlinear. The function below writes
these two derivatives into `du`, with `u = [x, v]` and `p = [μ, K, sc]`.
"""

# ╔═╡ 0a72128d-f117-5d09-8fac-0050ff3d6dea
function model!(du, u, p, t)
    x, v = u
    μ, K, sc = p
    du[1] = v
    du[2] = -K*x + (μ-sc^2*v^2)*v
    nothing
end

# ╔═╡ e1e36c53-069d-4c5b-958a-a80706a94015
md"""
Set the initial condition, time span, and parameters, then build and solve the
problem with `Tsit5()`. Edit these values directly to explore the system; this
longer time span shows the approach to sustained oscillation.
"""

# ╔═╡ 30daf2b7-9d99-596e-a0fe-822e9a7fcf59
u0 = [0.1, 0.0]  # x(0), v(0)

# ╔═╡ bf02afe1-302b-5060-abd2-0bf1eac1a3a3
tspan = (0.0, 100.0)

# ╔═╡ dfa42413-a9f5-5a28-b603-3e33e473504f
p = [0.3, 1.0, 1.0]  # μ, K, sc

# ╔═╡ 23966858-70a2-5835-bd94-373f629d0b02
prob = ODEProblem(model!, u0, tspan, p)

# ╔═╡ c4de1520-29c1-573a-a809-1d1f9c0b155b
sol = solve(prob, Tsit5(); abstol=1e-9, reltol=1e-7, saveat=0.02);

# ╔═╡ d2c7d6f2-92c1-4df7-8c93-cd01914c8d96
md"""
Solve the same equations from the opposite initial condition. The two
trajectories approach the same limit cycle with different phases.
"""

# ╔═╡ b3d660bc-4c81-4f85-a0f2-c229a6e4b231
opposite_u0 = -u0

# ╔═╡ 2f161ca1-ea43-42f6-86d1-b785195124ed
opposite_prob = ODEProblem(model!, opposite_u0, tspan, p)

# ╔═╡ 23f9eea3-2a5b-4b5a-a88f-9f5112751284
opposite_sol = solve(opposite_prob, Tsit5(); abstol=1e-9, reltol=1e-7, saveat=0.02);

# ╔═╡ 9bad4688-3b4c-55d1-8e8e-eedcc90151a5
begin
    reed_phase = plot(sol; idxs=(1, 2), label="u₀", xlabel="x", ylabel="v",
        linewidth=1.5, aspect_ratio=1, legend=:outerright,
        title="Rayleigh reed", size=(650, 450), margin=5 * Plots.mm)
    plot!(reed_phase, opposite_sol; idxs=(1, 2), label="−u₀", linewidth=1.5)
    scatter!(reed_phase, [u0[1], opposite_u0[1]], [u0[2], opposite_u0[2]];
        color=:black, markersize=3, label="initial conditions")
end

# ╔═╡ c53d29ae-8f69-4cab-9d14-45a66a923225
md"""
The vector field shows the direction of motion at each point in the phase
plane. Evaluate `model!` on a grid and use a common scale for the arrows so
their relative lengths still represent speed.
"""

# ╔═╡ ac399f96-f1c8-4b37-a502-fc1d40d5b329
begin
    x_grid = range(-1.0, 1.0; length=21)
    v_grid = range(-1.0, 1.0; length=21)
    field_states = [[x, v] for x in x_grid for v in v_grid]
    field_derivatives = [zeros(2) for state in field_states]
    for (du, state) in zip(field_derivatives, field_states)
        model!(du, state, p, 0.0)
    end
    arrow_scale = 0.08 / maximum(hypot(du...) for du in field_derivatives)
    reed_field = quiver(first.(field_states), last.(field_states);
        quiver=(arrow_scale .* first.(field_derivatives),
                arrow_scale .* last.(field_derivatives)),
        color=:steelblue, alpha=0.6, legend=false,
        xlabel="x", ylabel="v", xlims=(-1, 1), ylims=(-1, 1),
        aspect_ratio=1, size=(600, 600), title="Rayleigh reed: vector field")
    plot!(reed_field, sol; idxs=(1, 2), color=:black, linewidth=1.5)
end

# ╔═╡ 2b0bb6b3-393f-4cce-a1d3-8ad67ea2a596
md"""
The time traces show the initial growth and the eventual bounded oscillation.
For the oscillator energy,

```math
E=\frac{Kx^2+v^2}{2},\qquad
\dot{E}=\mu v^2-s_c^2v^4.
```

Energy input and dissipation balance over one period of the limit cycle.
Try a negative $\mu$: energy then decreases and the motion decays to rest.
"""

# ╔═╡ af23b3c0-fd6f-5eac-976b-ee99eb9ef52a
begin
    reed_displacement = plot(sol; idxs=(0, 1), legend=false, xlabel="t", ylabel="x")
    reed_velocity = plot(sol; idxs=(0, 2), legend=false, xlabel="t", ylabel="v")
    hline!(reed_displacement, [0.0]; color=:black, linewidth=0.5)
    hline!(reed_velocity, [0.0]; color=:black, linewidth=0.5)
    reed_timeseries = plot(reed_displacement, reed_velocity;
        layout=(2, 1), size=(900, 400), margin=5 * Plots.mm)
end

# ╔═╡ 0323f914-6959-42f0-9849-74b301ef8366
md"""
## Rubbed oscillator (bowed string)

A bowed string is another example of self-oscillation. A simple mechanical
analogy is a mass attached to a spring and resting on a moving conveyor belt:
the belt supplies energy through friction, just as a bow drives a string.
"""

# ╔═╡ e6c2a8a9-92c7-45d5-9341-eec3ca7824da
html"""
<div>
<img src="https://i.imgur.com/qW4INmr.png" width="300px"
     alt="A mass attached to a spring rests on a conveyor belt moving to the right.">
</div>
"""

# ╔═╡ 429f3efd-0809-4ec4-bcf1-2ad35024f9f8
md"""
At first, static friction keeps the mass moving with the belt while the spring
stretches. When the required force exceeds the maximum static friction, the
mass slips and the spring pulls it back. The mass can stick again when its
velocity matches the belt and the required friction is below the static
limit. Repetition of these stages produces **stick–slip motion**.

An idealized friction law distinguishes sticking from sliding. At zero slip,
static friction can take a range of values up to a threshold; during sliding,
the friction changes sign with the slip and its magnitude can decrease as
sliding becomes faster.
"""

# ╔═╡ 3ffb5784-e3f7-448a-bf3e-8a1c68c7c9fb
html"""
<div>
<img src="https://i.imgur.com/KrRu2Ub.png" width="200px"
     alt="Idealized friction characteristic: a static-friction interval at zero slip and decreasing sliding friction on either side.">
</div>
"""

# ╔═╡ 3d61c9a1-8c5b-4d98-a641-e4a8e22b2f62
md"""
The horizontal axis represents the **slip**, the difference between the mass
velocity and the belt velocity. Below we write it as $s=v-V$, where
$v=\dot{x}$ and $V$ is the constant belt or bow speed (the sketch labels that
speed $v$).

For a smooth ODE, we use the cubic approximation from the interactive notebook:

```math
C(s)=-s(1-s^2).
```

Its negative slope near zero slip provides excitation, while the cubic term
limits large velocities. This smooth law does not reproduce the discontinuity
above or impose exact sticking; a true stick–slip model needs a separate
static-friction rule.
"""

# ╔═╡ c15c3c18-42a6-47a3-9374-8a5785019220
friction(s) = -s * (1 - s^2)

# ╔═╡ af1580b7-6041-4fcb-96e3-d8f147b3d8b4
begin
    slip_grid = range(-1.5, 1.5; length=401)
    friction_plot = plot(slip_grid, friction.(slip_grid);
        xlabel="s = v − V", ylabel="C(s)", legend=false,
        title="Smooth cubic friction law", size=(600, 300), margin=5 * Plots.mm)
    hline!(friction_plot, [0.0]; color=:black, linewidth=0.5)
    vline!(friction_plot, [0.0]; color=:black, linewidth=0.5)
end

# ╔═╡ 6b395512-6c72-4ad2-b658-4c144cb05a82
md"""
With unit mass and spring stiffness, the equations are

```math
\begin{aligned}
\dot{x} &= v,\\
\dot{v} &= -x-\mu C(v-V).
\end{aligned}
```

Here $\mu>0$ scales the friction term, and $V$ sets the bow speed. As before,
`u = [x, v]`; this time the parameter vector is `p = [μ, V]`.
"""

# ╔═╡ 5d099649-1a70-4855-883c-0919985ff116
function bow!(du, u, p, t)
    x, v = u
    μ, V = p
    du[1] = v
    du[2] = -x - μ * friction(v - V)
    nothing
end

# ╔═╡ 6e985be6-b502-47cd-b61a-27a4b6d817c0
bow_u0 = [0.7, 0.0]  # x(0), v(0)

# ╔═╡ 27f839e9-8c58-4d3b-a31b-5eb367293d85
bow_tspan = (0.0, 250.0)

# ╔═╡ d3f29269-e630-4ee5-aeae-4f1e6e9427d6
bow_p = [0.02, 0.02]  # μ, V

# ╔═╡ 014ecf1e-3638-4ac9-ae3e-21218d0044e1
bow_prob = ODEProblem(bow!, bow_u0, bow_tspan, bow_p)

# ╔═╡ 5c797db7-2d6c-4e3d-b2ce-fcfd27f5e825
bow_sol = solve(bow_prob, Tsit5(); abstol=1e-9, reltol=1e-7, saveat=0.02);

# ╔═╡ bd085ecc-1f11-40b9-82b6-9f8177c59d37
md"""
The time traces show the transient and the sustained oscillation, while the
phase portrait approaches a closed orbit. Edit `bow_p` to explore how the
excitation strength and bow speed affect the motion.
"""

# ╔═╡ 96992066-915d-4b06-b07f-08d5ba8cb98c
begin
    bow_timeseries = plot(bow_sol; idxs=[1, 2], label=["x" "v"],
        xlabel="t", ylabel="state", title="Bowed string: time traces")
    bow_phase = plot(bow_sol; idxs=(1, 2), legend=false,
        xlabel="x", ylabel="v", aspect_ratio=1, arrow=true,
        title="Bowed string: phase portrait")
    bow_figure = plot(bow_timeseries, bow_phase;
        layout=(1, 2), size=(900, 450), margin=5 * Plots.mm)
end

# ╔═╡ 8a9c2b54-3e5b-532e-b367-4432f490b2ac
md"""
## Excitability and self-oscillation: FitzHugh–Nagumo

A fast activation variable $x$ and a slower recovery variable $y$ provide
a simple neuron analogy:

```math
\dot x=x-x^3/3-y+I,\qquad \dot y=(ax+b-y)/10.
```

The constant input $I$ shifts the cubic nullcline vertically; the recovery
nullcline is $y=ax+b$. A displaced stable equilibrium may generate one large
excursion before returning (excitability). With different parameters, an
attracting cycle gives repeated excursions. The factor $10$ separates time scales.
"""

# ╔═╡ c3a7a06c-a001-5f42-87ae-ce35c38d40f0
function course_fhn!(du,u,p,t)
    a,b,I=p
    x,y=u
    du[1]=x-x^3/3-y+I
    du[2]=(a*x+b-y)/10
    nothing
end

# ╔═╡ 43bb3e59-f8d5-53ef-8092-0a0d91c4d94f
course_current = 0.5  # I

# ╔═╡ fe44aef9-4ef4-55bb-8ea3-b6b74770e4fb
course_fhn_u0 = [-1.0,-0.5]

# ╔═╡ f2593d11-26a9-5931-b6e5-02188ca73616
course_fhn_tspan = (0.0,200.0)

# ╔═╡ a7232533-dc64-5b4c-89f1-45c660ae83e9
course_fhn_p = [1.25,0.875,course_current]

# ╔═╡ 57487c27-4608-546b-bd7f-704a3f268685
course_fhn_prob = ODEProblem(course_fhn!, course_fhn_u0, course_fhn_tspan, course_fhn_p)

# ╔═╡ c9de10b7-aea5-5846-9cd2-b5bee9d84f69
course_fhn_sol = solve(course_fhn_prob, Tsit5(); abstol=1e-9, reltol=1e-8);

# ╔═╡ 8109211d-c989-5b68-a299-473d28bd7fd9
let xs=range(-2.5,2.5;length=401)
    a,b,I=course_fhn_p
    phase=plot(course_fhn_sol;idxs=(1,2),label="trajectory",xlabel="x",ylabel="y")
    plot!(phase,xs,xs .-xs.^3 ./3 .+I;label="dx/dt = 0",color=:red)
    plot!(phase,xs,a.*xs .+b;label="dy/dt = 0",color=:green)
    trace=plot(course_fhn_sol;idxs=(0,1),xlabel="t",ylabel="activation x",legend=false)
    plot(phase,trace;layout=(1,2),size=(950,390),margin=5*Plots.mm)
end

# ╔═╡ b6827f62-16f2-5744-bc8b-948af517dc26
md"""
**Try:** compare $I=0$ and $I=0.5$, and perturb the initial activation.
Distinguish a long transient from repeated oscillation. For this family,

```math
J=\begin{pmatrix}1-x_*^2&-1\\a/10&-1/10\end{pmatrix},\quad
\operatorname{tr}J=0.9-x_*^2,\quad
\det J=(x_*^2+a-1)/10.
```

A possible Hopf onset requires zero trace **and positive determinant**,
as well as the usual nondegeneracy conditions. A zero trace alone is insufficient.
"""

# ╔═╡ Cell order:
# ╠═04c345c5-fd34-5830-97c4-edadbb05055c
# ╟─a57d36fd-758f-4de7-ae18-29214e5c3565
# ╟─f7505e1a-89cd-5006-b8f1-f42026ef3eda
# ╟─9d1e4712-71cc-54c1-bbaf-99e1520e5274
# ╠═e6eff0eb-9700-59ec-b0c1-9d71a5e7eaa9
# ╠═e34cbc3b-4ceb-5a58-a2f8-42f9a87c6154
# ╠═f14251bd-0578-5e4a-81bf-752d089a6620
# ╠═f8e42bcf-fd2a-5b3d-9629-8f510a549e89
# ╠═818ef04d-2ff4-507b-8690-ed8331c23b57
# ╠═bc511430-4a94-5514-8a23-56ad852237fc
# ╠═1ad45939-6c73-5879-bb5e-184115d89996
# ╠═5b8f7715-dd61-54ab-9c85-1720ef9f9ee7
# ╟─a885d36c-998e-5655-bb34-050a5aa40db7
# ╠═feb3747f-f9bd-5f1c-90d3-ef67c9a520d1
# ╠═82fefbf8-76a0-5110-bfd4-71ca3c54a9e9
# ╠═a84faf77-7103-5030-9fa3-01ae932b2f3b
# ╠═3e5f1e47-16a1-5b7e-b529-3a74c941bde6
# ╠═2ac8a1bb-710d-50e1-a77e-291dceaed81d
# ╠═6d50fecd-cee6-5bf8-8a78-c2e44ec027b9
# ╠═7246f732-7143-549b-8f4c-84a3ca01cf84
# ╟─ebabfb23-459f-5d5a-87f8-e9628d5d59cf
# ╟─79124cae-eeb5-5104-9c14-f317202308dc
# ╟─f7e71652-3503-5e50-9997-66693bd59334
# ╠═0a72128d-f117-5d09-8fac-0050ff3d6dea
# ╟─e1e36c53-069d-4c5b-958a-a80706a94015
# ╠═30daf2b7-9d99-596e-a0fe-822e9a7fcf59
# ╠═bf02afe1-302b-5060-abd2-0bf1eac1a3a3
# ╠═dfa42413-a9f5-5a28-b603-3e33e473504f
# ╠═23966858-70a2-5835-bd94-373f629d0b02
# ╠═c4de1520-29c1-573a-a809-1d1f9c0b155b
# ╟─d2c7d6f2-92c1-4df7-8c93-cd01914c8d96
# ╠═b3d660bc-4c81-4f85-a0f2-c229a6e4b231
# ╠═2f161ca1-ea43-42f6-86d1-b785195124ed
# ╠═23f9eea3-2a5b-4b5a-a88f-9f5112751284
# ╠═9bad4688-3b4c-55d1-8e8e-eedcc90151a5
# ╟─c53d29ae-8f69-4cab-9d14-45a66a923225
# ╠═ac399f96-f1c8-4b37-a502-fc1d40d5b329
# ╟─2b0bb6b3-393f-4cce-a1d3-8ad67ea2a596
# ╠═af23b3c0-fd6f-5eac-976b-ee99eb9ef52a
# ╟─0323f914-6959-42f0-9849-74b301ef8366
# ╟─e6c2a8a9-92c7-45d5-9341-eec3ca7824da
# ╟─429f3efd-0809-4ec4-bcf1-2ad35024f9f8
# ╟─3ffb5784-e3f7-448a-bf3e-8a1c68c7c9fb
# ╟─3d61c9a1-8c5b-4d98-a641-e4a8e22b2f62
# ╠═c15c3c18-42a6-47a3-9374-8a5785019220
# ╠═af1580b7-6041-4fcb-96e3-d8f147b3d8b4
# ╟─6b395512-6c72-4ad2-b658-4c144cb05a82
# ╠═5d099649-1a70-4855-883c-0919985ff116
# ╠═6e985be6-b502-47cd-b61a-27a4b6d817c0
# ╠═27f839e9-8c58-4d3b-a31b-5eb367293d85
# ╠═d3f29269-e630-4ee5-aeae-4f1e6e9427d6
# ╠═014ecf1e-3638-4ac9-ae3e-21218d0044e1
# ╠═5c797db7-2d6c-4e3d-b2ce-fcfd27f5e825
# ╟─bd085ecc-1f11-40b9-82b6-9f8177c59d37
# ╠═96992066-915d-4b06-b07f-08d5ba8cb98c
# ╟─8a9c2b54-3e5b-532e-b367-4432f490b2ac
# ╠═c3a7a06c-a001-5f42-87ae-ce35c38d40f0
# ╠═43bb3e59-f8d5-53ef-8092-0a0d91c4d94f
# ╠═fe44aef9-4ef4-55bb-8ea3-b6b74770e4fb
# ╠═f2593d11-26a9-5931-b6e5-02188ca73616
# ╠═a7232533-dc64-5b4c-89f1-45c660ae83e9
# ╠═57487c27-4608-546b-bd7f-704a3f268685
# ╠═c9de10b7-aea5-5846-9cd2-b5bee9d84f69
# ╠═8109211d-c989-5b68-a299-473d28bd7fd9
# ╟─b6827f62-16f2-5744-bc8b-948af517dc26
