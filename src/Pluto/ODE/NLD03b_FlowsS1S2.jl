### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 03c69f7b-1eeb-585d-a787-5f968b7eaa08
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots
end

# ╔═╡ 0e026a8c-edbd-521c-938f-421f82bb3ec1
md"""
# Circle dynamics and coupled phases

Integrate the unwrapped Adler phases and plot their sines to see whether the oscillators lock.
"""

# ╔═╡ 848b6ebe-0b70-5c80-a623-e4a08a858f1b
md"""
## From a fixed phase to a rotating phase

Adler's equation is a flow on a circle:

```math
\dot\theta=\omega-a\sin\theta,\qquad
T=\int_0^{2\pi}\frac{d\theta}{\omega-a\sin\theta}
  =\frac{2\pi}{\sqrt{\omega^2-a^2}}\quad(\omega>a\geq0).
```

For $0<\omega<a$ there are two fixed phases, one stable and one unstable.
At $\omega=a$ they collide on the circle; just above it the trajectory spends
most of each turn near the vanished pair. This is a **saddle-node on an
invariant circle** (SNIC, also called SNIPER). Its period diverges, unlike the
finite onset period of a generic Hopf bifurcation.

We integrate the unwrapped phase and wrap it modulo $2\pi$ only for display.
"""

# ╔═╡ 1e94ea76-b004-5707-8179-5d2325e5f49c
function course_adler!(du,u,p,t)
    ω,a=p
    du[1]=ω-a*sin(u[1])
    nothing
end

# ╔═╡ 7b07ca66-fa5a-5c97-9309-1513e602495c
course_omega = 1.05  # ω (a = 1)

# ╔═╡ 07d78039-374e-59c4-805f-cf4490f00b3f
course_adler_u0 = [0.0]

# ╔═╡ 8fd4ffb7-44cf-5cf1-b727-f6b83c08f18e
course_adler_tspan = (0.0,80.0)

# ╔═╡ 123e9200-8ab9-5bff-83b7-c524e073633e
course_adler_p = [course_omega,1.0]

# ╔═╡ ee536616-ef75-5151-a1c9-40ec1d64c63b
course_adler_prob = ODEProblem(course_adler!, course_adler_u0, course_adler_tspan, course_adler_p)

# ╔═╡ 58c9063e-5061-50dd-9012-6c40670dec7a
course_adler_sol = solve(course_adler_prob, Tsit5(); abstol=1e-9, reltol=1e-8);

# ╔═╡ 26b74603-d8f5-5391-8da9-299c2cd216bf
let ts=range(course_adler_tspan...;length=1601), θ=range(0,2pi;length=401)
    phases=course_adler_sol(ts;idxs=1).u
    a=plot(θ,course_omega .-sin.(θ);xlabel="θ",ylabel="dθ/dt",legend=false)
    hline!(a,[0];color=:gray)
    b=plot(ts,mod.(phases,2pi);xlabel="t",ylabel="θ mod 2π",legend=false)
    c=plot(ts,sin.(phases);xlabel="t",ylabel="sin θ",legend=false)
    plot(a,b,c;layout=(1,3),size=(1000,340),margin=5*Plots.mm)
end

# ╔═╡ ae6b22d5-4acd-5200-bb38-a52d2e8b8722
let frequencies=range(1.001,1.5;length=400)
    f=plot(frequencies,2pi ./sqrt.(frequencies.^2 .-1);
        xlabel="ω (a = 1)",ylabel="period T",label="exact rotating solution")
    if course_omega>1
        scatter!(f,[course_omega],[2pi/sqrt(course_omega^2-1)];label="selected ω")
    end
    f
end

# ╔═╡ 447b862a-6a18-542c-85ed-3785505a90b0
md"""
## Two coupled phases

Now add coupling and compare the two unwrapped phases.
"""

# ╔═╡ 1278bc7f-dbfe-5b07-906a-782c3dcc4ebf
md"""
```math
\begin{aligned}
\dot{\theta}_1 &= \omega_1-\sin\theta_1+\kappa\sin(\theta_2-\theta_1),\\
\dot{\theta}_2 &= \omega_2-\sin\theta_2+\kappa\sin(\theta_1-\theta_2).
\end{aligned}
```
"""

# ╔═╡ 1ffd3f2d-1e4f-5b37-8536-84b15e1277db
function model!(du, u, p, t)
    θ1, θ2 = u
    ω1, ω2, κ = p
    du[1] = ω1 - sin(θ1) + κ*sin(θ2-θ1)
    du[2] = ω2 - sin(θ2) + κ*sin(θ1-θ2)
    nothing
end

# ╔═╡ e45a6c49-f2f3-5598-a4d8-70c4504da728
u0 = [0.0, 0.2]

# ╔═╡ 1294d097-6879-58b3-98e8-20c5d7588479
tspan = (0.0, 80.0)

# ╔═╡ d155dc43-216e-59e4-9a14-6ab09f5ccc31
p = [1.2, 1.4, 0.3]  # ω1, ω2, κ

# ╔═╡ 2935661f-c3b9-5b21-b878-3ceca93d1491
prob = ODEProblem(model!, u0, tspan, p)

# ╔═╡ 42aba460-4bc1-508b-be55-db624b15da0d
sol = solve(prob, Tsit5(); abstol=1e-9, reltol=1e-7, saveat=0.02);

# ╔═╡ 0bbc70dd-3183-52b2-b5d5-d0438d473b56
plot(sol.t, [sin.(sol[1,:]) sin.(sol[2,:])];
    xlabel="t", ylabel="sin(θ)", label=["oscillator 1" "oscillator 2"])

# ╔═╡ fdf3fd3c-6b90-5302-89df-115bc59a08f6
md"""
**Try:** compare $\omega=0.95$, $1$, and $1.05$. Does an almost constant
phase necessarily mean a stable equilibrium? For the coupled phases, compare
the **unwrapped phase difference**: bounded variation is consistent with
1:1 locking; repeated $2\pi$ slips indicate drift. Two similar sine traces
over a short interval are not enough to establish locking.
"""

# ╔═╡ 158e8d19-ce87-5789-9c39-dadc6fc61cd8
let s=sol, ts=range(first(sol.t),last(sol.t);length=1201)
    delta=[s(t)[2]-s(t)[1] for t in ts]
    plot(ts,delta;xlabel="t",ylabel="θ₂ − θ₁ (unwrapped)",legend=false)
end

# ╔═╡ Cell order:
# ╠═03c69f7b-1eeb-585d-a787-5f968b7eaa08
# ╟─0e026a8c-edbd-521c-938f-421f82bb3ec1
# ╟─848b6ebe-0b70-5c80-a623-e4a08a858f1b
# ╠═1e94ea76-b004-5707-8179-5d2325e5f49c
# ╠═7b07ca66-fa5a-5c97-9309-1513e602495c
# ╠═07d78039-374e-59c4-805f-cf4490f00b3f
# ╠═8fd4ffb7-44cf-5cf1-b727-f6b83c08f18e
# ╠═123e9200-8ab9-5bff-83b7-c524e073633e
# ╠═ee536616-ef75-5151-a1c9-40ec1d64c63b
# ╠═58c9063e-5061-50dd-9012-6c40670dec7a
# ╠═26b74603-d8f5-5391-8da9-299c2cd216bf
# ╠═ae6b22d5-4acd-5200-bb38-a52d2e8b8722
# ╟─447b862a-6a18-542c-85ed-3785505a90b0
# ╟─1278bc7f-dbfe-5b07-906a-782c3dcc4ebf
# ╠═1ffd3f2d-1e4f-5b37-8536-84b15e1277db
# ╠═e45a6c49-f2f3-5598-a4d8-70c4504da728
# ╠═1294d097-6879-58b3-98e8-20c5d7588479
# ╠═d155dc43-216e-59e4-9a14-6ab09f5ccc31
# ╠═2935661f-c3b9-5b21-b878-3ceca93d1491
# ╠═42aba460-4bc1-508b-be55-db624b15da0d
# ╠═0bbc70dd-3183-52b2-b5d5-d0438d473b56
# ╟─fdf3fd3c-6b90-5302-89df-115bc59a08f6
# ╠═158e8d19-ce87-5789-9c39-dadc6fc61cd8
