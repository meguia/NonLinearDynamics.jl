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

# ╔═╡ c4bd2050-9f00-51fa-8723-78b9fe12e8ff
md"""
## Compare a time trace, a stroboscopic section, and recurrence

Sample at $t_n=nT$, with $T=2\pi/\omega$, after discarding a transient.
One point suggests forcing-period locking; a finite set suggests a multiple
of that period. A smooth closed curve is consistent with quasiperiodicity.
A more complicated set calls for additional checks before claiming chaos.

The same equations above are used for the comparison below, with
$\gamma=0.14$, $\beta=1$, $\omega=1$, and adjustable forcing amplitude.
Try $A=0.1$ and $A=0.27$. The transient and the displayed segment are solved
on one continuous time interval, preserving the forcing phase.
"""

# ╔═╡ 7d88428f-b092-56a3-8e37-077b15f88e9f
course_forcing_A = 0.27  # forcing amplitude A

# ╔═╡ f966ba37-7f28-5ce6-8bf6-2c771489ff61
course_forcing_u0 = [0.5,0.5]

# ╔═╡ 19be70d4-b519-539b-8a2c-a2269e5deabe
course_forcing_p = [0.14,1.0,course_forcing_A,1.0]

# ╔═╡ 9ad1daf9-64f1-59e4-94bd-69cdc829d0ed
course_forcing_period = 2pi/course_forcing_p[4]

# ╔═╡ 09c2eb3c-9850-5d87-a148-cdb1947278d5
course_forcing_tspan = (0.0,160course_forcing_period)

# ╔═╡ 25662eae-f1e9-5018-9442-e7908e6b129d
course_forcing_prob = ODEProblem(model!, course_forcing_u0, course_forcing_tspan, course_forcing_p)

# ╔═╡ 9812d1c0-0c73-510b-961e-f76e911189ab
course_forcing_sol = solve(course_forcing_prob, Tsit5(); abstol=1e-9, reltol=1e-8, dtmax=course_forcing_period/40);

# ╔═╡ 6c4f7788-447b-5100-9648-87a25f4db1fc
let period=course_forcing_period
    section_times=(60:160).*period
    points=Array(course_forcing_sol(section_times;idxs=[1,2]))
    trace=plot(course_forcing_sol;idxs=(0,1),tspan=(140period,160period),
        xlabel="t",ylabel="x",legend=false)
    section=scatter(points[1,:],points[2,:];xlabel="x(nT)",ylabel="v(nT)",
        markersize=2,markerstrokewidth=0,legend=false)
    plot(trace,section;layout=(1,2),size=(950,380),margin=5*Plots.mm)
end

# ╔═╡ 466fc6b7-095a-5c28-9e13-90a771d4c552
md"""
A recurrence plot asks when the system returns near a previous **full state**:

```math
R_{ij}=\begin{cases}1&\|\mathbf z_i-\mathbf z_j\|\leq\varepsilon,\\0&\text{otherwise},
\end{cases}
\qquad
\mathbf z(t)=(x(t),v(t),\cos\omega t,\sin\omega t).
```

The two phase coordinates identify the forcing state without a discontinuity
at $2\pi$. Comparing only $(x,v)$ could mark a return at a different forcing
phase. All coordinates here are dimensionless with unit scales; changing
scales or the threshold changes the notion of “near”.

Repeated diagonal bands indicate repeated evolution; broken bands indicate
less regular recurrence. Neither a recurrence image nor an irregular time
trace alone establishes chaos.
"""

# ╔═╡ 31a70409-d1cb-54f3-8b12-8fc7b3912d0b
course_epsilon = 0.2  # recurrence threshold ε

# ╔═╡ 8c7b3170-2501-53ab-b9ac-f74feec0d52b
course_recurrence_times = collect(range(120course_forcing_period,140course_forcing_period;length=401))

# ╔═╡ 8886e5c2-5f09-5a80-b42e-1953b603222f
course_recurrence_states = vcat(
    Array(course_forcing_sol(course_recurrence_times;idxs=[1,2])),
    permutedims(cos.(course_forcing_p[4].*course_recurrence_times)),
    permutedims(sin.(course_forcing_p[4].*course_recurrence_times)))

# ╔═╡ 668f8c74-0f95-5a83-8e11-f5050765e97d
let states=course_recurrence_states, times=course_recurrence_times
    # Each column is [x, v, cos(phase), sin(phase)] at one sample time.
    n = size(states, 2)
    R = [sum(abs2, states[:,i] - states[:,j]) <= course_epsilon^2
         for i in 1:n, j in 1:n]
    heatmap(times, times, Int.(R); color=cgrad([:white, :black]), clims=(0,1),
        colorbar=false, xlabel="tᵢ", ylabel="tⱼ", title="Recurrence: distance ≤ ε",
        aspect_ratio=1, size=(550,500), margin=5*Plots.mm)
end

# ╔═╡ 5d4a13f9-afe7-5085-89fd-8a5992659a69
md"""
**Try:** compare the two amplitudes with the same threshold, sampling window,
and discarded transient. Then change only $\varepsilon$. Which features
persist? Return to NLD01 and explain why a periodic orbit of a flow can
be studied as a fixed or periodic point of a map.
"""

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
# ╟─c4bd2050-9f00-51fa-8723-78b9fe12e8ff
# ╠═7d88428f-b092-56a3-8e37-077b15f88e9f
# ╠═f966ba37-7f28-5ce6-8bf6-2c771489ff61
# ╠═19be70d4-b519-539b-8a2c-a2269e5deabe
# ╠═9ad1daf9-64f1-59e4-94bd-69cdc829d0ed
# ╠═09c2eb3c-9850-5d87-a148-cdb1947278d5
# ╠═25662eae-f1e9-5018-9442-e7908e6b129d
# ╠═9812d1c0-0c73-510b-961e-f76e911189ab
# ╠═6c4f7788-447b-5100-9648-87a25f4db1fc
# ╟─466fc6b7-095a-5c28-9e13-90a771d4c552
# ╠═31a70409-d1cb-54f3-8b12-8fc7b3912d0b
# ╠═8c7b3170-2501-53ab-b9ac-f74feec0d52b
# ╠═8886e5c2-5f09-5a80-b42e-1953b603222f
# ╠═668f8c74-0f95-5a83-8e11-f5050765e97d
# ╟─5d4a13f9-afe7-5085-89fd-8a5992659a69
