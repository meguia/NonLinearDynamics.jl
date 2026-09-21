### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 4df18ed1-628e-546a-a9f0-44a7912b76d0
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots
end

# ╔═╡ f1b40b34-ca27-5718-9166-0958242f6280
md"""
# Population growth with harvesting

A constant harvest lowers the logistic equilibrium; this example stops at extinction instead of continuing to negative population.
"""

# ╔═╡ c7e1d611-664d-517d-b9d4-5447cca130c5
md"""
```math
\begin{aligned}
\dot{N} &= rN\left(1-\frac{N}{K}\right)-H.
\end{aligned}
```
"""

# ╔═╡ e148d1f4-92be-5e85-a78a-3896e4825940
function model!(du, u, p, t)
    N = only(u)
    r, K, H = p
    du[1] = r*N*(1-N/K) - H
    nothing
end

# ╔═╡ a409ae6b-93fe-5b85-b2f3-28f6fa8c1aa3
u0 = [4.0]

# ╔═╡ 40136758-4be9-5dae-b930-250f5463c65a
tspan = (0.0, 20.0)

# ╔═╡ 053562d4-fdc6-5ee2-a144-e2ec5e3c6743
p = [1.0, 10.0, 1.5]  # r, K, H

# ╔═╡ 68f2e8f6-dd9d-5f89-a80e-19b8f46a8213
prob = ODEProblem(model!, u0, tspan, p)

# ╔═╡ 5778074b-0427-5c24-9266-728bd6b3c481
extinction = ContinuousCallback((u,t,integrator) -> u[1], nothing, terminate!)

# ╔═╡ 2cbba37a-6d8d-584e-b5f4-2b5b1411f8f0
sol = solve(prob, Tsit5(); abstol=1e-9, reltol=1e-7, callback=extinction, saveat=0.02);

# ╔═╡ eb8a3047-1699-554c-85fa-0d56c56e3595
plot(sol; idxs=[1], xlabel="t", ylabel="state")

# ╔═╡ 2bba9cec-1794-59fb-b0a8-270767292247
md"""
## Harvest threshold and recovery time

For $r,K>0$, constant harvesting has a saddle-node at
$H_c=rK/4$. Below it, the larger equilibrium attracts and the smaller one is
a threshold for extinction:

```math
N_\pm=\frac K2\left(1\pm\sqrt{1-\frac{4H}{rK}}\right),\qquad
\tau_{\rm recovery}=\frac{1}{r\sqrt{1-H/H_c}}.
```

The recovery time is obtained by linearizing at $N_+$, so it describes small
perturbations. Its divergence explains **critical slowing down** before collapse.
"""

# ╔═╡ 0528df99-5afe-5f06-bae4-44615193c8ae
let r=p[1], K=p[2]
    if r==0
        md"With zero growth there is no positive sustainable constant harvest."
    else
        Hc=r*K/4
        hs=range(0,Hc;length=201)
        root=sqrt.(max.(0,1 .-hs./Hc))
        a=plot(hs,K/2 .* (1 .+root);label="stable",xlabel="H",ylabel="N*")
        plot!(a,hs,K/2 .* (1 .-root);label="unstable",linestyle=:dash)
        scatter!(a,[Hc],[K/2];label="saddle-node",color=:black)
        fractions=range(0,0.99;length=201)
        b=plot(fractions,1 ./ (r.*sqrt.(1 .-fractions));
            xlabel="H/Hc",ylabel="recovery time",legend=false)
        plot(a,b;layout=(1,2),size=(900,370),margin=5*Plots.mm)
    end
end

# ╔═╡ a2314024-d703-5103-8e9f-753c932e713b
md"""
## Outbreaks, two potential wells, and hysteresis

Saturating predation can allow both a small and a large population to persist:

```math
\dot N=rN(1-N/K)-P\frac{N^2}{1+N^2},\qquad
V(N)=-\frac{rN^2}{2}+\frac{rN^3}{3K}+P(N-\arctan N).
```

Here $-\partial_N V=\dot N$. Two minima separated by a maximum mean two
attractors separated by an unstable population threshold. This differs from
constant harvesting: $N=0$ is an equilibrium and positive solutions stay nonnegative.
"""

# ╔═╡ 0a6a94ef-e618-5093-9958-e0e3332e067a
function course_outbreak!(du,u,p,t)
    r,K,P=p
    N=u[1]
    du[1]=r*N*(1-N/K)-P*N^2/(1+N^2)
    nothing
end

# ╔═╡ e5696789-11d4-583f-b24b-9b68d6de7bf2
course_outbreak_u0 = [0.1]

# ╔═╡ 32b731c3-03ff-5df4-8d0a-6f61836a0b2f
course_outbreak_tspan = (0.0,300.0)

# ╔═╡ 4406ccd2-b281-5f5d-985b-9fb06c450435
course_outbreak_p = [0.5,8.0,1.0]

# ╔═╡ 1120e447-cea2-5bfc-879e-9e768beac13a
course_outbreak_prob = ODEProblem(course_outbreak!, course_outbreak_u0, course_outbreak_tspan, course_outbreak_p)

# ╔═╡ 20f63bad-231c-59f2-9557-78866223e602
course_outbreak_sol = solve(course_outbreak_prob, Tsit5(); abstol=1e-9, reltol=1e-8);

# ╔═╡ 69e82baa-47f5-5516-9039-cbaa4b62334f
course_outbreak_high = solve(remake(course_outbreak_prob;u0=[7.0]),Tsit5();abstol=1e-9,reltol=1e-8);

# ╔═╡ 915b4906-2f03-5510-87ca-bf15a34f3bbe
let f=plot(course_outbreak_sol;label="small initial population",xlabel="t",ylabel="N")
    plot!(f,course_outbreak_high;label="large initial population")
end

# ╔═╡ dee4ac6d-6d6d-50dd-84c8-19a0b46db4f0
let r=course_outbreak_p[1], K=course_outbreak_p[2], P=course_outbreak_p[3]
    ns=range(0,max(K,1);length=501)
    f=r.*ns.*(1 .-ns./K) .-P.*ns.^2 ./ (1 .+ns.^2)
    V=-r.*ns.^2 ./2 .+r.*ns.^3 ./(3K) .+P.*(ns .-atan.(ns))
    a=plot(ns,f;xlabel="N",ylabel="dN/dt",legend=false)
    hline!(a,[0];color=:gray)
    b=plot(ns,V;xlabel="N",ylabel="V(N)",legend=false)
    plot(a,b;layout=(1,2),size=(900,370),margin=5*Plots.mm)
end

# ╔═╡ 53e3034d-fe39-58ff-b537-ef9def251a10
md"""
For positive equilibria, solving for the growth parameter gives a branch
without root-finding:

```math
r(N)=\frac{PN}{(1+N^2)(1-N/K)},\quad 0<N<K.
```

The next plot uses $K=8$, $P=1$ to show both folds. Stability is evaluated
from $f'(N)=r(1-2N/K)-2PN/(1+N^2)^2$. A slow sweep follows an attracting
branch until it ends, then jumps to the other one; the jump back occurs at a
different fold. Arrows would describe this sweep, not the direction of time along
the equilibrium curve itself.
"""

# ╔═╡ c7d648d5-7285-5061-be5d-ac505ed27573
let K=8.0, P=1.0, ns=range(0.001,7.5;length=1200)
    rates=P.*ns ./((1 .+ns.^2).*(1 .-ns./K))
    slopes=rates.*(1 .-2 .*ns./K) .-2P.*ns./(1 .+ns.^2).^2
    stable=ifelse.(slopes.<0,ns,NaN)
    unstable=ifelse.(slopes.>=0,ns,NaN)
    f=plot(rates,stable;label="stable",xlabel="r",ylabel="N*",xlims=(0.3,0.8))
    plot!(f,rates,unstable;label="unstable",linestyle=:dash)
end

# ╔═╡ f0cfcae9-b14e-5e87-a4e2-128c57287c13
md"""
**Try:** use $r=0.5$, $K=8$, $P=1$ and compare initial populations $0.1$ and $7$.
Locate their destinations in the potential. Then vary $r$ slowly in both
directions, using the previous final population as the next initial condition.
Does the same $r$ always give the same final population?
"""

# ╔═╡ Cell order:
# ╠═4df18ed1-628e-546a-a9f0-44a7912b76d0
# ╟─f1b40b34-ca27-5718-9166-0958242f6280
# ╟─c7e1d611-664d-517d-b9d4-5447cca130c5
# ╠═e148d1f4-92be-5e85-a78a-3896e4825940
# ╠═a409ae6b-93fe-5b85-b2f3-28f6fa8c1aa3
# ╠═40136758-4be9-5dae-b930-250f5463c65a
# ╠═053562d4-fdc6-5ee2-a144-e2ec5e3c6743
# ╠═68f2e8f6-dd9d-5f89-a80e-19b8f46a8213
# ╠═5778074b-0427-5c24-9266-728bd6b3c481
# ╠═2cbba37a-6d8d-584e-b5f4-2b5b1411f8f0
# ╠═eb8a3047-1699-554c-85fa-0d56c56e3595
# ╟─2bba9cec-1794-59fb-b0a8-270767292247
# ╠═0528df99-5afe-5f06-bae4-44615193c8ae
# ╟─a2314024-d703-5103-8e9f-753c932e713b
# ╠═0a6a94ef-e618-5093-9958-e0e3332e067a
# ╠═e5696789-11d4-583f-b24b-9b68d6de7bf2
# ╠═32b731c3-03ff-5df4-8d0a-6f61836a0b2f
# ╠═4406ccd2-b281-5f5d-985b-9fb06c450435
# ╠═1120e447-cea2-5bfc-879e-9e768beac13a
# ╠═20f63bad-231c-59f2-9557-78866223e602
# ╠═69e82baa-47f5-5516-9039-cbaa4b62334f
# ╠═915b4906-2f03-5510-87ca-bf15a34f3bbe
# ╠═dee4ac6d-6d6d-50dd-84c8-19a0b46db4f0
# ╟─53e3034d-fe39-58ff-b537-ef9def251a10
# ╠═c7d648d5-7285-5061-be5d-ac505ed27573
# ╟─f0cfcae9-b14e-5e87-a4e2-128c57287c13
