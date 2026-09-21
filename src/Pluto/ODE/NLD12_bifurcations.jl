### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ ff555515-e183-5212-836c-842b5fb4e043
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using NonLinearDynamics: flow2d_nullclines
    using DifferentialEquations, Plots
    using LinearAlgebra
end

# ╔═╡ 3ea14c7f-45a9-5c8f-b6ef-3990e2b33b65
md"""
# From local bifurcations to codimension two

Follow the lecture progression: first equilibria and small cycles, then
finite-amplitude cycles and global connections, and finally two-parameter
unfoldings. In each example, predict the attractors from the equations,
compare phase and time plots, and only then read a bifurcation diagram.
"""

# ╔═╡ c963fcf0-dcd2-5c30-9690-66f889710723
md"""
## A saddle-node in the plane

A fast attracting direction can accompany an effectively one-dimensional
bifurcation:

```math
\dot x=a-y,\qquad \dot y=x^2-by,\qquad b>0.
```

For $a>0$, the nullclines intersect at $(x_*,y_*)=(\pm\sqrt{ab},a)$.
The Jacobian has trace $-b$ and determinant $2x_*$: the negative root is a
saddle, the positive root an attractor. They collide at $a=0$.
For large $b$, $y$ relaxes quickly near $y=x^2/b$; the slow drift along it
resembles the scalar saddle-node.
"""

# ╔═╡ 6250fe4a-012e-565b-8a9b-92a8e4d35eea
function course_saddle_node!(du,u,p,t)
    a,b=p
    x,y=u
    du[1]=a-y
    du[2]=x^2-b*y
    nothing
end

# ╔═╡ 2a59de4c-0b96-5f29-864f-afa618b324bd
course_sn_a = 0.05  # a (b = 2)

# ╔═╡ 1dbae155-1808-5f38-aba4-dc7ee0753b34
course_sn_u0 = [0.1,0.1]

# ╔═╡ f8a08c17-78c3-50f4-b661-7d946145892e
course_sn_tspan = (0.0,60.0)

# ╔═╡ 9c773740-8846-5051-ba58-2ce07dae318c
course_sn_p = [course_sn_a,2.0]

# ╔═╡ 384fa86c-319f-53c0-89a7-307cfc940fb4
course_sn_prob = ODEProblem(course_saddle_node!, course_sn_u0, course_sn_tspan, course_sn_p)

# ╔═╡ 8da01a4b-58d8-5ce0-917b-2de62c122d11
course_sn_sol = solve(course_sn_prob, Tsit5(); abstol=1e-9, reltol=1e-8, callback=ContinuousCallback((u,t,i)->6-maximum(abs,u),nothing,terminate!));

# ╔═╡ 3753cd45-778d-5846-ae5b-40a420846bb7
let f=flow2d_nullclines(course_saddle_node!,course_sn_p;xlims=[-1.0,1.0],
        ylims=[-0.3,1.0],regions=false,vectorfield=true,title="Planar saddle-node")
    plot!(f,course_sn_sol;idxs=(1,2),color=:black)
    if course_sn_a>=0
        scatter!(f,[-sqrt(2course_sn_a),sqrt(2course_sn_a)],[course_sn_a,course_sn_a];color=:black)
    end
    plot!(f;xlims=(-1,1),ylims=(-0.3,1))
end

# ╔═╡ 1de73558-0b7a-54f0-97af-faed1a511a90
md"""
## A local oscillatory onset: supercritical Hopf

The normal form separates amplitude and phase:

```math
\dot x=(\beta-r^2)x-\omega y,\qquad
\dot y=\omega x+(\beta-r^2)y,\qquad r^2=x^2+y^2,
```
```math
\dot r=\beta r-r^3,\qquad \dot\theta=\omega.
```

The origin has eigenvalues $\beta\pm i\omega$. For $\omega>0$, it loses
stability at $\beta=0$ and an attracting cycle with radius $\sqrt\beta$
appears for $\beta>0$, with period $2\pi/\omega$. At the bifurcation itself,
$\dot r=-r^3$ still attracts algebraically: purely imaginary eigenvalues do
not make the nonlinear equilibrium a center.
"""

# ╔═╡ a9618c4f-51d5-5264-b81c-ef0df1eb22d2
function course_hopf!(du,u,p,t)
    β,ω=p
    x,y=u
    r2=x^2+y^2
    du[1]=(β-r2)*x-ω*y
    du[2]=ω*x+(β-r2)*y
    nothing
end

# ╔═╡ e87eeab1-bdf5-5a12-8bb4-7319c9d8dabd
course_hopf_beta = 0.2  # β (ω = 1)

# ╔═╡ 08ec05dc-fc9a-552d-8e10-637b291687dd
course_hopf_u0 = [0.1,0.0]

# ╔═╡ 95721dd1-6eb5-5bf2-98cc-b441f8456d6a
course_hopf_tspan = (0.0,60.0)

# ╔═╡ 4e555c95-22ca-5caf-9d2a-7b65f284d967
course_hopf_p = [course_hopf_beta,1.0]

# ╔═╡ 53d62014-3e57-5e10-a25f-8f1c7506d19c
course_hopf_prob = ODEProblem(course_hopf!, course_hopf_u0, course_hopf_tspan, course_hopf_p)

# ╔═╡ 784f77f1-c42a-5e2a-8454-9dcbd85c3cd0
course_hopf_sol = solve(course_hopf_prob, Tsit5(); abstol=1e-9, reltol=1e-8);

# ╔═╡ ad60ddbb-da7e-57b3-98be-b8f210d028e4
let phase=plot(course_hopf_sol;idxs=(1,2),xlabel="x",ylabel="y",legend=false,aspect_ratio=1),
    amplitude=plot(course_hopf_sol.t,hypot.(course_hopf_sol[1,:],course_hopf_sol[2,:]);
        xlabel="t",ylabel="r(t)",label="trajectory")
    hline!(amplitude,[sqrt(max(0,course_hopf_beta))];label="limiting radius",linestyle=:dash)
    plot(phase,amplitude;layout=(1,2),size=(900,390),margin=5*Plots.mm)
end

# ╔═╡ 1405f5c9-d4df-5051-b072-8d7669afe8e0
md"""
### Hopf in a population model

Replace unlimited prey growth by logistic growth and let predation saturate.
If predator carrying capacity is proportional to prey density, the
dimensionless equations become

```math
\dot x=x(1-x)-\frac{axy}{x+c},\qquad
\dot y=by(1-y/x),\qquad x>0,\ a,b,c>0.
```

Here $a$ measures predation, $b$ is the predator's relative growth rate, and
$c$ sets the saturation scale. The positive equilibrium is
$x_*=y_*=(1-a-c+\sqrt{(1-a-c)^2+4c})/2$. At it,

```math
J_*=\begin{pmatrix}
1-2x_*-\dfrac{acx_*}{(x_*+c)^2}&-\dfrac{ax_*}{x_*+c}\\
b&-b
\end{pmatrix}.
```

With $a=1$, $c=0.1$, varying $b$ changes stability without moving the
nullclines or their intersection. Compare an attracting cycle here with
the family of center orbits in Lotka–Volterra.
"""

# ╔═╡ 74367db9-24e6-5bfc-ba41-95c3421ae7bb
function course_predator!(du,u,p,t)
    a,b,c=p
    x,y=u
    du[1]=x*(1-x)-a*x*y/(x+c)
    du[2]=b*y*(1-y/x)
    nothing
end

# ╔═╡ 8c5a9ecf-2abe-5cae-92c9-9465e7465c2b
course_predator_b = 0.2  # b (a = 1, c = 0.1)

# ╔═╡ cb931bdf-b3db-50e6-b98e-c1fb9ebb1a2e
course_predator_u0 = [0.3,0.25]

# ╔═╡ fcfbe757-5158-5870-9c50-7a853325dae3
course_predator_tspan = (0.0,500.0)

# ╔═╡ 17f99cf1-f837-5b9c-8bad-39df86064c35
course_predator_p = [1.0,course_predator_b,0.1]

# ╔═╡ 4aa91ad9-c8e9-58fa-b375-147b41c918aa
course_predator_prob = ODEProblem(course_predator!, course_predator_u0, course_predator_tspan, course_predator_p)

# ╔═╡ fff25d81-1c59-5d0f-a6c5-50ba024046b2
course_predator_sol = solve(course_predator_prob, Tsit5(); abstol=1e-9, reltol=1e-8, isoutofdomain=(u,p,t)->any(x->x<=0,u));

# ╔═╡ d5e77cd2-9870-50ab-bbb8-fabd756a95b6
let xs=range(0.005,1;length=401)
    a,b,c=course_predator_p
    phase=plot(course_predator_sol;idxs=(1,2),label="trajectory",xlabel="prey x",ylabel="predator y")
    plot!(phase,xs,(xs .+c).*(1 .-xs)./a;label="dx/dt = 0",color=:red)
    plot!(phase,xs,xs;label="dy/dt = 0",color=:green)
    trace=plot(course_predator_sol;label=["x" "y"],xlabel="t")
    plot(phase,trace;layout=(1,2),size=(950,390),margin=5*Plots.mm)
end

# ╔═╡ 0b13480a-b93d-5895-a209-099eafdae9f7
let
    a,b,c=course_predator_p
    xstar=(1-a-c+sqrt((1-a-c)^2+4c))/2
    J=[1-2xstar-a*c*xstar/(xstar+c)^2 -a*xstar/(xstar+c); b -b]
    (equilibrium=[xstar,xstar],trace=tr(J),determinant=det(J),
        eigenvalues=eigvals(J),b_at_zero_trace=1-2xstar-a*c*xstar/(xstar+c)^2)
end

# ╔═╡ e9556def-31c4-534e-8ada-8d877d5f4893
md"""
## Finite-amplitude onset: a fold of cycles

An additional radial term allows a stable equilibrium and an attracting
cycle to coexist:

```math
\dot r=r(\beta+r^2-r^4),\quad \dot\theta=\omega,
\qquad
\dot x=(\beta+r^2-r^4)x-\omega y,\quad
\dot y=\omega x+(\beta+r^2-r^4)y.
```

For $-1/4<\beta<0$, the origin attracts, the inner cycle repels, and the
outer cycle attracts. Their radii are
$r_\pm=\sqrt{(1\pm\sqrt{1+4\beta})/2}$. The two cycles collide at
$\beta=-1/4$; the inner cycle shrinks into the origin in a **subcritical
Hopf** at $\beta=0$. Initial conditions on different sides of the inner
cycle have different destinations, producing hysteresis under slow sweeps.
"""

# ╔═╡ cbf72d4f-2fed-5945-817f-a9e0a72946b7
function course_fold_cycles!(du,u,p,t)
    β,ω=p
    x,y=u
    r2=x^2+y^2
    growth=β+r2-r2^2
    du[1]=growth*x-ω*y
    du[2]=ω*x+growth*y
    nothing
end

# ╔═╡ c0939af7-1d87-53c2-a4a7-a96bc73b483a
course_fold_beta = -0.15  # β (ω = 1)

# ╔═╡ 5ef5fc81-0405-5892-b4b5-6c641efcdf6f
course_fold_u0 = [0.2,0.0]

# ╔═╡ 08bf666a-0988-5791-9a1d-893b366b1a02
course_fold_tspan = (0.0,100.0)

# ╔═╡ b31c4088-70b5-56d1-aefa-c4c4a880b8a7
course_fold_p = [course_fold_beta,1.0]

# ╔═╡ bf73880e-045d-5152-a514-6f2e27249c8a
course_fold_prob = ODEProblem(course_fold_cycles!, course_fold_u0, course_fold_tspan, course_fold_p)

# ╔═╡ cbc009d4-d297-5662-bf3f-4df66c383819
course_fold_sol = solve(course_fold_prob, Tsit5(); abstol=1e-9, reltol=1e-8);

# ╔═╡ 474897c6-6e43-5920-9c23-2a2f3b76223d
course_fold_outer = solve(remake(course_fold_prob;u0=[1.2,0.0]),Tsit5();abstol=1e-9,reltol=1e-8);

# ╔═╡ b3702979-fc55-5f09-8684-d67b9e48f85f
let phase=plot(course_fold_sol;idxs=(1,2),label="small initial radius",xlabel="x",ylabel="y",aspect_ratio=1)
    plot!(phase,course_fold_outer;idxs=(1,2),label="large initial radius")
    if -0.25<course_fold_beta<0
        θ=range(0,2pi;length=401)
        inner=sqrt((1-sqrt(1+4course_fold_beta))/2)
        plot!(phase,inner.*cos.(θ),inner.*sin.(θ);label="unstable cycle",linestyle=:dash,color=:black)
    end
    beta=range(-0.25,0.1;length=401)
    branch=plot(beta,sqrt.((1 .+sqrt.(1 .+4beta))./2);label="stable outer cycle",
        xlabel="β",ylabel="radius")
    negative=range(-0.25,0;length=301)
    plot!(branch,negative,sqrt.((1 .-sqrt.(1 .+4negative))./2);
        label="unstable inner cycle",linestyle=:dash)
    plot!(branch,[-0.3,0],[0,0];label="stable origin",color=:steelblue)
    plot!(branch,[0,0.1],[0,0];label="unstable origin",color=:darkorange,linestyle=:dash)
    vline!(branch,[course_fold_beta];label="selected β",color=:gray)
    plot(phase,branch;layout=(1,2),size=(1050,410),margin=5*Plots.mm)
end

# ╔═╡ ac1ea4d6-a2c6-58e2-b8ca-3a083faee849
md"""
### Distinguishing oscillatory onsets

| Mechanism | Amplitude near onset | Period near onset |
| --- | --- | --- |
| Supercritical Hopf | Shrinks like $\sqrt\delta$ | Finite |
| Fold of cycles | Finite | Finite |
| SNIC / SNIPER (NLD03b) | Finite on the circle | Diverges like $1/\sqrt\delta$ |
| Saddle homoclinic loop | Finite | Diverges like $-\log|\delta|$ |

Here $\delta$ is distance to the relevant bifurcation on the oscillating side;
the scalings assume generic cases. Near a homoclinic loop a trajectory spends
a long time near the saddle. Local equilibrium eigenvalues alone cannot locate
that global connection. The Bogdanov–Takens continuation below brings
saddle-node, Hopf, and homoclinic curves into the same parameter plane.
"""

# ╔═╡ 09a80719-4b5a-53b3-9d92-c5763b54ad79
md"""
## Two parameters: the cusp unfolding

Breaking pitchfork symmetry with an offset gives the scalar family

```math
\dot x=a+bx-x^3.
```

Fold points satisfy both $f=0$ and $f_x=0$, giving
$(a,b)=(-2x^3,3x^2)$. They meet at the cusp $(0,0)$.
There are three distinct equilibria inside $27a^2<4b^3$ ($b>0$), and one
outside. A **two-parameter plot** does not by itself make a bifurcation
codimension two: generic points on either fold curve are still codimension one.
"""

# ╔═╡ e2208589-f8bf-511b-b93a-5e85585a7bc5
function course_cusp!(du,u,p,t)
    a,b=p
    x=u[1]
    du[1]=a+b*x-x^3
    nothing
end

# ╔═╡ c96cb049-9781-5e09-a1cd-43d9ff3b63fc
course_cusp_a = 0.05  # a (b = 0.5)

# ╔═╡ 1cd54891-6057-5207-8e09-4f66a9198bf9
course_cusp_u0 = [-1.0]

# ╔═╡ 22a545d2-41c7-5cee-a79a-483d17cb356c
course_cusp_tspan = (0.0,40.0)

# ╔═╡ c34b5399-cd8a-5fe1-8754-2787684bc604
course_cusp_p = [course_cusp_a,0.5]

# ╔═╡ 50197c42-fdfb-5999-bf9b-942fd613542f
course_cusp_prob = ODEProblem(course_cusp!, course_cusp_u0, course_cusp_tspan, course_cusp_p)

# ╔═╡ 4ac6be46-5099-5d2f-a1ea-4a149f708e59
course_cusp_sol = solve(course_cusp_prob, Tsit5(); abstol=1e-9, reltol=1e-8);

# ╔═╡ bc904ac9-d3a9-5802-85bc-5a2b025dd542
let xs=range(-0.7,0.7;length=501), roots=range(-1.2,1.2;length=701), b=course_cusp_p[2]
    plane=plot(-2 .*xs.^3,3 .*xs.^2;xlabel="a",ylabel="b",label="fold curves")
    scatter!(plane,[course_cusp_a],[b];label="selected parameters")
    offset=roots.^3 .-b.*roots
    slope=b .-3 .*roots.^2
    branch=plot(offset,ifelse.(slope.<0,roots,NaN);label="stable",xlabel="a",ylabel="x*",xlims=(-0.4,0.4))
    plot!(branch,offset,ifelse.(slope.>=0,roots,NaN);label="unstable",linestyle=:dash)
    scatter!(branch,[course_cusp_a],[course_cusp_sol.u[end][1]];label="ODE final state")
    plot(plane,branch;layout=(1,2),size=(950,390),margin=5*Plots.mm)
end

# ╔═╡ 0d3358e4-6e2a-59a4-a0b4-b69b8500e9cc
md"""
## Bogdanov–Takens: two zero eigenvalues

The cusp concerns a scalar degeneracy. Bogdanov–Takens instead requires a
planar equilibrium with a **double zero eigenvalue and one eigenvector**.
Its unfolding organizes nearby saddle-node, Hopf, and homoclinic
bifurcations. It is a bifurcation, not a physical model.

For the cubic extension used below, the origin at $\mu_1=\mu_2=0$ has
$J=\left(\begin{smallmatrix}0&1\\0&0\end{smallmatrix}\right)$.
Use the following example to connect the local pictures with continuation
of equilibria and cycles. Time integration finds attracting behavior;
continuation can also follow unstable branches.
"""

# ╔═╡ 4b6c1283-065c-5492-8815-05773296f288
md"""
# Bifurcations: the Bogdanov–Takens unfolding

This two-parameter unfolding of the codimension-two Bogdanov–Takens bifurcation includes cubic terms, as in interactive NLD12.
"""

# ╔═╡ 72b97f9b-94e2-574a-9870-5a8c2d808db6
md"""
```math
\begin{aligned}
\dot{x} &= v,\\
\dot{v} &= \mu_1+\mu_2x+x^2-xv-x^3-x^2v.
\end{aligned}
```
"""

# ╔═╡ da860a78-a2d3-5d01-b0e2-3611026f6c7f
function model!(du, u, p, t)
    x, v = u
    μ1, μ2 = p
    du[1] = v
    du[2] = μ1 + x*(μ2-v+x*(1-x-v))
    nothing
end

# ╔═╡ 1fbfcdab-37de-5762-97b8-b287998ca1eb
u0 = [0.1, 0.0]

# ╔═╡ 87f60eb2-2360-5c27-8256-3831e3f095eb
tspan = (0.0, 1000.0)

# ╔═╡ 4a92bfa0-2d35-50b3-999f-d0a08d838804
p = [-0.02, -0.16]  # μ1, μ2

# ╔═╡ 930b0512-c8d5-57ce-9091-8c9ade3935dc
prob = ODEProblem(model!, u0, tspan, p)

# ╔═╡ 54c1cff6-bd44-50a5-a6ce-bca3b4e62e0e
sol = solve(prob, Tsit5(); abstol=1e-9, reltol=1e-7, saveat=0.02);

# ╔═╡ 023ada46-bf58-5ad0-b33a-fc6a8401a79d
plot(sol; idxs=[1, 2], xlabel="t", ylabel="state")

# ╔═╡ 2ad41245-d960-5c10-9ff3-cf349f6dbc88
plot(sol; idxs=(1, 2), legend=false)

# ╔═╡ Cell order:
# ╠═ff555515-e183-5212-836c-842b5fb4e043
# ╟─3ea14c7f-45a9-5c8f-b6ef-3990e2b33b65
# ╟─c963fcf0-dcd2-5c30-9690-66f889710723
# ╠═6250fe4a-012e-565b-8a9b-92a8e4d35eea
# ╠═2a59de4c-0b96-5f29-864f-afa618b324bd
# ╠═1dbae155-1808-5f38-aba4-dc7ee0753b34
# ╠═f8a08c17-78c3-50f4-b661-7d946145892e
# ╠═9c773740-8846-5051-ba58-2ce07dae318c
# ╠═384fa86c-319f-53c0-89a7-307cfc940fb4
# ╠═8da01a4b-58d8-5ce0-917b-2de62c122d11
# ╠═3753cd45-778d-5846-ae5b-40a420846bb7
# ╟─1de73558-0b7a-54f0-97af-faed1a511a90
# ╠═a9618c4f-51d5-5264-b81c-ef0df1eb22d2
# ╠═e87eeab1-bdf5-5a12-8bb4-7319c9d8dabd
# ╠═08ec05dc-fc9a-552d-8e10-637b291687dd
# ╠═95721dd1-6eb5-5bf2-98cc-b441f8456d6a
# ╠═4e555c95-22ca-5caf-9d2a-7b65f284d967
# ╠═53d62014-3e57-5e10-a25f-8f1c7506d19c
# ╠═784f77f1-c42a-5e2a-8454-9dcbd85c3cd0
# ╠═ad60ddbb-da7e-57b3-98be-b8f210d028e4
# ╟─1405f5c9-d4df-5051-b072-8d7669afe8e0
# ╠═74367db9-24e6-5bfc-ba41-95c3421ae7bb
# ╠═8c5a9ecf-2abe-5cae-92c9-9465e7465c2b
# ╠═cb931bdf-b3db-50e6-b98e-c1fb9ebb1a2e
# ╠═fcfbe757-5158-5870-9c50-7a853325dae3
# ╠═17f99cf1-f837-5b9c-8bad-39df86064c35
# ╠═4aa91ad9-c8e9-58fa-b375-147b41c918aa
# ╠═fff25d81-1c59-5d0f-a6c5-50ba024046b2
# ╠═d5e77cd2-9870-50ab-bbb8-fabd756a95b6
# ╠═0b13480a-b93d-5895-a209-099eafdae9f7
# ╟─e9556def-31c4-534e-8ada-8d877d5f4893
# ╠═cbf72d4f-2fed-5945-817f-a9e0a72946b7
# ╠═c0939af7-1d87-53c2-a4a7-a96bc73b483a
# ╠═5ef5fc81-0405-5892-b4b5-6c641efcdf6f
# ╠═08bf666a-0988-5791-9a1d-893b366b1a02
# ╠═b31c4088-70b5-56d1-aefa-c4c4a880b8a7
# ╠═bf73880e-045d-5152-a514-6f2e27249c8a
# ╠═cbc009d4-d297-5662-bf3f-4df66c383819
# ╠═474897c6-6e43-5920-9c23-2a2f3b76223d
# ╠═b3702979-fc55-5f09-8684-d67b9e48f85f
# ╟─ac1ea4d6-a2c6-58e2-b8ca-3a083faee849
# ╟─09a80719-4b5a-53b3-9d92-c5763b54ad79
# ╠═e2208589-f8bf-511b-b93a-5e85585a7bc5
# ╠═c96cb049-9781-5e09-a1cd-43d9ff3b63fc
# ╠═1cd54891-6057-5207-8e09-4f66a9198bf9
# ╠═22a545d2-41c7-5cee-a79a-483d17cb356c
# ╠═c34b5399-cd8a-5fe1-8754-2787684bc604
# ╠═50197c42-fdfb-5999-bf9b-942fd613542f
# ╠═4ac6be46-5099-5d2f-a1ea-4a149f708e59
# ╠═bc904ac9-d3a9-5802-85bc-5a2b025dd542
# ╟─0d3358e4-6e2a-59a4-a0b4-b69b8500e9cc
# ╟─4b6c1283-065c-5492-8815-05773296f288
# ╟─72b97f9b-94e2-574a-9870-5a8c2d808db6
# ╠═da860a78-a2d3-5d01-b0e2-3611026f6c7f
# ╠═1fbfcdab-37de-5762-97b8-b287998ca1eb
# ╠═87f60eb2-2360-5c27-8256-3831e3f095eb
# ╠═4a92bfa0-2d35-50b3-999f-d0a08d838804
# ╠═930b0512-c8d5-57ce-9091-8c9ade3935dc
# ╠═54c1cff6-bd44-50a5-a6ce-bca3b4e62e0e
# ╠═023ada46-bf58-5ad0-b33a-fc6a8401a79d
# ╠═2ad41245-d960-5c10-9ff3-cf349f6dbc88
