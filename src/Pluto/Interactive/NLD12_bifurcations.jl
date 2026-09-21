### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 5e1d221c-18f7-11ee-20bc-b5624af0581e
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using NonLinearDynamics: flow2d_nullclines
    using Plots, DifferentialEquations, Setfield, ForwardDiff, PlutoUI, IntervalRootFinding, StaticArrays
    import BifurcationKit as BK
    import HclinicBifurcationKit as HBK
    using LinearAlgebra
end

# ╔═╡ 3ef551ed-71c8-4c56-8a14-46e70adbae04
PlutoUI.TableOfContents()

# ╔═╡ 9b3c7bab-5036-5f1b-9e44-73254ea506b0
md"""
# From local bifurcations to codimension two

Follow the lecture progression: first equilibria and small cycles, then
finite-amplitude cycles and global connections, and finally two-parameter
unfoldings. In each example, predict the attractors from the equations,
compare phase and time plots, and only then read a bifurcation diagram.
"""

# ╔═╡ 385b660a-05e8-56ab-a144-61c90eac704a
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

# ╔═╡ c285bb10-108c-5a51-8ca4-27756600b444
function course_saddle_node!(du,u,p,t)
    a,b=p
    x,y=u
    du[1]=a-y
    du[2]=x^2-b*y
    nothing
end

# ╔═╡ 26ee2eab-8184-5def-bc9a-3c254d018f82
md"""
a (b = 2): $(@bind course_sn_a Slider(-0.1:0.01:0.1; default=0.05, show_value=true))
"""

# ╔═╡ 671a891a-34bf-5470-9eec-686511a63990
course_sn_u0 = [0.1,0.1]

# ╔═╡ 77f8aa93-cfba-5b38-83f0-b1822d352a52
course_sn_tspan = (0.0,60.0)

# ╔═╡ 290a3ea8-96a6-5584-985a-b2195a32f96e
course_sn_p = [course_sn_a,2.0]

# ╔═╡ 4c68fded-6f1a-51ca-83e7-45873f5f2208
course_sn_prob = ODEProblem(course_saddle_node!, course_sn_u0, course_sn_tspan, course_sn_p)

# ╔═╡ 5f962302-4dee-5326-a383-e3af9e5e8344
course_sn_sol = solve(course_sn_prob, Tsit5(); abstol=1e-9, reltol=1e-8, callback=ContinuousCallback((u,t,i)->6-maximum(abs,u),nothing,terminate!));

# ╔═╡ 399a4606-410e-5f3d-ba07-3d1b05d9cb1e
let f=flow2d_nullclines(course_saddle_node!,course_sn_p;xlims=[-1.0,1.0],
        ylims=[-0.3,1.0],regions=false,vectorfield=true,title="Planar saddle-node")
    plot!(f,course_sn_sol;idxs=(1,2),color=:black)
    if course_sn_a>=0
        scatter!(f,[-sqrt(2course_sn_a),sqrt(2course_sn_a)],[course_sn_a,course_sn_a];color=:black)
    end
    plot!(f;xlims=(-1,1),ylims=(-0.3,1))
end

# ╔═╡ eb629797-35e8-5b4d-870a-389e02eb5a83
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

# ╔═╡ f723d036-5c8b-5c58-9870-25195741cca9
function course_hopf!(du,u,p,t)
    β,ω=p
    x,y=u
    r2=x^2+y^2
    du[1]=(β-r2)*x-ω*y
    du[2]=ω*x+(β-r2)*y
    nothing
end

# ╔═╡ cb10ff3c-130a-50a6-beb7-10a201a8bd09
md"""
β (ω = 1): $(@bind course_hopf_beta Slider(-0.4:0.02:0.4; default=0.2, show_value=true))
"""

# ╔═╡ 328597eb-805f-59f4-9c6d-abc023133e94
course_hopf_u0 = [0.1,0.0]

# ╔═╡ 6596bdf7-e667-5176-9281-08e923e4cef7
course_hopf_tspan = (0.0,60.0)

# ╔═╡ cf7b5eaf-43a9-5474-985e-96f157055c4f
course_hopf_p = [course_hopf_beta,1.0]

# ╔═╡ 7406136b-b6e4-5e76-9a68-fa37dbf9bde2
course_hopf_prob = ODEProblem(course_hopf!, course_hopf_u0, course_hopf_tspan, course_hopf_p)

# ╔═╡ 1421e8c7-c81e-54f3-a52d-547881d2f819
course_hopf_sol = solve(course_hopf_prob, Tsit5(); abstol=1e-9, reltol=1e-8);

# ╔═╡ d2cd519a-da55-54f2-a2c5-7bf529c98a68
let phase=plot(course_hopf_sol;idxs=(1,2),xlabel="x",ylabel="y",legend=false,aspect_ratio=1),
    amplitude=plot(course_hopf_sol.t,hypot.(course_hopf_sol[1,:],course_hopf_sol[2,:]);
        xlabel="t",ylabel="r(t)",label="trajectory")
    hline!(amplitude,[sqrt(max(0,course_hopf_beta))];label="limiting radius",linestyle=:dash)
    plot(phase,amplitude;layout=(1,2),size=(900,390),margin=5*Plots.mm)
end

# ╔═╡ b206a81e-3460-5737-b159-03f6a9d2c0e6
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

# ╔═╡ cd2f71d3-d79b-5241-a1dd-33b8dc0d86ea
function course_predator!(du,u,p,t)
    a,b,c=p
    x,y=u
    du[1]=x*(1-x)-a*x*y/(x+c)
    du[2]=b*y*(1-y/x)
    nothing
end

# ╔═╡ ea97d7be-4714-5c60-a3df-ceb787e2a480
md"""
b (a = 1, c = 0.1): $(@bind course_predator_b Slider(0.1:0.01:0.5; default=0.2, show_value=true))
"""

# ╔═╡ 2ff23db5-491d-508e-8c30-5594a2298f1b
course_predator_u0 = [0.3,0.25]

# ╔═╡ 99c9a459-13d2-5c4a-898b-409ccc88e235
course_predator_tspan = (0.0,500.0)

# ╔═╡ 3b2de99d-f775-51ff-8d5b-c22d775e9036
course_predator_p = [1.0,course_predator_b,0.1]

# ╔═╡ bc965a43-02ca-5c83-a79e-e2bfe8f684b4
course_predator_prob = ODEProblem(course_predator!, course_predator_u0, course_predator_tspan, course_predator_p)

# ╔═╡ 76730df7-e671-534b-a62f-51376d9a3c9b
course_predator_sol = solve(course_predator_prob, Tsit5(); abstol=1e-9, reltol=1e-8, isoutofdomain=(u,p,t)->any(x->x<=0,u));

# ╔═╡ 42cc2966-ad31-52cb-9e80-69def66654e0
let xs=range(0.005,1;length=401)
    a,b,c=course_predator_p
    phase=plot(course_predator_sol;idxs=(1,2),label="trajectory",xlabel="prey x",ylabel="predator y")
    plot!(phase,xs,(xs .+c).*(1 .-xs)./a;label="dx/dt = 0",color=:red)
    plot!(phase,xs,xs;label="dy/dt = 0",color=:green)
    trace=plot(course_predator_sol;label=["x" "y"],xlabel="t")
    plot(phase,trace;layout=(1,2),size=(950,390),margin=5*Plots.mm)
end

# ╔═╡ 00847c09-a1ad-5a87-bb7c-93fcec283554
let
    a,b,c=course_predator_p
    xstar=(1-a-c+sqrt((1-a-c)^2+4c))/2
    J=[1-2xstar-a*c*xstar/(xstar+c)^2 -a*xstar/(xstar+c); b -b]
    (equilibrium=[xstar,xstar],trace=tr(J),determinant=det(J),
        eigenvalues=eigvals(J),b_at_zero_trace=1-2xstar-a*c*xstar/(xstar+c)^2)
end

# ╔═╡ 55f5802d-3f2a-5ade-8a59-dd3534240ded
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

# ╔═╡ f885398f-6fae-50f8-aace-b18898efe4d8
function course_fold_cycles!(du,u,p,t)
    β,ω=p
    x,y=u
    r2=x^2+y^2
    growth=β+r2-r2^2
    du[1]=growth*x-ω*y
    du[2]=ω*x+growth*y
    nothing
end

# ╔═╡ dba87d4c-5ba4-5c23-a31b-e23c98ef3521
md"""
β (ω = 1): $(@bind course_fold_beta Slider(-0.3:0.01:0.1; default=-0.15, show_value=true))
"""

# ╔═╡ 8dcfcf37-5f29-5769-aeee-1459b67730b6
course_fold_u0 = [0.2,0.0]

# ╔═╡ 2343f559-74e5-56dc-94cf-1899ea356c75
course_fold_tspan = (0.0,100.0)

# ╔═╡ bb137189-06c4-5d68-8172-64a068f631a7
course_fold_p = [course_fold_beta,1.0]

# ╔═╡ a01be38b-e781-5408-a2c9-1a24ac0e5503
course_fold_prob = ODEProblem(course_fold_cycles!, course_fold_u0, course_fold_tspan, course_fold_p)

# ╔═╡ b0128688-d827-5dd6-9b66-a869bd9739a0
course_fold_sol = solve(course_fold_prob, Tsit5(); abstol=1e-9, reltol=1e-8);

# ╔═╡ 602a7743-3002-5a7d-b4ae-d385a955ea7b
course_fold_outer = solve(remake(course_fold_prob;u0=[1.2,0.0]),Tsit5();abstol=1e-9,reltol=1e-8);

# ╔═╡ 0afbb419-a46f-516b-9de9-9c82b76eda44
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

# ╔═╡ 93e37887-478d-5e70-87f8-475519caa4cc
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

# ╔═╡ a9f3eb77-a4f6-5d1d-91e4-d82d0dbaf6d2
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

# ╔═╡ c17a6aa4-91a4-5e15-bec7-c14f5f5e599c
function course_cusp!(du,u,p,t)
    a,b=p
    x=u[1]
    du[1]=a+b*x-x^3
    nothing
end

# ╔═╡ 14401898-c1f6-5f2e-84ed-81ca7aa77256
md"""
a (b = 0.5): $(@bind course_cusp_a Slider(-0.3:0.01:0.3; default=0.05, show_value=true))
"""

# ╔═╡ 6d8fed6b-7f23-5b71-8e60-ed406cff9ad2
course_cusp_u0 = [-1.0]

# ╔═╡ f2d5e1ad-0f0f-51ac-bd8b-5e6cf413ae28
course_cusp_tspan = (0.0,40.0)

# ╔═╡ 01058bac-b766-50dc-b583-e2a28d762c50
course_cusp_p = [course_cusp_a,0.5]

# ╔═╡ 7fe71b73-6632-5c3c-a36f-f0e04e9812f9
course_cusp_prob = ODEProblem(course_cusp!, course_cusp_u0, course_cusp_tspan, course_cusp_p)

# ╔═╡ 9aba9f01-e96e-558a-bfd2-92d4ffc06d6c
course_cusp_sol = solve(course_cusp_prob, Tsit5(); abstol=1e-9, reltol=1e-8);

# ╔═╡ a5d0bbb6-ff06-5087-b733-612965cd0cec
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

# ╔═╡ a6a935c2-cb82-53f6-b13e-b76230860839
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

# ╔═╡ 77fd29b5-4cc6-4461-ab16-3f0720643dee
md"""
## Bogdanov–Takens unfolding with cubic terms

Bogdanov–Takens is a codimension-two bifurcation: at the critical equilibrium, the Jacobian has a double zero eigenvalue with one eigenvector. Here we study its two-parameter unfolding with additional cubic terms:

$\dot{x} = y$

$\dot{y} = \mu_1+\mu_2x+ x^2 -xy - x^3 -x^2y$ 

In this case by having cubic terms we will have in general one or three fixed points, as in the case of the cusp the fixed points go from 1 to 3 through saddle-node bifurcations that occur in pairs of distinct points. In fact the cubic terms introduce a cusp in addition to the Bogdanov-Takens.

Let us see how the varieties are organized from these new terms.

"""

# ╔═╡ 4380968c-9421-4c34-82bd-850c3027b676
function takens3!(du,u,p,t)
    du[1]=u[2]
    du[2]=p[1]+u[1]*(p[2]-u[2]+u[1]*(1-u[1]-u[2]))
    du
end    

# ╔═╡ e65f74c0-edf0-4378-ab3d-5e848bda5d2f
md"""
# Bifurcation Analysis

## Single parameter Hopf
"""

# ╔═╡ b4574225-ff15-44f5-82e2-5eb779a16c49
begin
	p = (μ1=0.1,μ2=-0.2)
	μ1 = @lens _.μ1
	μ2 = @lens _.μ2
	takens3(u,p) = takens3!(similar(u),u,p,0) # out of place method
	u0 = [0.9,-0.1] #initial condition
end	

# ╔═╡ f3ad61c9-33de-45b0-bee5-89e5f86110e7
opts_br = BK.ContinuationPar(p_min=-0.5,p_max=0.1, ds = -0.001, dsmax = 0.02, detect_bifurcation=3, n_inversion=8);

# ╔═╡ 6ed8af70-4918-4dc2-b419-d87c8f35cfb9
md"""
μ2 $(@bind mu_2 Slider(-0.4:0.005:0.4,default=-0.15;show_value=true)) 
"""

# ╔═╡ 90461200-bc32-4f01-ac06-a7e9b81966df
begin
	# define the bifurcation problem
	prob1 = BK.BifurcationProblem(takens3, u0, set(p,μ2,mu_2), μ1, record_from_solution = (x, p) -> (x = x[1], y = x[2]))
	# continuation options
	# continuation of equilibria
	br1 = BK.continuation(prob1, BK.PALC(tangent=BK.Bordered()),opts_br)
	scene = plot(br1,xlabel="\\mu_1",title=string("BT Cubica \\mu_2 = ",mu_2));
end

# ╔═╡ 9eff4af0-8cef-4a9a-9e6c-092a1e046b96
md"""
1. Definition of the Bifurcation Problem, giving the out of place function, initial condition, parameters (wher we fixed μ2 value with set), and then the parameter to be explored (μ1). With `record_from_solution = (x, p) -> (x = x[1], y = x[2])` we are recording only the variables x, y

2. Options for the continuation algorithm: starting from the higher value of μ1 (0.1) with a stable fixed point and going backwards (ds<0) down to μ1=-0.5. `detect_bifurcation=3` improves accuracy and `n_inversion=8` the convergence.

3. Continuation algorithm using Pseudo Arc-Length Continuation with a Bordered predictor. It is based in computing the tangent to the curve of solutions in an extended plane (x,p) by solving a bordered linear system $F_x dx + F_p dp = 0$ ; $\theta dx_0 dx + (1-\theta) dp_0 dp = 1$
"""

# ╔═╡ fa557b12-9017-4166-ac5b-ead0c3e108b0
for n = 1:(length(br1.specialpoint)-1)
	print(BK.get_normal_form(br1, n))
	println()
end	

# ╔═╡ 1b94418e-60b7-4da4-8720-4ac37155c909
md"""
The branch obtained contain information about the bifurcation points `br.specialpoint` and `get_normal_form` can compute the Normal Form at that points.
"""

# ╔═╡ 4b0d25fc-551b-48e6-b4db-154c12b8391f
md"""
## Continuation of the Fold (Saddle Node) > CUSP
"""

# ╔═╡ ffb2f68b-493e-4ae9-8626-3eb5a49b0e17
begin
	# we generate a new branch for a fixed value of μ2
	prob1b = BK.BifurcationProblem(takens3, u0, set(p,μ2,-0.1), μ1, record_from_solution = (x, p) -> (x = x[1], y = x[2]))
	br1b = BK.continuation(prob1b, BK.PALC(tangent=BK.Bordered()),opts_br)
	opts_br2 = BK.ContinuationPar(p_min=-0.4,p_max=0.4, ds = -0.001, dsmax = 0.02,n_inversion=6)
	# continuation of critical point 
	br2 = @time BK.continuation(br1b, 1, μ2, opts_br2,detect_codim2_bifurcation=2,bdlinsolver=BK.MatrixBLS(),update_minaug_every_step=1,start_with_eigen=true)
	scene2 = plot(br2);
end

# ╔═╡ ccf38f3c-d73e-4174-84bd-1399d64982f4
md"""
We first generate a new branch (br1b) for a fixed value of μ2=-0.1 as a starting point (-0.0841231,-0.1) in the parameter space (μ1,μ2), and then give that branch as an initial guess for a Fold point to continuation in order to calculate a curve of Fold points in parameter space based on a Minimally Augmented formulation. We have to indicate index of the critical point (1) as a second argument, the second parameter to be varied (μ2) and the options. 
In this case the options also include the minimum and maximum value for the parameter to be varied and the initial step ds (note that in this case we are also going backwards at first).
The Matrix based Bordered Linear Solver is used as the bordered linear solver for the constraint equation (recommended for ODEs)
`detect_codim2_bifurcation=2` serves to detect Cusp/Bogdanov-Takens/Bautin codimension 2 critical points precisely.
And finally `update_minaug_every_step=1` update vectors a, b in Minimally Formulation at each step.
"""

# ╔═╡ 29b9b64a-19fa-470f-9596-45f67b30b983
cusp = BK.get_normal_form(br2, 1)

# ╔═╡ 831b2733-fb18-451a-8261-1a30bf71298f
bt = BK.get_normal_form(br2, 2; lens=μ2, nev = 2, autodiff = false)

# ╔═╡ 41f4608f-fad2-4202-a9be-c0a3dc5a379b
md"""
## Continuation of the Hopf > BOGDANOV TAKENS
"""

# ╔═╡ b5e04981-da64-4e1d-b7d5-e9a34028b26b
begin
	# we generate a new branch for another fixed value of μ2
	prob1c = BK.BifurcationProblem(takens3, u0, set(p,μ2,-0.3), μ1, record_from_solution = (x, p) -> (x = x[1], y = x[2]))
	br1c = BK.continuation(prob1c, BK.PALC(tangent=BK.Bordered()),opts_br)
	opts_br3 = BK.ContinuationPar(p_min=-0.3,p_max=0.01, ds = 0.001, dsmax = 0.02,n_inversion=6)
	# continuation of critical point 3 (Hopf)
	br3 = @time BK.continuation(br1c, 3, μ2, opts_br3,detect_codim2_bifurcation=2,update_minaug_every_step = 1)
end

# ╔═╡ 7083a569-b724-4281-a7d9-0074a357f541
md"""
This case is similar to the previous one with the exception that now we start with a initial guess of the Hopf for μ2=-0.3 and go forward.
"""

# ╔═╡ b143e725-7d5f-45e1-ad32-523766b579bc
begin
	plot(br2,branchlabel="SN")
	plot!(br3,branchlabel="Hopf")
end	

# ╔═╡ 311c4ba1-be7f-465a-b38d-bbb635305ca1
br3.specialpoint[1]

# ╔═╡ 79d047e0-0660-4b27-b1c7-0c24f46ed1d4
md"""
## Crossing the Hopf for a fixed μ2

We fix μ2=-0.16, start the equilibrium branch at μ1=0.01, and continue the periodic orbits born at the Hopf point while decreasing μ1. Their period initially grows and then decreases on this slice. A periodic orbit approaching a homoclinic connection would instead have a diverging period.

We first compute the branch for this trajectory in parameter space as before
"""

# ╔═╡ c6175b02-6711-45c8-b499-8ef8a0a2b82e
begin
	p2 = (μ1=0.01, μ2=-0.16)
	prob1d = BK.BifurcationProblem(takens3, u0, p2, μ1, record_from_solution = (x, p) -> (x = x[1], y = x[2]))
	opts_br1d = BK.ContinuationPar(p_min=-0.2,p_max=p2[:μ1], ds = -0.001, dsmax = 0.02, detect_bifurcation=3, n_inversion=8)
	br1d = BK.continuation(prob1d, BK.PALC(tangent=BK.Bordered()),opts_br1d)
end	

# ╔═╡ 548518de-0251-4bee-be0c-06a96d952e43
md"""
## Continuation of Periodic Orbits using the Trapezoid Method

The options are the same than those used for the branch `br1d`, adding tolerance and maximum step options

"""

# ╔═╡ 77cdc1b9-d115-4bc5-aa00-5b28b3bdfae5
opts_po = setproperties(opts_br1d, max_steps = 150, tol_stability = 1e-8);

# ╔═╡ 6152a0c0-4900-473c-8528-54d04b08eb0f
md"""
Then we call continuation with the computed branch and the index corresponding to the Hopf critical point and a functional that is encoded in the composite type `PeriodicOrbitTrapProblem`
"""

# ╔═╡ e3039d9c-b06b-47ba-9262-014b7eba5b2a
br_pot = @time BK.continuation(br1d, 3, opts_po, BK.PeriodicOrbitTrapProblem(M = 150));

# ╔═╡ 8ce52366-7bea-4019-89d3-dec24648dde9
md"""
## Periodic orbits with Parallel Standard Shooting

For this method we need to provide the ODE problem because the PSS is based on finding a solution to the flow $\Phi^T(x)=x$ for some period T, where $\Phi^t(x)$ is the flow of the system at time t. 
"""

# ╔═╡ 5af26e47-3a0f-4c52-9a9b-803093652531
odeprob = ODEProblem(takens3!, copy(u0), (0., 1.), p2; abstol = 1e-10, reltol = 1e-9);

# ╔═╡ f9681d1d-ac9a-4186-b654-77c991bcb5fa
md"""
This is similar to the previous one with the exception of the functional and some additional parameters to ensure convergence.
"""

# ╔═╡ 660fba6e-f7e0-4f42-a572-3507afd24a7e
br_pos = @time BK.continuation(br1d, 3, opts_po, BK.ShootingProblem(35, odeprob, Tsit5(), parallel = true); δp = 0.0001);

# ╔═╡ f96a46fc-e2a0-413f-a04c-e083ae9b5d3c
begin
	plot(br_pot, vars=(:param,:period), label="Trapezoidal", xlabel="μ₁", ylabel="period")
	plot!(br_pos, vars=(:param,:period), label="Standard Shooting")
end

# ╔═╡ 96e1eac9-5956-48b5-a76c-cfd5c7644f12
md"""
Period growth alone does not establish a homoclinic bifurcation. Below we compute the homoclinic curve directly by varying both parameters.
"""

# ╔═╡ 699bb082-68c0-446f-b11e-7281861b7f13
md"""
## Branch of homoclinic orbits with Orthogonal Collocation
"""

# ╔═╡ dce93919-7291-46b4-83cb-10ac297a502c
#odeprob_hc = ODEProblem(takens3!, copy(u0), (0., 1.), set(p,μ2,-0.2); abstol = 1e-10, reltol = 1e-9);

# ╔═╡ c053ce00-e4bb-4f41-8323-b9ce9bc61e3f
prob_hc = BK.BifurcationProblem(takens3, u0, set(p,μ2,0.01), μ1, record_from_solution = (x, p) -> (x = x[1], y = x[2]))

# ╔═╡ 1c2a4ece-df85-4dbb-8cc3-4a293d724b61
opt_new = BK.NewtonPar(verbose = true, tol = 1e-9, max_iterations = 12)

# ╔═╡ 48829aee-b1ad-4382-b0be-bde0a9259241
opt_hc = BK.ContinuationPar(newton_options = opt_new, max_steps = 300, save_sol_every_step = 1, dsmax = 0.001, p_min = -0.157, ds = -0.0001, p_max=0.001, dsmin = 1e-5,
	detect_event = 0, detect_bifurcation = 0);

# ╔═╡ 3599952e-d26b-4a66-856b-85ea05759a1f
br_hom_c = @time BK.continuation(prob_hc, bt, BK.PeriodicOrbitOCollProblem(20, 3; meshadapt = true, K = 100), BK.PALC(tangent = BK.Bordered()), opt_hc, 
		 ϵ0 = 1e-5, amplitude = 0.01,  freeparams = ((@lens _.ϵ0), (@lens _.ϵ1)), update_every_step = 1)

# ╔═╡ 61088c17-f8fd-4404-80e2-c9e911d5a858
begin
	plot(br2,branchlabel="SN")
	plot!(br3,branchlabel="Hopf")
	plot!(br_hom_c,branchlabel="Homoclinic")
end	

# ╔═╡ 4bb8c192-7dc0-4780-acfb-34f189a3df4a
# Esto es para ensancha la caja por defecto
html"""
<style>
	main {
		margin: 0 auto;
		max-width: 1800px;
    	padding-left: max(160px, 10%);
    	padding-right: max(160px, 10%);
	}
input[type*="range"] {
	width: 40%;
}
</style>
"""

# ╔═╡ Cell order:
# ╠═5e1d221c-18f7-11ee-20bc-b5624af0581e
# ╠═3ef551ed-71c8-4c56-8a14-46e70adbae04
# ╟─9b3c7bab-5036-5f1b-9e44-73254ea506b0
# ╟─385b660a-05e8-56ab-a144-61c90eac704a
# ╠═c285bb10-108c-5a51-8ca4-27756600b444
# ╟─26ee2eab-8184-5def-bc9a-3c254d018f82
# ╠═671a891a-34bf-5470-9eec-686511a63990
# ╠═77f8aa93-cfba-5b38-83f0-b1822d352a52
# ╠═290a3ea8-96a6-5584-985a-b2195a32f96e
# ╠═4c68fded-6f1a-51ca-83e7-45873f5f2208
# ╠═5f962302-4dee-5326-a383-e3af9e5e8344
# ╠═399a4606-410e-5f3d-ba07-3d1b05d9cb1e
# ╟─eb629797-35e8-5b4d-870a-389e02eb5a83
# ╠═f723d036-5c8b-5c58-9870-25195741cca9
# ╟─cb10ff3c-130a-50a6-beb7-10a201a8bd09
# ╠═328597eb-805f-59f4-9c6d-abc023133e94
# ╠═6596bdf7-e667-5176-9281-08e923e4cef7
# ╠═cf7b5eaf-43a9-5474-985e-96f157055c4f
# ╠═7406136b-b6e4-5e76-9a68-fa37dbf9bde2
# ╠═1421e8c7-c81e-54f3-a52d-547881d2f819
# ╠═d2cd519a-da55-54f2-a2c5-7bf529c98a68
# ╟─b206a81e-3460-5737-b159-03f6a9d2c0e6
# ╠═cd2f71d3-d79b-5241-a1dd-33b8dc0d86ea
# ╟─ea97d7be-4714-5c60-a3df-ceb787e2a480
# ╠═2ff23db5-491d-508e-8c30-5594a2298f1b
# ╠═99c9a459-13d2-5c4a-898b-409ccc88e235
# ╠═3b2de99d-f775-51ff-8d5b-c22d775e9036
# ╠═bc965a43-02ca-5c83-a79e-e2bfe8f684b4
# ╠═76730df7-e671-534b-a62f-51376d9a3c9b
# ╠═42cc2966-ad31-52cb-9e80-69def66654e0
# ╠═00847c09-a1ad-5a87-bb7c-93fcec283554
# ╟─55f5802d-3f2a-5ade-8a59-dd3534240ded
# ╠═f885398f-6fae-50f8-aace-b18898efe4d8
# ╟─dba87d4c-5ba4-5c23-a31b-e23c98ef3521
# ╠═8dcfcf37-5f29-5769-aeee-1459b67730b6
# ╠═2343f559-74e5-56dc-94cf-1899ea356c75
# ╠═bb137189-06c4-5d68-8172-64a068f631a7
# ╠═a01be38b-e781-5408-a2c9-1a24ac0e5503
# ╠═b0128688-d827-5dd6-9b66-a869bd9739a0
# ╠═602a7743-3002-5a7d-b4ae-d385a955ea7b
# ╠═0afbb419-a46f-516b-9de9-9c82b76eda44
# ╟─93e37887-478d-5e70-87f8-475519caa4cc
# ╟─a9f3eb77-a4f6-5d1d-91e4-d82d0dbaf6d2
# ╠═c17a6aa4-91a4-5e15-bec7-c14f5f5e599c
# ╟─14401898-c1f6-5f2e-84ed-81ca7aa77256
# ╠═6d8fed6b-7f23-5b71-8e60-ed406cff9ad2
# ╠═f2d5e1ad-0f0f-51ac-bd8b-5e6cf413ae28
# ╠═01058bac-b766-50dc-b583-e2a28d762c50
# ╠═7fe71b73-6632-5c3c-a36f-f0e04e9812f9
# ╠═9aba9f01-e96e-558a-bfd2-92d4ffc06d6c
# ╠═a5d0bbb6-ff06-5087-b733-612965cd0cec
# ╟─a6a935c2-cb82-53f6-b13e-b76230860839
# ╟─77fd29b5-4cc6-4461-ab16-3f0720643dee
# ╠═4380968c-9421-4c34-82bd-850c3027b676
# ╟─e65f74c0-edf0-4378-ab3d-5e848bda5d2f
# ╠═b4574225-ff15-44f5-82e2-5eb779a16c49
# ╠═f3ad61c9-33de-45b0-bee5-89e5f86110e7
# ╠═6ed8af70-4918-4dc2-b419-d87c8f35cfb9
# ╠═90461200-bc32-4f01-ac06-a7e9b81966df
# ╟─9eff4af0-8cef-4a9a-9e6c-092a1e046b96
# ╠═fa557b12-9017-4166-ac5b-ead0c3e108b0
# ╟─1b94418e-60b7-4da4-8720-4ac37155c909
# ╟─4b0d25fc-551b-48e6-b4db-154c12b8391f
# ╠═ffb2f68b-493e-4ae9-8626-3eb5a49b0e17
# ╟─ccf38f3c-d73e-4174-84bd-1399d64982f4
# ╠═29b9b64a-19fa-470f-9596-45f67b30b983
# ╠═831b2733-fb18-451a-8261-1a30bf71298f
# ╟─41f4608f-fad2-4202-a9be-c0a3dc5a379b
# ╠═b5e04981-da64-4e1d-b7d5-e9a34028b26b
# ╟─7083a569-b724-4281-a7d9-0074a357f541
# ╠═b143e725-7d5f-45e1-ad32-523766b579bc
# ╠═311c4ba1-be7f-465a-b38d-bbb635305ca1
# ╟─79d047e0-0660-4b27-b1c7-0c24f46ed1d4
# ╠═c6175b02-6711-45c8-b499-8ef8a0a2b82e
# ╟─548518de-0251-4bee-be0c-06a96d952e43
# ╠═77cdc1b9-d115-4bc5-aa00-5b28b3bdfae5
# ╟─6152a0c0-4900-473c-8528-54d04b08eb0f
# ╠═e3039d9c-b06b-47ba-9262-014b7eba5b2a
# ╟─8ce52366-7bea-4019-89d3-dec24648dde9
# ╠═5af26e47-3a0f-4c52-9a9b-803093652531
# ╟─f9681d1d-ac9a-4186-b654-77c991bcb5fa
# ╠═660fba6e-f7e0-4f42-a572-3507afd24a7e
# ╠═f96a46fc-e2a0-413f-a04c-e083ae9b5d3c
# ╟─96e1eac9-5956-48b5-a76c-cfd5c7644f12
# ╟─699bb082-68c0-446f-b11e-7281861b7f13
# ╠═dce93919-7291-46b4-83cb-10ac297a502c
# ╠═c053ce00-e4bb-4f41-8323-b9ce9bc61e3f
# ╠═1c2a4ece-df85-4dbb-8cc3-4a293d724b61
# ╠═48829aee-b1ad-4382-b0be-bde0a9259241
# ╠═3599952e-d26b-4a66-856b-85ea05759a1f
# ╠═61088c17-f8fd-4404-80e2-c9e911d5a858
# ╟─4bb8c192-7dc0-4780-acfb-34f189a3df4a
