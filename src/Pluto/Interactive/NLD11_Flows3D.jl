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

# ╔═╡ 6a3d6236-a28d-4040-aa15-d5ef8ddbe3ed
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using PlutoUI, Plots, DifferentialEquations, ForwardDiff, IntervalRootFinding, StaticArrays
    using LinearAlgebra
end

# ╔═╡ 5b2eef3e-cd79-4989-8c94-e2dcc215bc2a
md"""
# Lorenz System

This is undoubtedly the most famous dynamical system featuring chaos. The origin goes back to an original work done by Edward Lorenz in 1963, proposing a reduction to three differential equations of an already simplified model of atmospheric convection (the original model had 12 equations).

In the model, the lower layer of the atmosphere is at a higher temperature than the upper layer. In the case where the temperature difference is small there is a temperature gradient (a linear variation of temperature with height) but above a certain critical value the hot air rises and the cold air falls forming convection rolls as seen in the figure.
"""

# ╔═╡ d7c7c3dd-8eab-4a9d-b372-804a8b2c9477
html"""
<div>
<img src="https://i.imgur.com/ukyDvpM.png" width="700px">
</div>
"""

# ╔═╡ c74fb54e-b49b-4687-8c0d-5ada8769ddbf
md"""
For a history of how Lorenz discovered the sensitivity to the initial conditions of this system from a truncation of the numerical simulation printed on paper and the true origin of the term "butterfly effect" see

$\dot{x} = \sigma(y-x)$

$\dot{y} = \rho x - y - xz$

$\dot{z} = xy - \beta z$

In this system, the variable $x$ corresponds to the intensity of convection (how fast the coils rotate), the variable $y$ to the temperature difference between the rising hot air stream and the falling cold one, and $z$ to the deviation of the linear temperature variation with height. The parameters also have physical significance, although there are not ver common outside the domain of fluid dynamics: $\sigma$ is the Prandtl number, $\rho$ the Rayleigh number and $\beta$ a geometric factor.

In order to obtain chaos, the values $\sigma=10$, $\rho=28$ and $\beta=8/3$ are traditionally used, although these values do not correspond to any particular physical system. For high Rayleigh number values $\rho>1$ the truncation made by Lorenz to only three modes is no longer valid. Therefore, although the model can explain the origin of the convection rolls (for $\rho=1$) the chaotic regime does not correspond to the physical model. However, being a relatively simple system to analyze and mainly for historical reasons, the Lorenz model became the "model of models" displaying chaotic behavior.

"""

# ╔═╡ 80630ff0-64d2-49c0-b879-cf5de0d8fad0
function lorenz!(du,u,p,t)
    (σ,ρ,β)=p
    du[1]=σ*(u[2]-u[1])
    du[2]=ρ*u[1]-u[2]-u[1]*u[3]
    du[3]=u[1]*u[2]-β*u[3]
    du
end    

# ╔═╡ 4de8b3eb-4a77-47d1-b7cc-7a2d5c5bfd8d
gr();

# ╔═╡ b1cce79e-495e-46c3-9d7e-cb0869cada3a
@bind parl (
	PlutoUI.combine() do bind
		md"""
		σ: $(bind(Slider(0:0.2:20.0,default=10.0;show_value=true))) \
		ρ: $(bind(Slider(0.0:0.2:30.0,default=28.0;show_value=true))) \
		β: $(bind(Slider(sort([collect(0:0.02:3.0);8/3]),default=8/3;show_value=true))) \
		"""
	end
)	

# ╔═╡ c2338e97-1c4c-4b61-bb83-794a54db455f
sol = solve(ODEProblem(lorenz!,[0.1,0.1,5.0],(0.0,100),parl));

# ╔═╡ 59b83e8b-75f6-48bd-80d0-e3eeedccbdfc
sol(10.0)

# ╔═╡ 40597d23-a13c-4954-8fca-d6423e6b36a0
plot(sol,idxs=(1,2,3),label="lorenz")

# ╔═╡ 012607ea-b56e-4a78-909b-fb982fe4a06c
20/0.04

# ╔═╡ 3592613a-538c-49cd-a27e-65d41b8b9f5d
begin 
	p1 = plot(sol,idxs=(1,2,3),legend=false,title="Lorenz")
	p2 = plot(sol,idxs=(0,1),label="x")
	p3 = plot(sol,idxs=(0,2),label="y")
	p4 = plot(sol,idxs=(0,3),label="z")
	plot(p1,p2,p3,p4,layout=@layout [a{0.5w} grid(3,1)])
end	

# ╔═╡ 5dc493e1-51b7-51ed-99f7-89a021057154
md"""
## Equilibria organize the three-dimensional flow

In the convection interpretation, $x$ measures circulation, $y$ a horizontal
temperature contrast, and $z$ a vertical temperature-profile correction.
The three-mode truncation is an instructive dynamical system, not a complete
weather model.

For $\sigma,\beta>0$, the origin is stable for $\rho<1$. For $\rho>1$,
two equilibria appear:

```math
C_\pm=(\pm\sqrt{\beta(\rho-1)},\ \pm\sqrt{\beta(\rho-1)},\ \rho-1),\qquad
J=\begin{pmatrix}
-\sigma&\sigma&0\\ \rho-z&-1&-x\\ y&x&-\beta
\end{pmatrix}.
```

For $\sigma>\beta+1$, their linear Hopf threshold is
$\rho_H=\sigma(\sigma+\beta+3)/(\sigma-\beta-1)$, about $24.74$ for
$\sigma=10$, $\beta=8/3$. This local threshold does **not** locate all
global transitions or guarantee chaos at every larger $\rho$.
"""

# ╔═╡ 653325b0-8bc0-592b-ad81-2f7218aee293
let
    σ,ρ,β=collect(parl)
    points=[[0.0,0.0,0.0]]
    if ρ>1 && β>0
        q=sqrt(β*(ρ-1))
        append!(points,[[q,q,ρ-1],[-q,-q,ρ-1]])
    end
    [(point=u,eigenvalues=eigvals([-σ σ 0;ρ-u[3] -1 -u[1];u[2] u[1] -β])) for u in points]
end

# ╔═╡ 87355113-e467-5e3a-9ea7-86868ccbf295
md"""
## Nearby initial states and the prediction horizon

Integrate the same equations from $\mathbf u_0$ and
$\mathbf u_0+(10^{-7},0,0)$ with the same tight tolerances. Compare the
full-state separation $d(t)=\|\mathbf u_1(t)-\mathbf u_2(t)\|$ on a log scale.
In a chaotic regime, an initial growth interval is followed by saturation at
the attractor's scale. A finite-time separation curve is not by itself a
Lyapunov-exponent calculation; also check tolerance sensitivity.
"""

# ╔═╡ ad03e43d-da9f-5b0c-9f60-8d912c3da54f
course_lorenz_u0 = [1.0,0.0,0.0]

# ╔═╡ 4b701899-d9a6-5996-af4c-cdd347a1d50d
course_lorenz_tspan = (0.0,60.0)

# ╔═╡ e382f9ac-3bb7-500d-866c-42a7f4b628cf
course_lorenz_p = collect(parl)

# ╔═╡ 53553469-e156-5b98-bd76-bd0aa0294bf1
course_lorenz_prob = ODEProblem(lorenz!, course_lorenz_u0, course_lorenz_tspan, course_lorenz_p)

# ╔═╡ 4b84aa2c-741d-5576-8682-3ff446e5c631
course_lorenz_sol = solve(course_lorenz_prob, Tsit5(); abstol=1e-9, reltol=1e-8);

# ╔═╡ f926d0d2-5be7-5629-abb3-ef1bafdc0b5e
course_lorenz_near = solve(remake(course_lorenz_prob;u0=course_lorenz_u0+[1e-7,0,0]),
    Tsit5();abstol=1e-9,reltol=1e-8);

# ╔═╡ 55005b74-bec9-587a-8979-3ac00279c38b
let ts=range(course_lorenz_tspan...;length=3001)
    distance=[norm(course_lorenz_sol(t)-course_lorenz_near(t)) for t in ts]
    a=plot(course_lorenz_sol;idxs=(0,1),label="reference",xlabel="t",ylabel="x")
    plot!(a,course_lorenz_near;idxs=(0,1),label="perturbed",alpha=0.7)
    b=plot(ts,max.(distance,eps());yscale=:log10,xlabel="t",ylabel="‖δu(t)‖",legend=false)
    plot(a,b;layout=(2,1),size=(900,500),margin=5*Plots.mm)
end

# ╔═╡ bdda09ed-86f3-54d9-a53f-2217f4116cab
md"""
## Reduce a flow to a return map

Record successive maxima of $z(t)$ after the transient. A maximum occurs
when $\dot z=xy-\beta z$ crosses zero **from positive to negative**.
The callback locates that event between solver steps; sampling a coarse
time grid can miss it. Plot $z_{n+1}$ against $z_n$.

At the classical chaotic parameters this projection gives an almost
one-dimensional folded relation. It illustrates stretching and folding,
but is a projection of a return map rather than an exact scalar evolution law.
"""

# ╔═╡ 2e90b391-e3b8-58f7-94f0-ed6498bdefec
course_lorenz_maxima = let times=Float64[], maxima=Float64[]
    callback=ContinuousCallback(
        (u,t,i)->u[1]*u[2]-i.p[3]*u[3],nothing,
        i->begin
            if i.t>30
                push!(times,i.t)
                push!(maxima,i.u[3])
            end
        end;save_positions=(false,false))
    solution=solve(remake(course_lorenz_prob;tspan=(0.0,200.0)),Tsit5();
        callback,abstol=1e-9,reltol=1e-8,dtmax=0.05,save_everystep=false)
    (times=times,values=maxima,retcode=solution.retcode)
end

# ╔═╡ 9c8d1746-d359-555e-a838-54af36c40cb7
let z=course_lorenz_maxima.values
    if length(z)>=3 && maximum(z)-minimum(z)>1e-5
        f=scatter(z[1:end-1],z[2:end];xlabel="zₙ (maximum)",ylabel="zₙ₊₁",
            markersize=2,markerstrokewidth=0,legend=false,size=(550,450),margin=5*Plots.mm)
        plot!(f,collect(extrema(z)),collect(extrema(z));color=:gray,linestyle=:dash)
    else
        md"No resolved sequence of distinct maxima at these parameters; try σ = 10, ρ = 28, β = 8/3."
    end
end

# ╔═╡ 85b1b39d-98ae-5c75-8b8f-02cf69182337
md"""
**Try:** compare $\rho=0.5$, $10$, and $28$ with $\sigma=10$, $\beta=8/3$.
Explain the final behavior using equilibria, the separation curve, and the
maxima map together. Near bifurcations, extend the transient before deciding
that irregular motion is a persistent attractor.
"""

# ╔═╡ Cell order:
# ╠═6a3d6236-a28d-4040-aa15-d5ef8ddbe3ed
# ╟─5b2eef3e-cd79-4989-8c94-e2dcc215bc2a
# ╟─d7c7c3dd-8eab-4a9d-b372-804a8b2c9477
# ╟─c74fb54e-b49b-4687-8c0d-5ada8769ddbf
# ╠═80630ff0-64d2-49c0-b879-cf5de0d8fad0
# ╠═c2338e97-1c4c-4b61-bb83-794a54db455f
# ╠═59b83e8b-75f6-48bd-80d0-e3eeedccbdfc
# ╠═4de8b3eb-4a77-47d1-b7cc-7a2d5c5bfd8d
# ╠═40597d23-a13c-4954-8fca-d6423e6b36a0
# ╠═b1cce79e-495e-46c3-9d7e-cb0869cada3a
# ╠═012607ea-b56e-4a78-909b-fb982fe4a06c
# ╠═3592613a-538c-49cd-a27e-65d41b8b9f5d
# ╟─5dc493e1-51b7-51ed-99f7-89a021057154
# ╠═653325b0-8bc0-592b-ad81-2f7218aee293
# ╟─87355113-e467-5e3a-9ea7-86868ccbf295
# ╠═ad03e43d-da9f-5b0c-9f60-8d912c3da54f
# ╠═4b701899-d9a6-5996-af4c-cdd347a1d50d
# ╠═e382f9ac-3bb7-500d-866c-42a7f4b628cf
# ╠═53553469-e156-5b98-bd76-bd0aa0294bf1
# ╠═4b84aa2c-741d-5576-8682-3ff446e5c631
# ╠═f926d0d2-5be7-5629-abb3-ef1bafdc0b5e
# ╠═55005b74-bec9-587a-8979-3ac00279c38b
# ╟─bdda09ed-86f3-54d9-a53f-2217f4116cab
# ╠═2e90b391-e3b8-58f7-94f0-ed6498bdefec
# ╠═9c8d1746-d359-555e-a838-54af36c40cb7
# ╟─85b1b39d-98ae-5c75-8b8f-02cf69182337
