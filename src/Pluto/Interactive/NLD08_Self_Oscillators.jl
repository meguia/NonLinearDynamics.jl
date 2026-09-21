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

# ╔═╡ 1097d874-4468-4e76-bdbc-c893a5dbfdc0
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using NonLinearDynamics: flow2d_vectorfield
    using Plots, PlutoUI, DifferentialEquations, ForwardDiff, StaticArrays, IntervalRootFinding
end

# ╔═╡ 18a83e00-0dd4-4fe7-a5be-644361f875d3
TableOfContents()

# ╔═╡ 39c95198-781e-5899-9155-1584dcbaf093
md"""
# Self-oscillation and excitability

Negative damping can supply energy at small amplitude while nonlinear
dissipation limits large motion. Follow the circuit oscillator into its
slow–fast description, then compare reed, friction, and excitable systems.
"""

# ╔═╡ bd7a55af-7de5-5ce6-a440-3e7c48e87465
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

# ╔═╡ 16bff63a-b6b8-5947-9292-86a40780b98b
function course_vdp!(du,u,p,t)
    μ=only(p)
    x,v=u
    du[1]=v
    du[2]=μ*(1-x^2)*v-x
    nothing
end

# ╔═╡ a154f530-73f1-5ef3-acc3-95b24e9ab59a
md"""
μ: $(@bind course_vdp_mu Slider(0.2:0.2:8.0; default=3.0, show_value=true))
"""

# ╔═╡ 690bf67b-504f-5ca9-8257-a93f2e7e8852
course_vdp_u0 = [0.1,0.0]

# ╔═╡ 31d6139e-6bf8-5879-8f36-1330e8016ccc
course_vdp_tspan = (0.0,80.0)

# ╔═╡ f0aee207-34e6-5bf0-9bd4-e27d2ce0c053
course_vdp_p = [course_vdp_mu]

# ╔═╡ 8ad23053-4c21-54b1-975d-b1e1f8eff977
course_vdp_prob = ODEProblem(course_vdp!, course_vdp_u0, course_vdp_tspan, course_vdp_p)

# ╔═╡ 8e2722b0-9301-51c7-95b2-e6617be43b81
course_vdp_sol = solve(course_vdp_prob, Tsit5(); abstol=1e-9, reltol=1e-8);

# ╔═╡ 1f810b57-70df-54ec-af77-ff24c9dfad8c
let a=plot(course_vdp_sol;idxs=(0,1),xlabel="t",ylabel="x",legend=false),
    b=plot(course_vdp_sol;idxs=(1,2),xlabel="x",ylabel="v",legend=false)
    plot(a,b;layout=(1,2),size=(900,370),margin=5*Plots.mm)
end

# ╔═╡ 42a2b8b7-965e-5967-b183-da428d6639d3
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

# ╔═╡ 1552efba-5bc0-5c43-a15d-9bd3ac90a5f5
function course_lienard!(du,u,p,t)
    μ=only(p)
    x,z=u
    du[1]=μ*(x-x^3/3-z)
    du[2]=x/μ
    nothing
end

# ╔═╡ d735e88b-0f3a-5cd3-b8f2-64257c5328b5
course_lienard_u0 = [course_vdp_u0[1],course_vdp_u0[1]-course_vdp_u0[1]^3/3-course_vdp_u0[2]/course_vdp_mu]

# ╔═╡ 64217aad-a2a5-5844-83ee-4706ede0b6c3
course_lienard_tspan = course_vdp_tspan

# ╔═╡ 9c403a3b-77ef-529a-9dea-131c4effce8c
course_lienard_p = course_vdp_p

# ╔═╡ 2538ca1c-b851-54f9-afc6-c90939e7809f
course_lienard_prob = ODEProblem(course_lienard!, course_lienard_u0, course_lienard_tspan, course_lienard_p)

# ╔═╡ 4ca69de6-e260-5269-ba5d-43c95d09f33a
course_lienard_sol = solve(course_lienard_prob, Tsit5(); abstol=1e-9, reltol=1e-8);

# ╔═╡ 416400d1-6545-55f0-ab43-3d466bd492d3
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

# ╔═╡ d5063a87-48df-59e6-9b40-ce55cc76435b
md"""
**Try:** compare $\mu=0.2$, $3$, and $8$. Where does the motion spend most of
its time? In the standard van der Pol family, $\mu=0$ removes the nonlinear
damping entirely and leaves a center. It is a degenerate onset, not the
generic small-amplitude Hopf bifurcation studied in NLD12; the limit-cycle
amplitude tends to about $2$ as $\mu\to0^+$.
"""

# ╔═╡ d263b7d9-a736-480c-b401-813e6dbb41ca
md"""

## Self-oscillator: Simple reed model (Rayleigh)


In the "Theory of Sound" (1877) Rayleigh proposed a simple model for the blowing of a clarinet reed, in which the dissipation coefficient ($\gamma$) is negative (i.e. acts in favor of the motion and delivers energy) for small velocities and positive (i.e. acts by slowing down the motion) for large velocities. 

$\gamma = s_c^2 v^2 - \mu$

This is a very general model of a self-oscillator that can be written by replacing the expression of the nonlinear dissipation in the harmonic oscillator.

$\dot{x} = v$

$\dot{v} = -(s_c^2v^2-\mu)v-Kx$

This is a **nonlinear system** because it includes a $v^3$ term. The scale $s_c=1$ gives the model introduced in the 2D flows notebook.
"""

# ╔═╡ 2cd075b8-0fd1-481b-b230-1bf3a51f2d5f
function reed!(du,u,p,t)
    (μ,K,sc) = p
    du[1] = u[2]
	du[2] = -K*u[1]+u[2]*(μ-sc^2*u[2]^2)
end

# ╔═╡ 7557f8d7-2654-41c2-a443-9d5a76622501
md"""
x(0) $(@bind x03 Slider(-1.0:0.01:1.0,default=0.1;show_value=true)) 
K : $(@bind K2 Slider(0.1:0.1:5.0,default=1.0;show_value=true)) \
μ : $(@bind μ Slider(-0.1:0.0001:2.0,default=0.1;show_value=true))
sc : $(@bind sc Slider(1.0:0.1:100.0,default=1.0;show_value=true))\
tmax : $(@bind tmax2 Slider(5.0:5.0:40.0,default=10.0;show_value=true))
"""

# ╔═╡ 6e5fe124-a757-4d25-a270-26941b71ece3
begin
	sol1 = solve(ODEProblem(reed!,[x03,0],(0,tmax2),[μ,K2,sc]))
	sol2 = solve(ODEProblem(reed!,[-x03,0],(0,tmax2),[μ,K2,sc]))
	plot(sol1,idxs=(1,2))
	plot!(sol2,idxs=(1,2))
end

# ╔═╡ 9522e566-d087-4ea4-8f37-438dd33ac250
flow2d_vectorfield(reed!,[x03,0],tmax2,[μ,K2,sc];title="Simple Reed Model")

# ╔═╡ c45e6d60-3086-42f5-ba1b-df6de93827ed
begin
	sol4 = solve(ODEProblem(reed!, [x03,0], (0,tmax2), [μ,K2,sc]))
	pa4 = plot(sol4,idxs=(0,1),legend=false,xlabel="t",ylabel="x")
	pb4 = plot(sol4,idxs=(0,2),legend=false,xlabel="t",ylabel="v")
	plot!(pa4,[0,50.0],[0,0],c=:black)
	plot!(pb4,[0,50.0],[0,0],c=:black)
	plot(pa4,pb4,layout=(2,1),size=(900,400))
end	

# ╔═╡ a6624901-991d-4053-9270-ba1b54ec48ec
md"""
## Rubbed Oscillator (bowed string)
Another system with simple self-oscillations was proposed (also by Rayleigh in 1877!) to model the slip & stick action of the bow against the violin string, but it can be applied to a lot of systems that generate self-oscillations from friction.
"""

# ╔═╡ e94dd082-b9f5-4bc4-a530-37459881dcd0
html"""
<div>
<img src="https://i.imgur.com/qW4INmr.png" width="300px">
</div>
"""

# ╔═╡ 2ba8419d-63bd-44cb-a6cf-b72e57c2494b
md"""
The proposed model was similar to the one shown in the figure. A mass attached to a spring is supported on a friction conveyor belt moving with constant velocity to the right. At first the static friction causes the mass to stick (stick momentum) to the belt and exerts a force equal to that of the spring. But the static friction has a maximum value and when the spring is very stretched it cannot overcome the elastic force and the mass is dragged by the spring and slides with dynamic friction on the belt to the left (slip moment) and can reach by inertia even to compress a little the spring until the mass is braked and is hooked again by the static friction and the process is repeated.

Although it is a simple system, the functional form of the friction (which has to be a function of the difference in velocity between the mass and the belt, i.e. whether it slides or not) cannot be something as simple as a quadratic or a cubic because it has to change sign abruptly, since the friction has to be maximal for low slides and decreasing for faster slides. The classical form is something like this:
"""

# ╔═╡ 9fc16286-ef4f-4808-8227-0e3f69bbbaec
html"""
<div>
<img src="https://i.imgur.com/KrRu2Ub.png" width="200px">
</div>
"""

# ╔═╡ 3f61ca09-63d8-45ee-8cc5-f4296d165f17
md"""
where $\dot{x}-v$  is the 'slip', i.e. the difference in velocities between the mass and the belt. When the mass is attached to the tape it can take all values in the vertical up to a maximum value to one side and to the other and then 'jumps' to the dynamic friction which is smaller as the slip is faster and faster.

The arc friction proposed by Rayleigh is written as follows:

$\dot{x}=v$

$\dot{v}=-\mu C(v-V)-x$

where $\mu$ is the friction coefficient, $V$ the bow velocity, and $C$ is the friction function. Here we explore the smooth cubic approximation $C(s)=-s(1-s^2)$; it does not impose exact sticking at zero slip.
"""


# ╔═╡ 1f53d7e2-569d-4b90-94ea-e317228a6efd
#friction(x) = atan(x/0.05)*exp(-2*abs(x))
friction(x) = - x *(1-x*x)

# ╔═╡ ad89473b-73a5-4c3e-a3f9-f645e19b386b
function bow!(du,u,p,t)
    du[1]=u[2]
    du[2]=-p[1]*friction(u[2]-p[2])-u[1]
    du
end

# ╔═╡ 82365a14-9870-43bc-b214-b1d248c5cd9d
md"""
x(0) $(@bind x0b Slider(-1.0:0.02:1.0,default=0.7;show_value=true)) 
y(0) : $(@bind v0b Slider(-1.0:0.01:1.0,default=0.0;show_value=true)) \
μ : $(@bind μb Slider(-0.3:0.01:1.0,default=0.02;show_value=true)) 
V : $(@bind V Slider(-1:0.001:1,default=0.02;show_value=true))
"""

# ╔═╡ 8c710c72-6c9b-40ea-bd7c-dcc566a317e0
begin
	solbow = solve(ODEProblem(bow!, [x0b; v0b], (0, 250.0), [μb,V]));
	p1 = plot(solbow,legend=false)
	p2 = plot(solbow,idxs=(1,2),legend=false,arrow=true)
	plot(p1,p2,layout=(1,2),size = (900,450),title="Bowed String")
end	

# ╔═╡ f480e291-b67c-4a2f-9f44-914a340df81e
html"""
<style>
input[type*="range"] {
	width: 30%;
}
</style>
"""

# ╔═╡ 2179cbc6-ecb3-510b-94ac-ea9c3249560d
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

# ╔═╡ 0a86deda-f0ec-5d8f-b071-d685bbde0231
function course_fhn!(du,u,p,t)
    a,b,I=p
    x,y=u
    du[1]=x-x^3/3-y+I
    du[2]=(a*x+b-y)/10
    nothing
end

# ╔═╡ 5fce1bb1-61bc-52ab-9019-ef7f8b72bae7
md"""
I: $(@bind course_current Slider(0.0:0.02:1.2; default=0.5, show_value=true))
"""

# ╔═╡ 9af01f47-7298-5189-bc7f-bb108691e33c
course_fhn_u0 = [-1.0,-0.5]

# ╔═╡ 5c8b2673-65dd-5d65-b05a-6de15689f124
course_fhn_tspan = (0.0,200.0)

# ╔═╡ 63aab325-264c-55d2-a4cb-4c115473f3cd
course_fhn_p = [1.25,0.875,course_current]

# ╔═╡ ea9883ba-e626-5ca4-9236-3286f05feb26
course_fhn_prob = ODEProblem(course_fhn!, course_fhn_u0, course_fhn_tspan, course_fhn_p)

# ╔═╡ 30a39caa-5b12-5526-b020-2d95c38d3e28
course_fhn_sol = solve(course_fhn_prob, Tsit5(); abstol=1e-9, reltol=1e-8);

# ╔═╡ 92778bfd-947c-539a-aca6-72d9c3a98264
let xs=range(-2.5,2.5;length=401)
    a,b,I=course_fhn_p
    phase=plot(course_fhn_sol;idxs=(1,2),label="trajectory",xlabel="x",ylabel="y")
    plot!(phase,xs,xs .-xs.^3 ./3 .+I;label="dx/dt = 0",color=:red)
    plot!(phase,xs,a.*xs .+b;label="dy/dt = 0",color=:green)
    trace=plot(course_fhn_sol;idxs=(0,1),xlabel="t",ylabel="activation x",legend=false)
    plot(phase,trace;layout=(1,2),size=(950,390),margin=5*Plots.mm)
end

# ╔═╡ 80aeb7a0-8dcf-5de8-b4de-b5fd933d1f7d
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
# ╠═1097d874-4468-4e76-bdbc-c893a5dbfdc0
# ╟─18a83e00-0dd4-4fe7-a5be-644361f875d3
# ╟─39c95198-781e-5899-9155-1584dcbaf093
# ╟─bd7a55af-7de5-5ce6-a440-3e7c48e87465
# ╠═16bff63a-b6b8-5947-9292-86a40780b98b
# ╟─a154f530-73f1-5ef3-acc3-95b24e9ab59a
# ╠═690bf67b-504f-5ca9-8257-a93f2e7e8852
# ╠═1f810b57-70df-54ec-af77-ff24c9dfad8c
# ╠═31d6139e-6bf8-5879-8f36-1330e8016ccc
# ╠═f0aee207-34e6-5bf0-9bd4-e27d2ce0c053
# ╠═8ad23053-4c21-54b1-975d-b1e1f8eff977
# ╠═8e2722b0-9301-51c7-95b2-e6617be43b81
# ╟─42a2b8b7-965e-5967-b183-da428d6639d3
# ╠═1552efba-5bc0-5c43-a15d-9bd3ac90a5f5
# ╠═d735e88b-0f3a-5cd3-b8f2-64257c5328b5
# ╠═64217aad-a2a5-5844-83ee-4706ede0b6c3
# ╠═9c403a3b-77ef-529a-9dea-131c4effce8c
# ╠═2538ca1c-b851-54f9-afc6-c90939e7809f
# ╠═4ca69de6-e260-5269-ba5d-43c95d09f33a
# ╠═416400d1-6545-55f0-ab43-3d466bd492d3
# ╟─d5063a87-48df-59e6-9b40-ce55cc76435b
# ╟─d263b7d9-a736-480c-b401-813e6dbb41ca
# ╠═2cd075b8-0fd1-481b-b230-1bf3a51f2d5f
# ╟─7557f8d7-2654-41c2-a443-9d5a76622501
# ╠═6e5fe124-a757-4d25-a270-26941b71ece3
# ╠═9522e566-d087-4ea4-8f37-438dd33ac250
# ╠═c45e6d60-3086-42f5-ba1b-df6de93827ed
# ╟─a6624901-991d-4053-9270-ba1b54ec48ec
# ╟─e94dd082-b9f5-4bc4-a530-37459881dcd0
# ╟─2ba8419d-63bd-44cb-a6cf-b72e57c2494b
# ╟─9fc16286-ef4f-4808-8227-0e3f69bbbaec
# ╟─3f61ca09-63d8-45ee-8cc5-f4296d165f17
# ╠═1f53d7e2-569d-4b90-94ea-e317228a6efd
# ╠═ad89473b-73a5-4c3e-a3f9-f645e19b386b
# ╠═8c710c72-6c9b-40ea-bd7c-dcc566a317e0
# ╟─82365a14-9870-43bc-b214-b1d248c5cd9d
# ╟─f480e291-b67c-4a2f-9f44-914a340df81e
# ╟─2179cbc6-ecb3-510b-94ac-ea9c3249560d
# ╠═0a86deda-f0ec-5d8f-b071-d685bbde0231
# ╟─5fce1bb1-61bc-52ab-9019-ef7f8b72bae7
# ╠═9af01f47-7298-5189-bc7f-bb108691e33c
# ╠═5c8b2673-65dd-5d65-b05a-6de15689f124
# ╠═63aab325-264c-55d2-a4cb-4c115473f3cd
# ╠═ea9883ba-e626-5ca4-9236-3286f05feb26
# ╠═30a39caa-5b12-5526-b020-2d95c38d3e28
# ╠═92778bfd-947c-539a-aca6-72d9c3a98264
# ╟─80aeb7a0-8dcf-5de8-b4de-b5fd933d1f7d
