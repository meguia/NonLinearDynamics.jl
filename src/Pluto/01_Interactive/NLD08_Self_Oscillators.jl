### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 1097d874-4468-4e76-bdbc-c893a5dbfdc0
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using Plots, PlutoUI, DifferentialEquations, ForwardDiff, StaticArrays, IntervalRootFinding
    using NonLinearDynamics
end

# ╔═╡ 18a83e00-0dd4-4fe7-a5be-644361f875d3
TableOfContents()

# ╔═╡ d263b7d9-a736-480c-b401-813e6dbb41ca
md"""

# Self-oscillator: Simple reed model (Rayleigh)


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

# ╔═╡ 3817cc7f-e32f-4323-8ffb-5ef275f76f9d
function rossler!(du,u,p,t)
    (a,b,c,sc) = p
	(x,y,z) = [sc*u[1],sc*u[2],sc*u[3]]
	du[1] = -y-z*2
	du[2] = x + a*y
	du[3] = b + 2*z*(x-c)
end

# ╔═╡ 07a719c4-7a00-4260-ae4b-5d0f18a525ec
begin
	sol3 = solve(ODEProblem(rossler!,[0.03,0,0],(0,2),[0.2,0.2, 5.0,200.0]))
	sol3b = solve(ODEProblem(rossler!,[0.01,0.02,0.01],(0,2),[0.2,0.2, 5.0,200.0]))
	plot(sol3,idxs=(1,2,3))
	plot!(sol3b,idxs=(1,2,3))
	plot!(xlims=(-0.05,0.05),ylims=(-0.05,0.05),zlims=(0.0,0.05))
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
# Rubbed Oscillator (bowed string)
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

# ╔═╡ 5bbb333e-47a3-464b-82f5-cb8053c5c733
function plot_state_flow2D!(p1,u,t;alpha=0.5)
	title = "t=$(round(t,digits=4))  x(t) = $(round(u[1],digits=6)) y(t) = $(round(u[2],digits=6))"
	scatter!(p1,u[1:1],u[2:2],c=:red,alpha=alpha, markerstrokewidth = 0)
	title!(p1,title)
end;

# ╔═╡ 64022ca2-6323-4428-b938-843e228d4922
function plot_state_flow2D(u,xrange,yrange,t)
	title = "t=$(round(t,digits=4))  x(t) = $(round(u[1],digits=6)) y(t) = $(round(u[2],digits=6))"
	p1 = plot(xlims=xrange,ylims=yrange,size=(800,800),title=title,legend=false)
	plot!(p1,[xrange[1],xrange[2]],[0,0],c=:black)
	plot!(p1,[0,0],[yrange[1],yrange[2]],c=:black)
end;

# ╔═╡ Cell order:
# ╠═1097d874-4468-4e76-bdbc-c893a5dbfdc0
# ╟─18a83e00-0dd4-4fe7-a5be-644361f875d3
# ╠═d263b7d9-a736-480c-b401-813e6dbb41ca
# ╠═2cd075b8-0fd1-481b-b230-1bf3a51f2d5f
# ╟─7557f8d7-2654-41c2-a443-9d5a76622501
# ╠═6e5fe124-a757-4d25-a270-26941b71ece3
# ╠═3817cc7f-e32f-4323-8ffb-5ef275f76f9d
# ╠═07a719c4-7a00-4260-ae4b-5d0f18a525ec
# ╠═9522e566-d087-4ea4-8f37-438dd33ac250
# ╠═c45e6d60-3086-42f5-ba1b-df6de93827ed
# ╟─a6624901-991d-4053-9270-ba1b54ec48ec
# ╠═e94dd082-b9f5-4bc4-a530-37459881dcd0
# ╟─2ba8419d-63bd-44cb-a6cf-b72e57c2494b
# ╟─9fc16286-ef4f-4808-8227-0e3f69bbbaec
# ╠═3f61ca09-63d8-45ee-8cc5-f4296d165f17
# ╠═1f53d7e2-569d-4b90-94ea-e317228a6efd
# ╠═ad89473b-73a5-4c3e-a3f9-f645e19b386b
# ╠═8c710c72-6c9b-40ea-bd7c-dcc566a317e0
# ╠═82365a14-9870-43bc-b214-b1d248c5cd9d
# ╟─f480e291-b67c-4a2f-9f44-914a340df81e
# ╟─5bbb333e-47a3-464b-82f5-cb8053c5c733
# ╟─64022ca2-6323-4428-b938-843e228d4922
