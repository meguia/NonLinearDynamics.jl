### A Pluto.jl notebook ###
# v0.19.9

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

# ╔═╡ 8992e664-fbdd-4b25-9a24-2c72a34d9e10
using Pkg;Pkg.add("Plots");Pkg.add("PlutoUI");Pkg.add("DifferentialEquations");

# ╔═╡ 1097d874-4468-4e76-bdbc-c893a5dbfdc0
using Plots, PlutoUI, DifferentialEquations

# ╔═╡ 18a83e00-0dd4-4fe7-a5be-644361f875d3
TableOfContents()

# ╔═╡ 4ad5df20-f85d-11ec-3803-a7ed7ed8f2f3
md"""
# Dynamical Systems. Flows in 1D

A **flow** is a dynamical system that evolves continuosly over time. 

A **trajectory** follows a continuous path in the space of states.

For one-dimensional flows the trajectory is given by the evolution of the single variable as a function of time: $x(t)$

The deterministic rule for the time evolution is given by the **instantaneous rate of change in time** of the variable, expressed as the time derivative: $\frac{dx}{dt}$, also abbreviated as $\dot{x}$.

Therefore as the evolution rule depends only on the state of the system $x$, a one-dimensional flow is fully characterized by the equation:

$\dot{x}=f(x)$

The function $f$ is known as the **vector field** of the flow (the reason for this name will become clear later).
"""

# ╔═╡ 157c155f-241e-4d39-93f7-35997d83191d
md"""
# Simplest example: Continuous population with continuous time

As in the continuous population discrete time example, the variable $x$ corresponds to the size of the population, expressed as a continuous value after some normalization, but now the time evolution is also continuous. Therefore, the growth rate is given by its instantaneous value $r$, and the equation for the flow is written as:

$\dot{x}=rx$

In the animated graphic below we can set the initial value of the population $x(0)$ and the instantaneous rate of growth $r$ with the sliders. The thirds slider control the refresh rate of the graphics.
"""

# ╔═╡ 9490e7f8-668e-4acc-8efd-dcc4d355acb9
md"""
$(@bind ticks Clock(1))
$(@bind reset Button("Reset"))
"""

# ╔═╡ d9bb5aa8-2708-40c3-8d42-3b21e6d634c0
md"""
x($0$) : $(@bind x0 Slider(0.0:0.02:1.0,default=0.1;show_value=true)) \
r : $(@bind R Slider(-0.1:0.001:0.1,default=-0.001;show_value=true)) \
t\_refresh : $(@bind tail Slider(0.001:0.001:1.0,default=0.1;show_value=true)) \
"""

# ╔═╡ dd1af6a9-c048-4845-8858-8051d23f6bee
md"""
The upper graphics displays snapshots (with time intervals *t_refresh*) of the states of the dynamical system given by the variable $x$ as red dots in the state space (black line). 

For small values of *t_refresh* the dots draw a continuos red line that corresponds to the **trajectory** of the dynamical system. 

The lower graph displays the value of the variable $x$ as a function of time $t$ as a continuous blue line. This is the **solution** of the system.

Starting from $x(0)$ at time $t=0$, the population evolves in time either growing without bounds (for $r>0$) or decreasing toward zero (for $r<0$). A marginal case occurs when $r=0$ and the population remains fixed at the initial value.

"""

# ╔═╡ 75d2969a-b83c-406d-90d5-42e4236042e4
function linear(u,p,t)
	p[1]*u[1]
end;	

# ╔═╡ 6a97c988-3ec8-4271-b61c-47c40b1deee6
md"""
# Logistic Equation (Flows)

Analogously to what we did with the maps, we can limit the population growth by considering an instantaneous growth rate that takes a maximum positive value ($r=r_0$) for zero population, becomes zero ($r=0$) for an optimal population size ($x=1$) and becomes negative for larger values of $x$:

$r = r_0(1-x)$

and the dynamical system (flow) is then defined then by the logistic equation with continuous time:

$\dot{x}=r_0x(1-x)$
"""

# ╔═╡ 50153c10-d204-4e8b-8bcf-007e22f287ee
function logistic(u,p,t)
	p[1]*u[1]*(1.0-u[1])
end;	

# ╔═╡ 7e7239c2-2edf-4005-94be-90c85b282b01
md"""
$(@bind ticks2 Clock(1))
$(@bind reset2 Button("Reset"))
"""

# ╔═╡ 1115e7e3-06f2-4f51-a71d-85ada76ce3a9
md"""
x($0$) : $(@bind x02 Slider(0.0:0.02:2.0,default=0.1;show_value=true)) \
r₀ : $(@bind r0 Slider(0.01:0.01:4.0,default=1.0;show_value=true)) \
t\_refresh : $(@bind tail2 Slider(0.001:0.001:1.0,default=0.1;show_value=true)) \
"""

# ╔═╡ 5a3db641-7fce-46a6-8953-59f50081556c
md"""
The effect of the decreasing growth rate is to stabilize the population at some intermediate value $x_*=1$. This is a **fixed point** of the flow and also an **attractor**.

All initial condition (therefore all trajectories) converge to this attractor and 
the only effect of $r$ is only on the rate of this convergence.

For one dimensional flows (in contrast to maps) fixed points are the only possible limit set. 
"""

# ╔═╡ 54d209de-d485-44f8-911d-d988792b1057
md"""
# Newton's Cooling Law

In this case the variable is the temperature of the object $T$ and $k$ is its heat exchange rate (which depends on the specific heat and the contact area). $T_{eq}$ is the ambient temperature. The temperature evolution is given by:

$\dot{T}=-k(T-T_{eq})$
"""

# ╔═╡ 8fef25bf-0633-4291-b27d-c3b6414147da
function cooling(u,p,t)
	-p[1]*(u[1]-p[2])
end;	

# ╔═╡ 4ff8085d-cc86-4cdc-bbc1-3d2c5ec52614
md"""
$(@bind ticks3 Clock(1))
$(@bind reset3 Button("Reset"))
"""

# ╔═╡ fff09537-0806-4179-8620-cc022b04f40a
md"""
T($0$) : $(@bind T0 Slider(0:1:100,default=80;show_value=true)) \
k : $(@bind k Slider(0.01:0.01:4.0,default=1.0;show_value=true)) \
Tₑ : $(@bind Teq Slider(0:1:100,default=20;show_value=true)) \
t\_refresh : $(@bind tail3 Slider(0.001:0.001:1.0,default=0.1;show_value=true)) \
"""

# ╔═╡ 36a0fb9c-1821-47d4-8dae-9a93b232ee1d
md"""
As with the logistic growth, all initial conditions converge to a single fixed point at $T_*=T_{eq}$. The parameter $k$ controls the rate at witch the object converge (either heating up or cooling down) to the ambient temperature.
"""

# ╔═╡ 3044f1e4-96ae-4482-8e05-1e5d7f7b2387
md"""
## Many Initial Conditions (Ensemble of Trajectories)

The set (or ensemble) of all possible trajectories in the space of states is called the **phase portrait** of the system (the origin of the name is also historical).

For the case of the cooling law or the logistic growth (if we limit the space of possible states to $R\geq 0$) all trajectories overlap to only two possible trajectories one converging to the fixed point from above and other from below.
"""

# ╔═╡ efc6ca91-1452-41b8-9270-30960e2091d8
begin
	function prob_func(prob,i,repeat)
  		remake(prob,u0=u0_arr[i])
	end
	prob = ODEProblem(cooling, [0.0], (0,10.0), [k,Teq])
	u0_arr = -10:10:100
	ensemble_prob = EnsembleProblem(prob,prob_func=prob_func)
	sole = solve(ensemble_prob,EnsembleThreads(),trajectories=length(u0_arr))
end;

# ╔═╡ 7891d31e-9873-4f60-9646-eb9c90ca45c1
begin
	p1a = plot(legend=false,ylims=(-0.1,0.1),yaxis=false,yticks=false)
	for s in sole.u
		plot!(p1a,s.u,zero(s.u),alpha=0.2,arrows=true,)
	end	
	p2a = plot(sole,vars=(0,1),xlabel="t",ylabel="T")
	plot(p1a,p2a,size=(900,400),layout=grid(2,1,heights=[0.2 ,0.8]))
end	

# ╔═╡ 3a6821d7-392a-4dd4-b5c8-c247edf9132f
mutable struct Flow
	f
	x0 
	t
	xv 
end	

# ╔═╡ 9b3a5d92-6f53-486d-aa14-d648c7330869
begin
	reset
	f1=Flow(linear,x0,0.0,[])
	p0b = plot(legend=false);
end;

# ╔═╡ f480e291-b67c-4a2f-9f44-914a340df81e
html"""
<style>
input[type*="range"] {
	width: 60%;
}
</style>
"""

# ╔═╡ e1e02db0-de4b-4349-af8e-acac7bf046f5
function plot_state_flow(x,xrange,t)
	title = "t=$(round(x,digits=4))  x(t) = $(round(x,digits=8))"
	p1 = plot(xlims=xrange,ylims=(-0.1,0.1),size=(800,100),yaxis=false,yticks=false,title=title,legend=false)
	plot!(p1,[xrange[1],xrange[2]],[0,0],c=:black)
	#scatter!(p1,[x],[0],c=:red)
end;

# ╔═╡ a204c1f9-8fd9-4feb-971d-4cf008730641
p0 = plot_state_flow(f1.x0,[0,2.0],0.0);

# ╔═╡ 5dfedd16-77be-4350-9e75-afae14684fe2
begin
	reset2
	f2=Flow(logistic,x02,0.0,[])
	p1b = plot(legend=false)
	p1 = plot_state_flow(f2.x0,[0,2.0],0.0)
end;

# ╔═╡ 1d6b10c9-76f0-4ddd-8d50-893cc3618aae
begin
	reset3
	f3=Flow(cooling,T0,0.0,[])
	p2b = plot(legend=false);
	p2 = plot_state_flow(f3.x0,[0,100.0],0.0);
end;

# ╔═╡ f1d2afd0-e6b4-4272-91f7-1b225759534f
function plot_state_flow!(p1,x,t;alpha=1.0)
	title = "t=$(round(x,digits=4))   x(t) = $(round(x,digits=8))"
	scatter!(p1,[x],[0],c=:red,alpha=alpha)
	title!(p1,title)
end;


# ╔═╡ 936816f9-5a79-474c-8847-96a1f0cb9a6a
begin
	ticks
	sol = solve(ODEProblem(f1.f,f1.x0,(f1.t,f1.t+tail),[R]))
	f1.x0 = sol.u[end]
	f1.t = sol.t[end]
	plot_state_flow!(p0,f1.x0,f1.t;alpha=0.3)
	plot!(p0b,sol.t,getindex.(sol.u,1),c=:blue)
	plot(p0,p0b,layout=grid(2,1,heights=[0.2 ,0.8]))
end	

# ╔═╡ 807c6201-c455-473d-96e1-9d99dfaa7c55
begin
	ticks2
	sol2 = solve(ODEProblem(f2.f,f2.x0,(f2.t,f2.t+tail2),[r0]))
	f2.x0 = sol2.u[end]
	f2.t = sol2.t[end]
	plot_state_flow!(p1,f2.x0,f2.t;alpha=0.3)
	plot!(p1b,sol2.t,getindex.(sol2.u,1),c=:blue)
	plot(p1,p1b,layout=grid(2,1,heights=[0.2 ,0.8]))
end	

# ╔═╡ c8cd81a0-256a-4d58-90a6-6d09842f7224
begin
	ticks3
	sol3 = solve(ODEProblem(f3.f,f3.x0,(f3.t,f3.t+tail3),[k,Teq]))
	f3.x0 = sol3.u[end]
	f3.t = sol3.t[end]
	plot_state_flow!(p2,f3.x0,f3.t;alpha=0.3)
	plot!(p2b,sol3.t,getindex.(sol3.u,1),c=:blue)
	plot(p2,p2b,layout=grid(2,1,heights=[0.2 ,0.8]))
end	

# ╔═╡ Cell order:
# ╠═8992e664-fbdd-4b25-9a24-2c72a34d9e10
# ╠═1097d874-4468-4e76-bdbc-c893a5dbfdc0
# ╟─18a83e00-0dd4-4fe7-a5be-644361f875d3
# ╟─4ad5df20-f85d-11ec-3803-a7ed7ed8f2f3
# ╟─157c155f-241e-4d39-93f7-35997d83191d
# ╟─9490e7f8-668e-4acc-8efd-dcc4d355acb9
# ╟─d9bb5aa8-2708-40c3-8d42-3b21e6d634c0
# ╟─936816f9-5a79-474c-8847-96a1f0cb9a6a
# ╟─dd1af6a9-c048-4845-8858-8051d23f6bee
# ╟─75d2969a-b83c-406d-90d5-42e4236042e4
# ╟─a204c1f9-8fd9-4feb-971d-4cf008730641
# ╟─9b3a5d92-6f53-486d-aa14-d648c7330869
# ╟─6a97c988-3ec8-4271-b61c-47c40b1deee6
# ╠═50153c10-d204-4e8b-8bcf-007e22f287ee
# ╟─7e7239c2-2edf-4005-94be-90c85b282b01
# ╟─1115e7e3-06f2-4f51-a71d-85ada76ce3a9
# ╠═807c6201-c455-473d-96e1-9d99dfaa7c55
# ╟─5a3db641-7fce-46a6-8953-59f50081556c
# ╟─5dfedd16-77be-4350-9e75-afae14684fe2
# ╟─54d209de-d485-44f8-911d-d988792b1057
# ╠═8fef25bf-0633-4291-b27d-c3b6414147da
# ╟─4ff8085d-cc86-4cdc-bbc1-3d2c5ec52614
# ╟─fff09537-0806-4179-8620-cc022b04f40a
# ╟─c8cd81a0-256a-4d58-90a6-6d09842f7224
# ╟─36a0fb9c-1821-47d4-8dae-9a93b232ee1d
# ╟─3044f1e4-96ae-4482-8e05-1e5d7f7b2387
# ╟─7891d31e-9873-4f60-9646-eb9c90ca45c1
# ╟─efc6ca91-1452-41b8-9270-30960e2091d8
# ╟─1d6b10c9-76f0-4ddd-8d50-893cc3618aae
# ╟─3a6821d7-392a-4dd4-b5c8-c247edf9132f
# ╟─f480e291-b67c-4a2f-9f44-914a340df81e
# ╠═e1e02db0-de4b-4349-af8e-acac7bf046f5
# ╟─f1d2afd0-e6b4-4272-91f7-1b225759534f
