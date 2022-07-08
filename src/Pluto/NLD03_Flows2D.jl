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

# ╔═╡ 48c61fc3-94fb-44d4-8a64-1697c0ff0fbd
using Pkg;Pkg.add("Plots");Pkg.add("PlutoUI");Pkg.add("DifferentialEquations");Pkg.add("ForwardDiff"); Pkg.add("StaticArrays"); Pkg.add("IntervalRootFinding")

# ╔═╡ 1097d874-4468-4e76-bdbc-c893a5dbfdc0
using Plots, PlutoUI, DifferentialEquations, ForwardDiff, StaticArrays, IntervalRootFinding

# ╔═╡ d511c3c8-c596-48a5-8182-4dbcaa607eb6
include("../NLD_utils.jl")

# ╔═╡ 18a83e00-0dd4-4fe7-a5be-644361f875d3
TableOfContents()

# ╔═╡ 4ad5df20-f85d-11ec-3803-a7ed7ed8f2f3
md"""
# Dynamical Systems. Flows in 2D

The **space of states** is now bi-dimensional $(x,y)$. Each point of this space corresponds to a unique state of the system. 

A two-dimensional flow is defined by assigning to each state $(x,y)$ an instantaneous rate of change in time **in each of the directions**, expressed by its time derivatives $(\dot{x},\dot{y})$. And, as the deterministic rule depend only on the state of the system, the flow is fully characterized by the pair of equations:

$\dot{x}=f(x,y)$

$\dot{y}=g(x,y)$

As with 1D flows where we could interpret the time derivative as the velocity of the trajectory along the line, we can now associate the pair of time derivatives $(\dot{x},\dot{y})$ as the horizontal and vertical components of the velocity vector of the trajectory in the plane.  

At each point $(x,y)$ this velocity vector points in the direction of the time evolution of the trajectory. The pair of functions ($f,g$) in right side are the **vector field** of the dynamical system.
"""

# ╔═╡ 157c155f-241e-4d39-93f7-35997d83191d
md"""
# The harmonic oscillator

We will take the mass hanging in a spring as the physical model for the HO. The dynamical state of the system is fully characterized by its position with respect to the resting place $x$, and its velocity $v$.

Assuming that the friction of system is very small and that it can be neglected, the vector field can be written as:

$\dot{x}=v$

$\dot{v}=-kx$

The first equation simply expresses the relationship between the two quantities. In fact, the definition of velocity is the time derivative of the position. 

The second equation expresses that the variation of the velocity, also known as acceleration, is proportional to the displacement and with opposite sign. 

What does this mean? That the more the mass is displaced from its equilibrium position, the greater will be the tendency for it to return to equilibrium.  

The physical cause of this is the elastic force of the spring (in fact this equation is none other than the famous second law of Newton). The negative sign is what generates the oscillation. 

"""

# ╔═╡ a04beba2-07be-4113-a578-048c12a34d82
function oscharm!(du,u,p,t)
	du[1]=u[2]
	du[2]=-p[1]*u[1]
	du
end

# ╔═╡ ecbab468-a04e-45df-91a8-4e76efd92163
md"""
In the graph below we display snapshots of the states of the system at different times, along the continuous trajectory in ($x,v$), starting from the initial condition ($x(0),y(0)$).
"""

# ╔═╡ d16430c4-6c92-4d09-b36d-73fd4164123f
md"""
$(@bind ticks Clock(1))
$(@bind reset Button("Reset"))
"""

# ╔═╡ d9bb5aa8-2708-40c3-8d42-3b21e6d634c0
md"""
x(0) $(@bind x0 Slider(-1.0:0.02:1.0,default=0.7;show_value=true)) 
v(0) : $(@bind y0 Slider(-1.0:0.01:1.0,default=0.0;show_value=true)) \
k : $(@bind k Slider(0.1:0.01:2.0,default=0.7;show_value=true)) 
t\_refresh : $(@bind tail Slider(0.001:0.001:0.1,default=0.02;show_value=true))
"""

# ╔═╡ 7cf3899a-a21c-4804-ba12-6543b1916cb1
md"""
## Vector field 

We can also display the **vector field** ($f,g$) as vectors (arrows) attached to a grid of points in the state space ($x,v$). At each point, vectors are headed toward the future direction of the trajectory passing through that point, and are tangent to this curve. 

In the graph below a scaled version of the vector field is displayed with blue arrows, along with a single trajectory in black, starting from ($x(0),v(0)$). Also, the time evolution of both variables is shown in the two graphs below.
"""

# ╔═╡ b5dea9d8-ebae-4a8a-aaec-180f8c83fdfd
flow2d_vectorfield(oscharm!,[x0,y0],30.0,[k];title="Harmonic Oscillator")

# ╔═╡ 3f64c1ef-4bcb-49ed-a2a9-bc041950d514
begin
	sol = solve(ODEProblem(oscharm!, [x0,y0], (0,30.0), [k]))
	pa = plot(sol,vars=(0,1),legend=false,xlabel="t",ylabel="x")
	pb = plot(sol,vars=(0,2),legend=false,xlabel="t",ylabel="v")
	plot!(pa,[0,30],[0,0],c=:black)
	plot!(pb,[0,30],[0,0],c=:black)
	plot(pa,pb,layout=(2,1),size=(900,400))
end	

# ╔═╡ 9ca2ec27-8600-4d6f-8c9f-05e08ed137f3
md"""
# Damped harmonic oscillator

Now we will take into account the dissipation in the harmonic oscillator. 

In its simplest form the dissipation, or friction, or damping (slightly different flavors of the same phenomenon) can be put as a force that always acts against the movement. If the velocity is zero there are no friction forces and, as the velocity increases, the friction grows linearly with the speed and always oposite in direction, so we write the force as $-\gamma v$ where $\gamma$ is the damping coefficient. 

Then incorporating this term in the second equation (which is the movement equation realting the forces to the changes in the velocity) the damped harmonic oscillator can be written as:

$\dot{x} = v$

$\dot{v} = -\gamma v - k x$

This is a **linear** system because both $f(x,v)$ and $g(x,v)$ are linear functions of the variables, i.e. do not involve powers higher than one, such as $x^2$, $xv$, $v^3$, etc.
"""	

# ╔═╡ 8fef25bf-0633-4291-b27d-c3b6414147da
function oscharmdamp!(du,u,p,t)
	du[1]=u[2]
	du[2]=-p[1]*u[1]-p[2]*u[2]
	du
end

# ╔═╡ bf009c87-d4cd-4001-8bd2-c95e3581fc8d
md"""
x(0) $(@bind x02 Slider(-0.5:0.02:0.5,default=0.0;show_value=true)) 
v(0) : $(@bind y02 Slider(-0.5:0.01:0.5,default=0.1;show_value=true)) \
k : $(@bind k2 Slider(0.1:0.01:2.0,default=0.1;show_value=true)) 
γ : $(@bind γ Slider(0.0:0.001:0.2,default=0.1;show_value=true))
"""

# ╔═╡ 5ee8744a-7b12-48e1-8c56-d641aed7e0c5
flow2d_vectorfield(oscharmdamp!,[x02,y02],250.0,[k2,γ];title="Damped Harmonic Oscillator")

# ╔═╡ a7440fa6-8ae6-425c-93e0-7c3dfa53416e
begin
	sol2 = solve(ODEProblem(oscharmdamp!, [x02,y02], (0,50.0), [k2,γ]))
	pa2 = plot(sol2,vars=(0,1),legend=false,xlabel="t",ylabel="x")
	pb2 = plot(sol2,vars=(0,2),legend=false,xlabel="t",ylabel="v")
	plot!(pa2,[0,50.0],[0,0],c=:black)
	plot!(pb2,[0,50.0],[0,0],c=:black)
	plot(pa2,pb2,layout=(2,1),size=(900,400))
end	

# ╔═╡ e7e3910b-bd29-4bf9-b416-65dd5cfc2502
md"""
## Stable Fixed Point (Attractor)
"""

# ╔═╡ e0b5fcfb-0f0d-41f6-b5f8-280c2bf2f44c
md"""
For $\gamma>0$ all initial conditions follow a spiraling trajectory that evolves towards the central point: 

($x=0,v=0$). 

This is a **fixed point** of the flow because for this point the vector field is null ($f=0,g=0$). And therefore there is no time evolution, or the fixed point evolves into itself. It is also called an **equilibrium point**. 

This fixed point is **stable** because if we perturb the system from the equilibrium position it returns to it. In fact, all trajectories converge to the fixed point and this makes it an **atractor**.

For $\gamma=0$ there is no more friction and we return to the classic harmonic oscillator. The point ($x=0,v=0$) is still a fixed point but no longer an attractor. 
"""

# ╔═╡ 7c0d16f6-242d-4ae7-b5f3-444dd227bd3b
md"""

# Self-oscillator: Simple reed model (Rayleigh)


In the "Theory of Sound" (1877) Rayleigh proposed a simple model for the blowing of a clarinet reed, in which the dissipation coefficient ($\gamma$) is negative (i.e. acts in favor of the motion and delivers energy) for small velocities and positive (i.e. acts by slowing down the motion) for large velocities. 

$\gamma = v^2 - \mu$ 

This is a very general model of a self-oscillator that can be written by replacing the expression of the nonlinear dissipation in the harmonic oscillator.

$\dot{x} = v$

$\dot{v} = -(v^2-\mu)v-kx$

This is a **nonlinear system** because it includes a $v^3$ term
"""

# ╔═╡ 2cd075b8-0fd1-481b-b230-1bf3a51f2d5f
function reed!(du,u,p,t)
    (μ,K) = p
    du[1] = u[2]
	du[2] = u[2]*(μ-u[2]*u[2])-K*u[1]          
end

# ╔═╡ 7557f8d7-2654-41c2-a443-9d5a76622501
md"""
x(0) $(@bind x03 Slider(-1.0:0.01:1.0,default=0.1;show_value=true)) 
K : $(@bind K2 Slider(1.0:0.1:5.0,default=2.0;show_value=true)) \
μ : $(@bind μ Slider(-0.1:0.001:0.3,default=0.1;show_value=true))
tmax : $(@bind tmax2 Slider(5.0:5.0:40.0,default=10.0;show_value=true))
"""

# ╔═╡ cb6d88cd-d6fb-4730-95a6-35a59224f35c
flow2d_vectorfield(reed!,[x03,0],tmax2,[μ,K2];title="Simple Reed Model")

# ╔═╡ c45e6d60-3086-42f5-ba1b-df6de93827ed
begin
	sol4 = solve(ODEProblem(reed!, [x03,0], (0,tmax2), [μ,K2]))
	pa4 = plot(sol4,vars=(0,1),legend=false,xlabel="t",ylabel="x")
	pb4 = plot(sol4,vars=(0,2),legend=false,xlabel="t",ylabel="v")
	plot!(pa4,[0,50.0],[0,0],c=:black)
	plot!(pb4,[0,50.0],[0,0],c=:black)
	plot(pa4,pb4,layout=(2,1),size=(900,400))
end	

# ╔═╡ 2c1199ef-6e95-43b8-92e8-690b3ea1d45c
md"""

## Unstable Fixed Point (Repulsor)
For $\mu<0$ we have a behavior similar to the case of the damped harmonic oscillator: the origin ($0,0$) is a fixed point (the vector field is null) and is an attractor. 

However, for $\mu>0$ there is a **qualitative change in the dynamics**. the origin is still a fixed point but now it is unstable (try placing a initial condition close to it). The trajectories spiral away from the center. 

This is the effect of negative dissipation and corresponds, in the physical system, to the energy delivery by the blowing pressure on the reed. An unstable fixed point from which the trajectories move away no matter how close they start from it is known as a **repulsor**

These changes of stability that happen in nonlinear sytems are known as **bifurcations**

## Stable Limit Cycle

The trajectory that spirals away from the repulsor for $\mu>0$ corresponds, in the physical system, to the reed that is set in vibration by the blowing energy, oscillating with increasing amplitude until that energy is compensated by the nonlinear dissipation (which is the form proposed by Rayleigh). Then the aplitude stops growing and stabilizes to an oscillation of well-defined amplitude. 

This corresponds in the state space to a periodic orbit to which all trajectories tend. The initial conditions that start with higher amplitude also decay to the same limit set. 

This new limit set is known as a **limit cycle**. And as all near initial conditions converge to it, it is a *stable limit cycle*. 

Note that, unlike the harmonic oscillator where each initial conditions gave rise to an oscillation with different amplitude, in this system all initial conditions converge to an equilibrium oscillation with a definite amplitude and frequency.


"""

# ╔═╡ 877f1288-8c71-4a96-8cec-b2f5e9c4c688
md"""
# Phase portrait (multiple trajectories)

The set (or ensemble) of all possible trajectories in the space of states is called the **phase portrait** of the system.

We can evolve a grid of initial conditions and observe the "flow" of the dynamical system in the state space. This representation is complementary to that of the vector field (in the context of calculus, they are called integral and differential respectively).
"""

# ╔═╡ 940e5c21-2d93-4fa1-a8b2-c8b1eb20dd45
md"""
k : $(@bind k3 Slider(0.3:0.01:2.0,default=0.5;show_value=true)) 
γ : $(@bind γ3 Slider(0.0:0.1:2.0,default=0.1;show_value=true))
"""

# ╔═╡ 48a9bf5d-5c44-4c4e-b9e4-a1ee53a4a978
begin
	function prob_func(prob,i,repeat)
  		remake(prob,u0=u0_arr[i])
	end
	prob = ODEProblem(oscharmdamp!, [0.0,0.0], (0,20.0), [k3,γ3])
	u0_arr = vec([[0.2*i-1.0,0.2*j-1.0] for i=0:10, j=0:10])
	ensemble_prob = EnsembleProblem(prob,prob_func=prob_func)
	sole1 = solve(ensemble_prob,EnsembleThreads(),trajectories=length(u0_arr))
end;

# ╔═╡ 13b4722a-ab98-4d09-9d0f-998a0596887a
begin
	plot(sole1,vars=(1,2),xlabel="x",ylabel="y",c=:blue,alpha=0.5)
	plot!([-1.0,1.0],[0,0],c=:black,title="Damped Harmonic Oscillator")
	plot!([0,0],[-1.0,1.0],c=:black,size = (800,800))
	xlims!((-1,1))
	ylims!((-1,1))
end	

# ╔═╡ 50929479-b635-4cd3-9491-038c555d3e60
md"""
K : $(@bind K4 Slider(1.0:0.1:5.0,default=2.0;show_value=true)) 
μ : $(@bind μ4 Slider(-0.5:0.1:0.5,default=0.5;show_value=true))
"""

# ╔═╡ 7de645ef-ab42-47cd-b8be-f7533c94e066
begin
	prob2 = ODEProblem(reed!, [0,0], (0,20.0), [μ4,K4])
	ensemble_prob2 = EnsembleProblem(prob2,prob_func=prob_func)
	sole2 = solve(ensemble_prob2,EnsembleThreads(),trajectories=length(u0_arr))
end;

# ╔═╡ 0cfeb6ad-605c-4e24-85ae-45be1339a79a
begin
	plot(sole2,vars=(1,2),xlabel="x",ylabel="y",c=:blue,alpha=0.3)
	plot!([-1.0,1.0],[0,0],c=:black,title="Simple Reed Model")
	plot!([0,0],[-1.0,1.0],c=:black,size = (800,800))
	xlims!((-1,1))
	ylims!((-1,1))
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

# ╔═╡ 842dff44-ed12-42fd-b0e4-e4d33025f8ae
mutable struct Flow
	f::Function
	u0::Array 
	t::Number
	xv::Array 
end	

# ╔═╡ c26b6f12-fd3d-4a96-9ca0-ddcf267e5a3c
begin
	reset
	f1=Flow(oscharm!,[x0,y0],0.0,[])
	p0 = plot_state_flow2D([x0,y0],[-1.0,1.0],[-1.0,1.0],0.0);
	p0b = plot(legend=false);
end;

# ╔═╡ 37c6bdde-e965-4d09-8a6d-3beb3366cbfb
begin
	ticks
	sol0 = solve(ODEProblem(f1.f,f1.u0,(f1.t,f1.t+tail),[k]))
	f1.u0 = sol0.u[end]
	f1.t = sol0.t[end]
	plot_state_flow2D!(p0,f1.u0,f1.t;alpha=0.25)
	#plot!(p0b,sol0.t,getindex.(sol0.u,1),c=:blue)
	#plot!(p0b,[0,sol0.t[end]],[0,0],c=:black)
	plot(p0,size=(600,600))
end	

# ╔═╡ Cell order:
# ╠═48c61fc3-94fb-44d4-8a64-1697c0ff0fbd
# ╠═1097d874-4468-4e76-bdbc-c893a5dbfdc0
# ╠═d511c3c8-c596-48a5-8182-4dbcaa607eb6
# ╟─18a83e00-0dd4-4fe7-a5be-644361f875d3
# ╟─4ad5df20-f85d-11ec-3803-a7ed7ed8f2f3
# ╟─157c155f-241e-4d39-93f7-35997d83191d
# ╠═a04beba2-07be-4113-a578-048c12a34d82
# ╟─ecbab468-a04e-45df-91a8-4e76efd92163
# ╟─d16430c4-6c92-4d09-b36d-73fd4164123f
# ╟─d9bb5aa8-2708-40c3-8d42-3b21e6d634c0
# ╟─37c6bdde-e965-4d09-8a6d-3beb3366cbfb
# ╟─7cf3899a-a21c-4804-ba12-6543b1916cb1
# ╟─b5dea9d8-ebae-4a8a-aaec-180f8c83fdfd
# ╟─3f64c1ef-4bcb-49ed-a2a9-bc041950d514
# ╟─c26b6f12-fd3d-4a96-9ca0-ddcf267e5a3c
# ╟─9ca2ec27-8600-4d6f-8c9f-05e08ed137f3
# ╠═8fef25bf-0633-4291-b27d-c3b6414147da
# ╟─bf009c87-d4cd-4001-8bd2-c95e3581fc8d
# ╟─5ee8744a-7b12-48e1-8c56-d641aed7e0c5
# ╟─a7440fa6-8ae6-425c-93e0-7c3dfa53416e
# ╟─e7e3910b-bd29-4bf9-b416-65dd5cfc2502
# ╟─e0b5fcfb-0f0d-41f6-b5f8-280c2bf2f44c
# ╟─7c0d16f6-242d-4ae7-b5f3-444dd227bd3b
# ╠═2cd075b8-0fd1-481b-b230-1bf3a51f2d5f
# ╟─7557f8d7-2654-41c2-a443-9d5a76622501
# ╟─cb6d88cd-d6fb-4730-95a6-35a59224f35c
# ╟─c45e6d60-3086-42f5-ba1b-df6de93827ed
# ╟─2c1199ef-6e95-43b8-92e8-690b3ea1d45c
# ╟─877f1288-8c71-4a96-8cec-b2f5e9c4c688
# ╟─940e5c21-2d93-4fa1-a8b2-c8b1eb20dd45
# ╟─13b4722a-ab98-4d09-9d0f-998a0596887a
# ╟─48a9bf5d-5c44-4c4e-b9e4-a1ee53a4a978
# ╟─50929479-b635-4cd3-9491-038c555d3e60
# ╟─0cfeb6ad-605c-4e24-85ae-45be1339a79a
# ╟─7de645ef-ab42-47cd-b8be-f7533c94e066
# ╟─f480e291-b67c-4a2f-9f44-914a340df81e
# ╟─5bbb333e-47a3-464b-82f5-cb8053c5c733
# ╟─64022ca2-6323-4428-b938-843e228d4922
# ╟─842dff44-ed12-42fd-b0e4-e4d33025f8ae
