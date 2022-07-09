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

# ╔═╡ 1878f6d9-d9a7-4bf5-9525-2b021ebd1521
using Pkg;Pkg.add("Plots");Pkg.add("PlutoUI");Pkg.add("DifferentialEquations");Pkg.add("ForwardDiff"); Pkg.add("StaticArrays"); Pkg.add("IntervalRootFinding")

# ╔═╡ 794bf502-f9e9-11ec-03c5-ad174a6ccd59
using Plots, PlutoUI, DifferentialEquations, ForwardDiff, StaticArrays, IntervalRootFinding

# ╔═╡ a830882a-8a52-4171-a0da-0449d4f5c901
include("../NLD_utils.jl")

# ╔═╡ 69c47e84-0bb9-425a-8cc8-3c50e2ffa2b3
TableOfContents()

# ╔═╡ 9969e088-322f-455b-8f3c-edb0ba8faed1
md"""
# Bouncing Ball
This system is simpler than the previous one because the acceleration is constant and negative. The system for a ball falling by the action of gravity can be written simply as follows

$\dot{h} = v$

$\dot{v} = -g -\gamma v$

where $h$ is the heigth of the ball from the ground, $v$ its vertical velocity and $g$ is the gravity acceleration. We also add dissipation, a drag force or air resistance that is always present and acting against the velocity with controled by the parameter $\gamma$ (we can turn it off later). 

However, we have to add a condition of bouncing on the ground. Let us consider two cases: Elastic rebound (Ball Soft) and rigid rebound (Ball Hard).

For the elastic rebound case we assume that when the height reaches zero it collides with a spring of elastic constant $k$ that makes it go back up (in this case the heigth of the ball can take negative values because of the deformation), this implies that we have to evaluate the sign of $h$ and apply the vector field written above only for $h>0$ and add a term with the elastic force for $h\leq 0$.

### Ball Soft


$\dot{h} = v$

$\dot{v} = -g  -\gamma v\quad [if \quad h>0]$

$\dot{v} = -g - kh -\gamma v \quad [if \quad h\leq 0]$

"""

# ╔═╡ 6527dc0e-6c30-4d84-8f41-5ec3e42eaf86
function ballsoft!(du,u,p,t)
    (g,K,γ) = p
    du[1] = u[2]
    if (u[1]>0)
        du[2] = -g - γ*u[2]
    else
        du[2] = -g-K*u[1] - γ*u[2]
    end      
end

# ╔═╡ 65bafd81-f8ff-49d4-a882-3007c8c8a71a
md"""
h(0) $(@bind h0 Slider(0.5:1:50,default=30.0;show_value=true)) 
K : $(@bind K Slider(10.0:5.0:50.0,default=20.0;show_value=true)) \
γ : $(@bind γ2 Slider(0.0:0.001:0.2,default=0.1;show_value=true))
tmax : $(@bind tmax1 Slider(0.5:0.5:40.0,default=10.0;show_value=true))
"""

# ╔═╡ 507ddd8e-9609-40d0-bd4f-682f72475f93
flow2d_vectorfield(ballsoft!,[h0,0],tmax1,[9.8,K,γ2];title="Ball Soft",xlims=[-10,50.1],ylims=[-30,30])

# ╔═╡ 6e3ebd5c-b797-44f9-82ba-8ad9a18152eb
begin
	sol3 = solve(ODEProblem(ballsoft!, [h0,0], (0,tmax1), [9.8,K,γ2]))
	pa3 = plot(sol3,vars=(0,1),legend=false,xlabel="t",ylabel="x")
	pb3 = plot(sol3,vars=(0,2),legend=false,xlabel="t",ylabel="v")
	plot!(pa3,[0,tmax1],[0,0],c=:black)
	plot!(pb3,[0,tmax1],[0,0],c=:black)
	plot(pa3,pb3,layout=(2,1),size=(900,400))
end	

# ╔═╡ 8826d6ea-a2cb-430a-9e03-5fd459df24b6
md"""
# Bouncing Ball in a Bowl

This is one of the simplest physical systems that exhibits chaos:

It is a ball bouncing inside a circular cavity. Unlike the original bouncing ball that occurred in one dimension, this physical system is two dimensional, so we need four variables to describe it: the position ($x,y$) , and the velocity components ($v_x,v_y$). The equations are the same as the original bouncing ball without dissipation in the coordinate and in fact the system seems overly simple:

$\dot{x}=v_x$

$\dot{y}=v_y$

$\dot{v_x}=0$

$\dot{v_y}=-g$


The main complication here comes from the boundary conditions, we will consider a "hard" collision, which if we assume the origin of coordinates in the center of the cavity and the radius of the same equal to will be when the condition is met:

$x^2+y^2=1$

In that case we suppose a perfectly elastic collision first and the effect of the reflection is to **invert the normal component of the velocity** and to **conserve the tangential one**.
Written in cartesian components, this seems a bit complicated but it is not difficult to compute.

$v_x \rightarrow v_x(y^2-x^2)-2xyv_y$

$v_y \rightarrow v_y(x^2-y^2)-2xyv_x$


We write the equations of the system, the condition for the collision and the function bounce! to change the direction of the velocity using the formulae above.
"""

# ╔═╡ 85011323-6e1b-4f28-8e94-9c9d2e3ae686
function ballcirclehard!(du,u,p,t)
  du[1] = u[3]
  du[2] = u[4]
  du[3] = 0  
  du[4] = -p[1]
  du  
end

# ╔═╡ 295c9623-7657-4ec4-ab35-af8d068adacc
function collision(u,t,integrator)
  1.0-sqrt(u[1]*u[1]+u[2]*u[2])
end

# ╔═╡ 55a5f7f0-e0fe-41dd-891a-9d733862dc9f
function bounce!(integrator)
  (x,y,vx,vy) = integrator.u
  integrator.u[3] = vx*(y*y-x*x)-2*x*y*vy
  integrator.u[4] = vy*(x*x-y*y)-2*x*y*vx
end

# ╔═╡ d18add15-daa9-4a46-9a1a-9304e74459fe
md"""
δ : $(@bind δ Slider(-8:0.5:-2,default=-3;show_value=true)) 
tmax : $(@bind tmax Slider(0.1:0.1:40,default=10;show_value=true)) \
"""

# ╔═╡ fa92921c-0d13-4e2f-ae97-37af89acb500
begin
	u1 = [0.5,0.85,0.0,0.0]
	u2 = [0.5-10^δ,0.85,0.0,0.0]
	tspan = (0.0,tmax)
	g = 1.0
	sol1 = solve(ODEProblem(ballcirclehard!,u1,tspan,[g]), callback=ContinuousCallback(collision,bounce!);)
	sol2 = solve(ODEProblem(ballcirclehard!,u2,tspan,[g]), callback=ContinuousCallback(collision,bounce!);)
end;

# ╔═╡ 3a3ad2f0-bfa8-47cb-9d07-1d4215dcad18
begin
	plot(sol1,vars=(1,2),plotdensity=5000,xlims=(-1,1),ylims=(-1,1),label="u1")
	plot!(sol2,vars=(1,2),plotdensity=5000,xlims=(-1,1),ylims=(-1,1),label="u2")
	plot!(cos.(2*pi*(0:0.001:1)),sin.(2*pi*(0:0.001:1)),label="",size=(600,600))
	scatter!(sol1.u[end][1:1],sol1.u[end][2:2], 
	color=:blue,markersize=10,alpha=0.5,label="1")
	scatter!(sol2.u[end][1:1],sol2.u[end][2:2], color=:red,markersize=10,alpha=0.5,label="2")
end	

# ╔═╡ a5fcf13e-7ff9-4cd4-920d-e892e5f48a0c
begin
	ts = range(0, stop=tmax, length=5000)
	dd = (sol1(ts,idxs=2)-sol2(ts,idxs=2)).^2+(sol1(ts,idxs=1)-sol2(ts,idxs=1)).^2
	plot(ts,0.5*log10.(dd),size=(800,300),legend=false,xlabel="time",ylabel="log₁₀  distance")
end	

# ╔═╡ Cell order:
# ╠═1878f6d9-d9a7-4bf5-9525-2b021ebd1521
# ╠═794bf502-f9e9-11ec-03c5-ad174a6ccd59
# ╠═a830882a-8a52-4171-a0da-0449d4f5c901
# ╟─69c47e84-0bb9-425a-8cc8-3c50e2ffa2b3
# ╟─9969e088-322f-455b-8f3c-edb0ba8faed1
# ╠═6527dc0e-6c30-4d84-8f41-5ec3e42eaf86
# ╟─65bafd81-f8ff-49d4-a882-3007c8c8a71a
# ╟─507ddd8e-9609-40d0-bd4f-682f72475f93
# ╟─6e3ebd5c-b797-44f9-82ba-8ad9a18152eb
# ╟─8826d6ea-a2cb-430a-9e03-5fd459df24b6
# ╠═85011323-6e1b-4f28-8e94-9c9d2e3ae686
# ╠═295c9623-7657-4ec4-ab35-af8d068adacc
# ╠═55a5f7f0-e0fe-41dd-891a-9d733862dc9f
# ╟─d18add15-daa9-4a46-9a1a-9304e74459fe
# ╟─3a3ad2f0-bfa8-47cb-9d07-1d4215dcad18
# ╟─fa92921c-0d13-4e2f-ae97-37af89acb500
# ╟─a5fcf13e-7ff9-4cd4-920d-e892e5f48a0c
