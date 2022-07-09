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

# ╔═╡ 977dfa56-983b-4b6d-b4dd-69786b2d0b1b
using Pkg;Pkg.add("Plots");Pkg.add("PlutoUI");Pkg.add("SIMD");

# ╔═╡ 1097d874-4468-4e76-bdbc-c893a5dbfdc0
using Plots, PlutoUI, SIMD

# ╔═╡ 18a83e00-0dd4-4fe7-a5be-644361f875d3
TableOfContents()

# ╔═╡ 4ad5df20-f85d-11ec-3803-a7ed7ed8f2f3
md"""
# Dynamical System

A dynamical system is a system evolves in time following a **fixed deterministic** rule.

At a given time $t$ the **state** of the dynamical system is completely characterized by the value of one or more **variables**.

We start studying systems with a single variable $x$. The **space of states** (all posible states) is the real number line. 
"""

# ╔═╡ 157c155f-241e-4d39-93f7-35997d83191d
md"""
## Discrete Time: $x_t \rightarrow x_{t+1}$

If the time evolution is discrete:

$t=0,1,2,3,4,5,6\dots$

we have a **map** and the dynamical system is defined by a rule that assigns a next state of the system at time $t+1$ ($x_{t+1}$) as a function of the present state at time $t$ ($x_t$) only, and this is denoted by:

$x_{t+1}=f(x_t)$

# Simplest example: Continuous Population with Discrete Time

Take for example the number of individuals of a population at each generation, and assume that the rate of growth $R$ is fixed. 
Also, assume that since the population is large enough it can be approximated by a continuous value  $x$ (can be also obtained taking some kind of average).
Then the rule of our dynamical system can be written as:

$x_{t+1}=R x_t$

Then, at each time step the population number is multiplied by $R$.

In the animated graphic below we can set the initial value of the population $x_0$ and the rate of growth $R$ with the sliders, and also start the iteration of the map controlling the number of seconds for each step (tick)

"""

# ╔═╡ 9490e7f8-668e-4acc-8efd-dcc4d355acb9
md"""
$(@bind n Clock(1))
$(@bind reset Button("Reset"))
"""

# ╔═╡ d9bb5aa8-2708-40c3-8d42-3b21e6d634c0
md"""
x₀ : $(@bind x0 Slider(0.0:0.02:1.0,default=0.1;show_value=true)) \
R : $(@bind R Slider(0.05:0.01:1.5,default=1.0;show_value=true)) \
"""

# ╔═╡ dd1af6a9-c048-4845-8858-8051d23f6bee
md"""
The upper graphics displays the states of the dynamical system given by the variable $x$ (average number of individuals in the population) for each generation (time step $t$) as red dots in the state space (black line). 

The sequence of red dots ($x_0,x_1,x_2,x_3,x_4,x_5,\dots$) is a **trajectory** of the dynamical system. 

The lower graph displays the value of the variable $x$ as a function of time $t$ as blue dots. This is the **solution** of the system.
"""

# ╔═╡ 226bcb3d-fdfa-484d-80da-2eb1e778017b
md"""
Each initial condition $x_0$ leads to an unique trajectory in the dynamical system.

On the other hand, $R$ is an **external parameter**  that determines the rule of the evolution and does not change over time.

Starting from $x_0$ for generation $t=0$, the population evolves in time either growing without bounds (for $R>1$) or decreasing toward zero (for $R<1$). A marginal case occurs when $R=1$ and the population remains fixed at the initial value. 

Therefore, the parameter $R$ determines the qualitative behavior of the system (population growth or decrease) for all initial conditions, while the initial condition correspond to a particular history of the system (trajectory).

The qualitative theory of dynamical systems is mainly concerned with the possible qualitative changes in the behavior of the systems as the parameters are varied.
"""

# ╔═╡ 6a97c988-3ec8-4271-b61c-47c40b1deee6
md"""
# Logistic Map

The growth of a population with a rate greater than one ($R>1$) is not sustainable over time.
Therefore, one of the first approximations is to consider a drop in the growth rate with the population, which takes its maximum value when the population is very small and becomes zero for a maximum population value. The logistic equation (the name is historical) proposes a growth rate: $R(x)=r(1-x)$ that becomes zero when $x=1$ (maximum population). If we replace in the previous map $R$ by $r(1-x)$ we obtain the logistic map.

$x_{t+1} = r(1-x_t)x_t$
"""

# ╔═╡ a31644b5-3aed-4ff8-9f46-0f8ff407f1b7
md"""
$(@bind n2 Clock(1))
$(@bind reset2 Button("Reset"))
"""

# ╔═╡ 1115e7e3-06f2-4f51-a71d-85ada76ce3a9
md"""
x₀ : $(@bind x02 Slider(0.0:0.02:1.0,default=0.1;show_value=true)) \
r : $(@bind R2 Slider(0.05:0.01:4.0,default=1.0;show_value=true))
"""

# ╔═╡ 25875787-048f-4baa-8ccc-6dececb097f9
md"""
For $r<1.0$ the logistic map displays a qualitative similar behavior to the simplest growth map: all initial conditions converge to $x_*=0$. This is a **fixed point** and an **atractor** of the map (we will use the $_*$ subindex to denote fixed points)

It is a **fixed point** because an an initial condition starting at $x=0$ remains there. The fixed point **maps into itself**. The fixed points are also **limit sets** in the sense that trajectories converge to them in the limit of $t \rightarrow \infty$

It is also an **attractor** because all neighboring initial conditions approach its value as time evolves (converge).

For $r$ between 1 and 3 all initial conditions (except $x_*=0$ which is still a fixed point but it is not an attractor anymore) converge to an intermediate value. This is also a fixed point and attractor. For example, for $r=2$ the fixed point is at $x_*=0.5$ (check it out!).

For values of $r$ slightly above 3 a new kind of behavior arises. The trajectory, in the long term, tend to **oscillate** between two (or four, eight, sixteen, ...) values. For example, for $r=3.2$ these two values are approximately $x_1\approx 0.51135$ and $x_2 \approx 0.80115$. This limit set is called a **periodic orbit** of the system and is also an **attractor** because all neighboring orbits converge to these points.

Finally for $r>3.56995...$ the logistic map displays a new kind of behavior that  neither converge to a fixed point nor approaches a periodic orbit. 
Except for some island of periodic behavior (check for example $r=3.83$), the regime in this region is **chaotic** and the limit set that all trajectories tend to is called a **chaotic attractor**.
"""

# ╔═╡ 1399a183-b4a0-46bf-934e-5e52a553cc67
md"""
## Bifurcation Diagram.

An extremelly useful representation of the different kind of behaviors of the logistic map as the parameter $r$ is varied is the **bifurcation diagram** of the system (the name will become clearer after we define a bifurcation).

In this representation, the **state space** is placed in the vertical axis and we perform multiple realizations of the logistic map for different $r$ values along the horizontal line. We then evolve an initial condition and wait until it gets close enough to the limit set (which can be a fixed point, a periodic orbit or a chaotic attractor) and plot the result.

In the graph below we display the bifurcation diagram of the logistic map for values of $r$ between $rmin$ and $rmax$ (note that the button "submit query" must be pushed in order to refresh the plot). All the limit sets for the differenta values of $r$ are displayed in blue, while the particular realization for the value of $r$ chosen by the sliders above are displayed in red.
"""

# ╔═╡ 9879a686-0e3f-4fe4-8f74-9a7613e89ccf
@bind output confirm(
	PlutoUI.combine() do bind
		md"""
		rmin $(bind(Slider(0.0:0.001:4.0,default=0.0;show_value=true))) \
		rmax $(bind(Slider(0.0:0.001:4.0,default=4.0;show_value=true)))\
		#points $(bind(Slider(100:100:4000,default=2000;show_value=true)))\
		#iterations $(bind(Slider(500:500:4000,default=1000;show_value=true)))\
		"""
	end
)

# ╔═╡ b2cb61e0-21c0-404b-8418-59c995d8d852
function bifur!(rs, x, xs, maxiter)
    x[:,1] .= xs
    for j in 1:maxiter-1
        @inbounds @simd for k in 1:length(rs) 
            x[k, j+1] =  rs[k] * x[k, j] * (1 - x[k,j])
        end
    end
    return x
end;

# ╔═╡ 15d2ab3e-2818-4823-8544-9adc4410cc1a
begin
	xs = 0.5
	rini,rend,Nrs,maxiter = output
	trans = Int(maxiter/2)
	rs = LinRange(rini, rend, Nrs+1) 
	x = zeros(length(rs), maxiter);
	bifur!(rs, x, xs, maxiter);
	p3 = plot(rs[:], x[:, trans:maxiter],seriestype = :scatter, markercolor=:blue, markersize = 1, markerstrokewidth = 0, legend = false, markeralpha = 0.1,size=(900,700))
end;	

# ╔═╡ bcd5ac86-828f-4c65-9885-e5ff5b6790f0
begin
	nr = findfirst(rs .> R2)
	if !isnothing(nr) && nr>1
		plot!(p3,rs[nr:nr],x[nr:nr,1:maxiter],seriestype = :scatter, markercolor=:red, markersize = 3, markerstrokewidth = 0, markeralpha = 0.3)
	else
		p3
	end
end	

# ╔═╡ 3a6821d7-392a-4dd4-b5c8-c247edf9132f
mutable struct Map
	x 
	t
	xv 
end	

# ╔═╡ ac175af3-0d5c-43ab-a41f-a5e6dfc1df72
begin
	reset
	tmax=1000
	m = Map(x0/R,-1,Array{Float64,1}(undef,tmax));
end;

# ╔═╡ 7c3fb7dc-47f6-4a6c-b42a-76f74a5bc575
begin
	reset2
	tmax2=2000
	m2 = Map(x02,0,Array{Float64,1}(undef,tmax2));
end;

# ╔═╡ 54918524-a078-497d-85ce-2bff7908841f
function plot_state(x,xrange,t)
	title = "t=$t   xₜ = $(round(x,digits=8))"
	p1 = plot(xlims=xrange,ylims=(-0.1,0.1),size=(800,100),yaxis=false,yticks=false,title=title,legend=false)
	plot!(p1,[xrange[1],xrange[2]],[0,0],c=:black)
	#scatter!(p1,[x],[0],c=:red)
end;

# ╔═╡ d41c5669-adab-4954-9563-b0312dc5077d
p1 = plot_state(m.x,(0,2.0),m.t);

# ╔═╡ 486aa75c-816a-47f6-8e68-aa260decc5d7
p2 = plot_state(m2.x,(0,1.0),m2.t);

# ╔═╡ 2f3f5e94-434e-46bc-94f2-0c4f8eda284d
function plot_state!(p1,x,t;alpha=1.0)
	title = "t=$t   xₜ = $(round(x,digits=8))"
	scatter!(p1,[x],[0],c=:red,alpha=alpha)
	title!(p1,title)
end;

# ╔═╡ 936816f9-5a79-474c-8847-96a1f0cb9a6a
begin
	n
	m.x = R * m.x
	m.t += 1
	if m.t>tmax
		m.t = 0
		reset
	end	
	m.xv[m.t+1]=m.x
	p1b = plot(m.xv[1:m.t+1],marker=4,size=(800,300),xlabel="t",ylabel="xₜ",legend=false)
	plot_state!(p1,m.x,m.t)
end

# ╔═╡ 4c273264-9c50-4f57-8136-b0756e116de0
p1b

# ╔═╡ cc09e3b4-8eca-40e0-a2ec-f8d58fe8cc68
begin
	n2
	m2.x = R2 * m2.x * (1.0-m2.x)
	m2.t += 1
	if m2.t>tmax2
		m2.t = 0
		reset2
	end	
	m2.xv[m2.t]=m2.x
	p2b = plot(m2.xv[1:m2.t],marker=4,size=(800,300),xlabel="t",ylabel="xₜ",legend=false)
	plot_state!(p2,m2.x,m2.t;alpha=0.2)
	plot(p2,p2b,layout=grid(2,1,heights=[0.2 ,0.8]))
end

# ╔═╡ f480e291-b67c-4a2f-9f44-914a340df81e
html"""
<style>
input[type*="range"] {
	width: 60%;
}
</style>
"""

# ╔═╡ Cell order:
# ╠═977dfa56-983b-4b6d-b4dd-69786b2d0b1b
# ╠═1097d874-4468-4e76-bdbc-c893a5dbfdc0
# ╟─18a83e00-0dd4-4fe7-a5be-644361f875d3
# ╟─4ad5df20-f85d-11ec-3803-a7ed7ed8f2f3
# ╟─157c155f-241e-4d39-93f7-35997d83191d
# ╟─9490e7f8-668e-4acc-8efd-dcc4d355acb9
# ╟─d9bb5aa8-2708-40c3-8d42-3b21e6d634c0
# ╟─936816f9-5a79-474c-8847-96a1f0cb9a6a
# ╟─dd1af6a9-c048-4845-8858-8051d23f6bee
# ╟─4c273264-9c50-4f57-8136-b0756e116de0
# ╟─d41c5669-adab-4954-9563-b0312dc5077d
# ╟─ac175af3-0d5c-43ab-a41f-a5e6dfc1df72
# ╟─226bcb3d-fdfa-484d-80da-2eb1e778017b
# ╟─6a97c988-3ec8-4271-b61c-47c40b1deee6
# ╟─a31644b5-3aed-4ff8-9f46-0f8ff407f1b7
# ╟─1115e7e3-06f2-4f51-a71d-85ada76ce3a9
# ╟─cc09e3b4-8eca-40e0-a2ec-f8d58fe8cc68
# ╟─25875787-048f-4baa-8ccc-6dececb097f9
# ╟─1399a183-b4a0-46bf-934e-5e52a553cc67
# ╟─bcd5ac86-828f-4c65-9885-e5ff5b6790f0
# ╟─9879a686-0e3f-4fe4-8f74-9a7613e89ccf
# ╟─486aa75c-816a-47f6-8e68-aa260decc5d7
# ╟─15d2ab3e-2818-4823-8544-9adc4410cc1a
# ╟─b2cb61e0-21c0-404b-8418-59c995d8d852
# ╟─7c3fb7dc-47f6-4a6c-b42a-76f74a5bc575
# ╟─3a6821d7-392a-4dd4-b5c8-c247edf9132f
# ╟─54918524-a078-497d-85ce-2bff7908841f
# ╟─2f3f5e94-434e-46bc-94f2-0c4f8eda284d
# ╟─f480e291-b67c-4a2f-9f44-914a340df81e
