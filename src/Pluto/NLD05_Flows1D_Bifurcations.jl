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

# ╔═╡ e886ebcd-7cbc-4fda-a208-13bebd50bed3
using Pkg;Pkg.add("Plots");Pkg.add("PlutoUI");Pkg.add("DifferentialEquations");Pkg.add("ForwardDiff"); Pkg.add("StaticArrays"); Pkg.add("IntervalRootFinding");

# ╔═╡ 1097d874-4468-4e76-bdbc-c893a5dbfdc0
using Plots, PlutoUI, DifferentialEquations, ForwardDiff, StaticArrays, IntervalRootFinding

# ╔═╡ d511c3c8-c596-48a5-8182-4dbcaa607eb6
include("../NLD_utils.jl")

# ╔═╡ 18a83e00-0dd4-4fe7-a5be-644361f875d3
TableOfContents()

# ╔═╡ 4ad5df20-f85d-11ec-3803-a7ed7ed8f2f3
md"""
# Flows 1D
## Definition and numerical solution

A 1D Flow is defined asigning the instantaneous rate of change of the state, i. e. the time derivative of the variable $\dot{x}$ to a vector field $f(x)$ which is a function of the variable only (the of the position in state space):

$\dot{x} = f(x)$

which has a particular solution $x(t)$ (or trajectory on state space) for a given initial condition $x(0)$. 

We can compute numerically this solution as in the cell below for the flow:

$\dot{x} = -2 x$
"""

# ╔═╡ a8b84ccb-977e-4cf7-b384-f2f68c48122a
function linear!(du,u,p,t)
	du[1]=-2*u[1]
end

# ╔═╡ d9bb5aa8-2708-40c3-8d42-3b21e6d634c0
md"""
x(0) $(@bind x0 Slider(-1.0:0.02:1.0,default=0.7;show_value=true)) 
"""

# ╔═╡ 1f8ae008-a195-42d8-b29e-66c78c81c8d6
begin
	sol = solve(ODEProblem(linear!,[x0],(0,3.0),[]))
	p1 = plot([-1,1],[0.0,0.0],c=:black,xlims=(-1,1),ylims=(-0.1,0.1),size=(800,100),yaxis=false,yticks=false,legend=false,xlabel="x")
	scatter!(p1,[x0],[0],c=:red)
	plot!(p1,[sol.u[1][1],sol.u[end][1]],[0.0,0.0],arrow=true,c=:red)
	p2 = plot(sol,ylabel="x",size=(800,300),ylims=(-1,1),legend=false)
	plot(p1,p2,layout=grid(2,1,heights=[0.2 ,0.8]))
end	

# ╔═╡ f317741a-f3ae-4450-b57d-8403cab18335
md"""
The solution (an exponential decay) $x(t)$ starting from the initial condition $x(0)=x_0$ is shown in the lower graph as a function of time, and as a trajectory (in red) in the state space (black line) above.

However we will be more interested in the topological structure of the flow that can be derived from the analysis of the vector field $f(x)$.

# Topological Structure of the Flow

The vector field $f(x)$ is defined assigning a vector (arrow) equal to the value of the function $f(x)$ to each of the points of the state space. For the case of the 1D Flow all regions with the vectors pointing in a given direction correspond to the same trajectory. 

The **fixed points** $x_*$ correspond to those states that evolve into themselves, i. e. where the vector field (and the time derivative) are null (zero).

**All trajectories in a 1D Flow start or end in a fixed point or at $\pm \infty$**

The topological structure of the flow is given by all the trajectories and fixed points of the system. 

In order to calculate this, there are two methods, one more graphical, the other more numerical. 

"""

# ╔═╡ d3d28638-8c26-461c-b804-d76a923515ff
md"""
### Graphical method

We plot $f(x)$ in the vertical axis superimposed with the state space in the $x$.

Then we determine the intervals of **positivity** ($f(x)>0$ green in the graph below) and **negativity** ($f(x)<0$ red in the graph below) of the function 

We draw trajectories going to the right in the intervals of positivity ($\dot{x}>0$) and to the left in the intervals of negativity ($\dot{x}<0$).

The **fixed points** are those where the flow is neither to the left nor to the right $f(x)=0$. 

A fixed point can be **stable** (or **atractor**) if the trajectories converge to it from the left and the right, or **unstable** (i. e. **repulsor**) if the trajectories diverge to both sides. all aother cases are called **neutral**.

"""

# ╔═╡ 7af7eecd-eb9e-4ab2-b2f8-890ad5969237
html"""
<div>
<img src="https://i.imgur.com/gIgFZOL.png" width="700px">
</div>
"""

# ╔═╡ c0aa4384-0807-42f6-86fd-e6182121e465
md"""
We also note that whenever the function crosses the horizontal axis (fixed point) with a negative slope the fixed point is an attractor (since a perturbation to the right goes back to the left and vice versa), while if the slope is positive the fixed point is a repulsor (a perturbation to the right keeps moving to the right and the same to the other side, the perturbation is amplified). This is the basis of the other method.
"""

# ╔═╡ 01b880a3-2e24-4177-b9b0-2931f028e6e9
md"""
### Numerical Method

This is not a numerical method to solve the differential equation but rather a method to analyze the vector field and determine the topological structure without the need to graph the entire function.

The first step is to determine all the fixed points $x_*$ that satisfies the equation $f(x_*)=0$. These are also called the **zeroes** of the function $f(x)$, the fixed points are the zeroes of the vector field $f(x)$. This can be done by hand if the function is simple. For example for the previous equation with $f(x)=-2x$ it is clear that the only $x_*$ that satifies $f(x_*)=0$ is $x_*=0$. If the vector field is not so simple there are many numerical methods to find the zeroes of a function. 

The second step is to determine the stability of the fixed points. This can be done by finding the slope of the function at the fixed point (or equivalently the **derivative of f with respect to x** at that point which we will denote using a tilde instead of a point $f'(x_*)$):

- If the slope is **positive** at the fixed point (or the zeros) the function (from left to right) is going from the negatives (the flow at the left is going to the left) to the positives (the flow at the right is going to the right), therefore the fixed point is **unstable**

- If the slope is **negative** at the fixed point (or the zeros) the function is going from the positives (the flow at the left is going to the right) to the negatives (the flow at the right is going to the left), therefore the fixed point is **stable**

- If the slope is **zero** we have a neutral fixed point. In this case it may happen, but not always, that the function tangentially touches the horizontal axis, in that case any slight perturbation of the shape of the function (later we will define what that means in a more formal way) either generates two fixed points (one stable and one unstable), or causes there to be no fixed point. Later we will see that this is linked to a type of change in the topology of the flow called repulsor-attractor bifurcation. 

Deriving the function and evaluating it at the fixed point is known as **linearization**. Why? Because for any function $f(x)$ the best linear approximation in a point environment is equal to the function evaluated at that point $f(x_*)$ (which is zero because we are at a fixed point) plus the derivative (slope of the tangent line to the function at that point that is the geometric interpretation of the derivative) multiplied by the deviation from the point $f'(x_*)(x-x_*)$, where we use the prime to denote the derivative with respect to the variable. Note that $f'(x_*)$ it is a constant because it is the derivative evaluated at the fixed point. 
In other words, we can approximate our nonlinear dynamical system as a linear system very close to the fixed point:

$\dot{x}=f'(x_*)(x-x_*)$

This means that in an environment of an attractor fixed point (if $f'(x_*)<0$ )  trajectories converge as a decreasing exponential to the attractor, while in a repulsor fixed point ($f'(x_*)<0$) the trajectories diverge exponentially, just as in a linear system.

If the derivative at the fixed point is zero the linealization is no longer valid.

We can now study the Logistic Equation using the graphical method
"""

# ╔═╡ 112c2bc7-352d-4025-8455-003d2543cd5c
md"""
# Logistic Equation

The logistic equation is defined as:

$\dot{x}=r_0x(1-x)$

Since we are not interested now in the parameters we fix $r_0=2$.

In the graph below we show in the left panel the vector field $f(x)$ in blue, superimposed witch the state space (in black) and a particular trajectory in red, and in the left panel the solution correspondig to this trajectory as a function of time.

The vector field is then $f(x) = 2 x (1-x)$. We must find the zeroes of these function that is already factorized. There are two posible solutions:

$x_* = 0$for this fixed point the slope is positive as can be seen in the graph (this can aso be computed analitycally from the derivatve function). Therefore this point is a repulsor.

$x_* = 1$ for this fixed point the slope is negative as can be seen in the graph (this can aso be computed analitycally from the derivatve function). This fixed point corresponds to the ideal population or capacity and is an atractor. All positive initial conditions tends towards it.


"""

# ╔═╡ d9808200-5320-49c1-9ebf-5c8dc128c6a3
function logistic(u,p,t)
	p[1]*u[1]*(1.0-u[1])
end;	

# ╔═╡ bf009c87-d4cd-4001-8bd2-c95e3581fc8d
md"""
x(0) $(@bind x02 Slider(0.0:0.02:2.0,default=0.3;show_value=true)) 
"""

# ╔═╡ 0e3c83d9-5268-44e4-8cee-dd2421943748
flow1D(logistic,x02,4.0,[2.0];xlims=[-0.1,2.0],ylims=[0,2.0])

# ╔═╡ e581cc6f-1a68-4b1b-a0e2-e0f3327334e6
md"""
# Parameter variation and structural stability

Now we can reintroduce the parameter concept. It is a value that remains constant during the time evolution, i.e. it does not vary, but it characterizes the flow and defines a family of "neighboring" flows in the space of all possible flows. It is very important to understand that although they cannot vary during the time evolution, we can be interested in **how the dynamical system varies** when we modify this parameter. In particular what will interest us most is if there are changes in the topological structure of the flow. Sometimes the parameters are called control parameters because they can be modified externally to obtain different flows.

Let's go to the example of the logistic equation and now let's incorporate two parameters: the growth rate $r$ and the population capacity $K$ (we will understand this name once we have analyzed the solutions and the flow structure).

This more "real" logistic equation can be written as

$\dot{x} = r_0 \left(1- \frac{x}{K}\right)x$

Let's see wha happens with the fixed points and their stability, which determine the topological structure of the flow, as we change the parameters.

"""

# ╔═╡ c8478515-c670-4d8f-bd22-336ebbe8b053
function logistic2(u,p,t)
	p[1]*u[1]*(1.0-u[1]/p[2])
end;	

# ╔═╡ ad65b408-4705-4bcc-b9bb-1009177117bf
md"""
x(0) $(@bind x03 Slider(0.0:0.02:2.0,default=0.3;show_value=true)) \
r0 $(@bind r0 Slider(0.02:0.02:2.0,default=0.3;show_value=true))
K $(@bind K Slider(0.02:0.02:2.0,default=1.0;show_value=true))
"""

# ╔═╡ 4b3851da-0bbf-46c8-99ec-b4fc7115982c
flow1D(logistic2,x03,20.0,[r0,K];xlims=[-0.1,2.0],ylims=[0,2.0])

# ╔═╡ 96df68ec-3d01-4eab-b270-b59d79e02f3a
md"""
Note that we vary the parameters in a realistic range, both the maximum growth rate $r_0$ and the population capacity $K$ are positive.

As can be seen, varying $r_0$ varies the rate at which we converge to the single attractor of the system, and varying $K$ varies the position of the attractor, which is precisely at $x_*=K$.

But what is more important, is that in all this range of variation **the topological structure of the flow does not change**. We always have a repeller on the left and an attractor on the right and that gives me only one possible structure, regardless of whether the position of the attractor changes.

This feature that the topological characteristics of the flows (the parametric family of flows) do not vary as the parameters vary is known as **structural stability**. In this case the logistic equation is structurally stable over the entire parameter domain, but we can also have structural stability in a range. 

And the most interesting thing happens when there is a loss of structural stability, because it implies that there will be a change in the topological structure of the flow, and this phenomenon is known as **bifurcation**. 
"""

# ╔═╡ ae8f5c80-2d38-4b6d-93ba-1ae45ddcfc86
md"""
# Qualitative Changes in the TS of the Flow																														
The 1D flow is characterized qualitatively by its fixed points and stability. Can these change? What is the most general form? 

Earlier we talked about structural stability and referred to the invariance of the topological structure to perturbations in the parameters. We had also wondered what happened when it tangentially grazed the horizontal axis at a fixed point and the linearization ceased to be valid. It is not difficult to see that when the latter happens we are faced with the imminence of a change in the topological structure of the flow. 

We can think of the function f(x) as a snaking rope that varies with the parameters. The topological structure of the flow will not change as long as we have the same number of crossings (number of fixed points) and the sequence of crossings up, down, up, etc. does not change. What is the generic way it can change? When a "loop" that is on the positive or negative side crosses the axis as seen in the sequence of drawings below in which we have three different system parameter values:																									
"""

# ╔═╡ 09c6c2ef-0b7f-460a-904c-f71ecf0f8fd2
html"""
<div>
<img src="https://i.imgur.com/dmDiLts.png" width="700px">
</div>
"""

# ╔═╡ e0bfc074-d0b0-42f3-8cf6-2ae5871f38c0
md"""
In this case at the exact point where the change occurs (center diagram) two things happen simultaneously:

- Two fixed points are merged into one. Any perturbation, however small it may be, on one side makes them disappear and on the other side generates two fixed points. Moreover, these two fixed points are always an attractor and a repeller.
- The derivative of the vector field is zero at that point so the linearization is no longer valid.

The conclusion is that there is no generic way to change the stability of the fixed points without creating or destroying a pair of fixed points: an attractor and a repeller. This is the attractor-repeller bifurcation (which is known as saddle node SN in systems of more than one dimension). In other words, a fixed point cannot go from being stable to unstable without passing through zero and this cannot be done in a generic way without collapsing with another fixed point of complementary stability (let's think that we try to invert the slope of one of the stretches of the black string above where it cuts the horizontal axis. Here "generic form" means that no additional conditions are imposed. If we force the origin to be always a fixed point or there is some symmetry condition there are other possible cases as we will see below.

Let us define a bifurcation more formally: **Any change in the number of fixed points or their stability (i.e. the topological structure of the flow) by varying the parameters of the system is known as a local bifurcation**. Local refers to occurring at a given point in phase space.
"""

# ╔═╡ a25549e0-8875-4b5c-aa1d-7e47e8b207e0
md"""

# Saddle Node Bifurcation

Now we will present the most general or "simple" model that presents a certain bifurcation that occurs at the origin and for a parameter value equal to zero. Particular systems that have a local bifurcation at some other point for some other parameter value can be brought to the normal form by a change of coordinates.

In the case of the repulsor attractor bifurcation (or SN in one dimension) it is fairly intuitive, both from the above drawing and from what we discussed when we wondered what happens to the linearization if the derivative of $f(x)$ becomes zero at the fixed point, that the simplest and most generic way to have a bifurcation of SN is by a (quadratic) parabola cutting the horizontal axis at the origin for the zero parameter value. That is:

$\dot{x} = a - x^2$

this is the version with the inverted parabola but we could equivalently have another equation with the positive quadratic term.

For this system we have two fixed points (the repeller on the left and the attractor on the right) when $a>0$ and no fixed point for $a<0$. In $a=0$ the parabola "kisses" the horizontal axis and the two fixed points collapse into one (which is attractor on the right and repulsor on the left).

The normal form gives us another important information that has to do with scaling. In a bifurcation environment, when the two fixed points with complementary stability are created they are located at $x_*=\pm\sqrt{a}$ . That is, they move away from the coalescence (or emergence) point with a scaling that goes as the square root of the parameter. When we find attractor repulsor (or saddle-node in higher dimensions) bifurcations there will always be an environment, for parameter values close to the value of the bifurcation, where this scaling is also observed, therefore we can anticipate the location of the fixed points without the need to find the zeros of the vector field.
"""

# ╔═╡ 6cd940fd-94ab-49c4-a29a-accdb2b81ec9
md"""
# Bifurcation Diagram

A representation that we are going to find very frequently and that it is extremely useful to characterize a given dynamic system is the bifurcation diagram in the product space of variables and parameters. 

In the case previos dynamical system we have a variable ($x$) and a parameter ($a$) and it is usual practice to represent the parameter on the horizontal axis (because it is in fact the independent control variable of our problem) and the variable on the vertical axis. That is to say that now the phase space will appear vertically. It is also usual not to represent the whole flow but only the fixed points. The fixed points are branches in the bifurcation diagram that give for each parameter value the location of the parameter (we could even write them as a function of the parameter, for example). Finally we can use blue color for the branches corresponding to the stable fixed points (or attractors) and red for the unstable ones (or repulsors).

Taking all this into account the bifurcation diagram for the normal form of the SN bifurcation is as follows:
"""

# ╔═╡ 3d502689-789e-408b-a2b2-bfeee64e5e4e
html"""
<div>
<img src="https://i.imgur.com/lmMKrvz.png" width="400px">
</div>
"""

# ╔═╡ 6c14fbbc-1ec3-4a80-884a-3e6121fa57da
md"""
The saddle-node bifurcation occurs at the point $a=0$, $x=0$ of the diagram (green point) from which arise to the right ($a>0$) a stable branch for the attractor (blue) and an unstable one for the repeller (red) that move away from the bifurcation point as \pm\sqrt{a} (the coordinates of the fixed points).
"""

# ╔═╡ 7557f8d7-2654-41c2-a443-9d5a76622501
md"""
x(0) $(@bind x04 Slider(-1.0:0.01:1.0,default=0.1;show_value=true)) 
a: $(@bind a Slider(-0.5:0.01:0.5,default=0.1;show_value=true)) 
"""

# ╔═╡ 3c64b9fd-c2b0-440f-84ac-b15c6574d2df
begin
	nodosilla(x,p,t)=p[1]-x*x
	flow1D(nodosilla,x04,2.0,[a],(u)->(u<-1.0);xlims=[-1.0,1.0],title="Saddle Node")
end	

# ╔═╡ f480e291-b67c-4a2f-9f44-914a340df81e
html"""
<style>
input[type*="range"] {
	width: 30%;
}
</style>
"""

# ╔═╡ Cell order:
# ╟─e886ebcd-7cbc-4fda-a208-13bebd50bed3
# ╠═1097d874-4468-4e76-bdbc-c893a5dbfdc0
# ╠═d511c3c8-c596-48a5-8182-4dbcaa607eb6
# ╟─18a83e00-0dd4-4fe7-a5be-644361f875d3
# ╟─4ad5df20-f85d-11ec-3803-a7ed7ed8f2f3
# ╠═a8b84ccb-977e-4cf7-b384-f2f68c48122a
# ╟─1f8ae008-a195-42d8-b29e-66c78c81c8d6
# ╟─d9bb5aa8-2708-40c3-8d42-3b21e6d634c0
# ╟─f317741a-f3ae-4450-b57d-8403cab18335
# ╟─d3d28638-8c26-461c-b804-d76a923515ff
# ╟─7af7eecd-eb9e-4ab2-b2f8-890ad5969237
# ╟─c0aa4384-0807-42f6-86fd-e6182121e465
# ╟─01b880a3-2e24-4177-b9b0-2931f028e6e9
# ╟─112c2bc7-352d-4025-8455-003d2543cd5c
# ╠═d9808200-5320-49c1-9ebf-5c8dc128c6a3
# ╠═0e3c83d9-5268-44e4-8cee-dd2421943748
# ╟─bf009c87-d4cd-4001-8bd2-c95e3581fc8d
# ╟─e581cc6f-1a68-4b1b-a0e2-e0f3327334e6
# ╠═c8478515-c670-4d8f-bd22-336ebbe8b053
# ╠═4b3851da-0bbf-46c8-99ec-b4fc7115982c
# ╟─ad65b408-4705-4bcc-b9bb-1009177117bf
# ╟─96df68ec-3d01-4eab-b270-b59d79e02f3a
# ╟─ae8f5c80-2d38-4b6d-93ba-1ae45ddcfc86
# ╟─09c6c2ef-0b7f-460a-904c-f71ecf0f8fd2
# ╟─e0bfc074-d0b0-42f3-8cf6-2ae5871f38c0
# ╟─a25549e0-8875-4b5c-aa1d-7e47e8b207e0
# ╟─6cd940fd-94ab-49c4-a29a-accdb2b81ec9
# ╟─3d502689-789e-408b-a2b2-bfeee64e5e4e
# ╟─6c14fbbc-1ec3-4a80-884a-3e6121fa57da
# ╟─3c64b9fd-c2b0-440f-84ac-b15c6574d2df
# ╟─7557f8d7-2654-41c2-a443-9d5a76622501
# ╟─f480e291-b67c-4a2f-9f44-914a340df81e
