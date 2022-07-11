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

# ╔═╡ c061212a-b7e5-420a-9fd2-350ff3ee807d
# ╠═╡ disabled = true
#=╠═╡
using Pkg;Pkg.add("Plots");Pkg.add("PlutoUI");Pkg.add("DifferentialEquations");Pkg.add("ForwardDiff"); Pkg.add("StaticArrays"); Pkg.add("IntervalRootFinding")
  ╠═╡ =#

# ╔═╡ 8601d8d7-d4df-473f-b65d-0f03aeb8f5f4
using PlutoUI, Plots, DifferentialEquations, ForwardDiff, IntervalRootFinding, StaticArrays

# ╔═╡ 09b22c30-bb22-4633-9939-2e97bb0beb5e
include("../NLD_utils.jl")

# ╔═╡ 1f49c325-fc99-4b15-81ee-dc1c5bbe6f08
gr();

# ╔═╡ 2153a4e0-5c90-487d-b35d-6b113dc95b78
md"""
# General form of the flow
A 2D flow in the plane of general form is given by the differential equations:

$\dot{x}=f(x,y)$

$\dot{y}=g(x,y)$

where $f$ and $g$ are two functions of two variables which can also be seen as a vector function or a vector field. The latter name comes from the fact that to each point in the phase space ($x,y$) (the field) corresponds a pair of values that form a tangent vector to the trajectory at that point, i.e. they point in the direction of the flow.

The name "flow" takes on a more literal sense in two dimensions, because the trajectories can be seen as the flow lines of a fluid moving in two dimensions and the vector field as the velocity field of that fluid. Of course this is a mental picture because real fluids obey their own equations while $f$ and $g$ are arbitrary functions. But continuing with the fluid analogy, the attracting fixed points could be interpreted as sinks of flow and the repellers as sources. And we will see that the attractor directions of the saddle points function as "watersheds" separating different "basins of attraction" leading to different boundary sets or different destinations.

## New features of 2D Flows

One-dimensional flows, although they can bifurcate, have a very limited repertoire of behaviors (basically diverging to infinity or converging to a fixed point).

Flows in 2D incorporate a new behavior and a new limit set (or equilibrium) in addition to the fixed points. The new behavior is the oscillation and the new equilibrium associated with the oscillation is called the **limit cycle**. This, as can be imagined, greatly expands the repertoire of possible behaviors and adds a new type of bifurcation that is not present in 1D flows (Hopf bifurcation).

Another new ingredient incorporated in 2D flows is that in addition to attracting and repelling fixed points, there are fixed points that are attractors in one direction and repellers in another, which is why they are called **saddle points** in allusion to the saddle shape that attracts in the front-back direction but repels to the sides. On the other hand, these directions in which they approach or move away from the saddle points will **function as organizers of the global flow**

The organization of the global flow is a final element that is added in two-dimensional flows, because topological changes can occur in this organization without local changes. All the bifurcations we have seen so far occur due to changes in the number and stability of fixed points and are known as **local bifurcations**. Bifurcations associated with changes in the topological organization of the flow that are not reducible to changes in fixed points are known as **global bifurcations**, and are much more difficult to study.

To illustrate this we are going to look at some examples of 2D dynamic systems without doing any analysis, just to see some ways of representing the flow in phase space.
"""

# ╔═╡ c9f6e916-0c30-4ae6-b56a-cafda01db350
md"""
# Duffing Oscillator

We are going to see a system with a slightly more varied repertoire of behaviors and that later (by adding a forced term) is going to be our guide to enter the chaotic systems so we are going to study it in detail.

This system is obtained by adding quadratic or cubic terms to the restoring force of an oscillator.

Starting from the linear equation of the Harmonic Oscilator

$\dot{x}=y$

$\dot{y}=-\gamma y - Kx$

we replace the linear force $-Kx$ by the nonlinear force $\beta x - x^3$.
In 1918 George Duffing published a paper systematically studying the influence of the different nonlinear terms in the oscillator, so this system is known as Duffing's oscillator.

A physical system (albeit somewhat artificial) that has this behavior can be assembled with a flexible metal tab placed between two magnets:

"""

# ╔═╡ a51cd197-9b85-4a17-91c2-44b99ec0d45b
html"""
<div>
<img src="https://i.imgur.com/7THzWtE.png" width="300px">
</div>
"""

# ╔═╡ dfd8837d-d437-4d2b-91ea-97c56c3ab150
md"""
If the power of the magnets does not exceed the rigidity of the tongue we have the situation of a single attractor in the middle (although it is not a free oscillation due to the presence of magnets). If the rigidity decreases (or the magnets get closer, or stronger magnets are placed) there are two situations of stable equilibrium (pointing to one or the other magnet) separated by an unstable point (which as we will see later is a saddle point).

The Duffing oscillator can be written by replacing $-Kx$ in the equation of the Harmonic Oscillator by the proposed force $\beta x-x^3$:

$\dot{x} = y$

$\dot{y} = -\gamma y + \beta x - x^3$

"""

# ╔═╡ f654e114-f3f7-4202-9451-eff4778d753d
md"""
In the case of the tab with the magnets $x$ is the horizontal position of the tip, $y$ the speed (hence the first differential equation), $\gamma$ is the dissipation of air and friction (as in the case of the spring is a force that opposes the speed) and in $\beta$ is summarized the relationship between the strength of the magnets and the stiffness of the tab. 

If $\beta$ is positive the magnets win (two attractors) and if $\beta$ is negative the stiffness wins (one attractor). Note that the cubic term is the one that always ends up winning, very far from equilibrium the force is always attractive, therefore the system will not explode.
"""

# ╔═╡ 34617d33-99c3-4b66-8440-1257ca2c0ce7
function duffing!(du,u,p,t)
    (γ,β)=p
    du[1] = u[2]
    du[2] = -γ*u[2]+u[1]*(β-u[1]*u[1])
    du
end;  

# ╔═╡ 7cfeb364-59cd-4c66-aa76-6a27ace2e68e
md"""
First we will explore the behavior of the Duffing oscillator for different initial conditions and parameter values and then we will move on to the analysis of the topological organization of the flow in the plane and its possible changes (bifurcations) with the parameters.
"""

# ╔═╡ 9b7ba5e5-0143-44c8-81a9-db10fa32d7e2
md"""
x(0) $(@bind x0 Slider(-1.0:0.02:1.0,default=0.7;show_value=true)) 
y(0) : $(@bind y0 Slider(-1.0:0.01:1.0,default=0.0;show_value=true)) \
γ : $(@bind γ Slider(0.0:0.01:1.0,default=0.7;show_value=true)) 
β : $(@bind β Slider(-1.0:0.01:1.0,default=0.7;show_value=true)) 
"""

# ╔═╡ eeb50040-b5fe-475f-954e-af70bfb48770
flow2d_vectorfield(duffing!,[x0,y0],50.0,[γ,β]; title="Duffing Oscillator. Vector Field")

# ╔═╡ 970330b1-e37e-4a8d-91f5-a91d365a82b0
md"""

# The plan of the qualitative analysis

The plan of the theory of dynamical systems for the qualitative analysis of 2D flows can be summarized in the following steps

1. determination of fixed points of the system
2. determine the stability class of the fixed points: (a) attractor node or focus, (b) repulsor node or focus, (c) saddle point.
3. calculate the stable varieties of the saddle points that organize the flow.
4. determine the existence of limit cycles
5. All the above conforms the topological organization of the flow in the plane (also called phase portrait), the most important part is to study how this topological organization can undergo changes as we vary one or more parameters.

In general there are analytical methods for points 1 and 2. In this case we will work graphically in most of the cases.
"""

# ╔═╡ fad1f89d-c3cd-4927-b3bb-c23814007665
md"""

# Fixed Points 

The fixed points will be those values ($x_*,y_*$) which make both equations equal to zero. The first equation is zero only for $y_*=0$. 

In order to make zero the second equation, knowing that $y$ must be zero, we need that:

$\beta x_* - x_*^3=0$

This is a cubic curve and in general we have one or three fixed points depending on the value of $\beta$
"""

# ╔═╡ 2def64ae-091e-4868-87ee-987cd6ab456f
html"""
<div>
<img src="https://i.imgur.com/4IMhywy.png" width="500px">
</div>
"""

# ╔═╡ 52c86644-3a24-4ea8-88ec-c259f6c1069a
md"""
It is therefore sufficient to solve

$\beta x_* = x_*^3$

which has a trivial solution $x^*_1=0$ 
 and then, only for the case of two symmetric solutions (which are obtained by dividing both members of the above equation by $x$ since it is different from 0) into
$x^*_{2,3}=\pm\sqrt{\beta}$
"""

# ╔═╡ 61908961-31f5-43e8-8330-bac3ce3fd274
md"""
# Nullclines and Flow Organization

The fixed points of a general 2D flow:

$\dot{x}=f(x,y)$

$\dot{y}=g(x,y)$

are going to be those points of the plane ($x,y$) that simultaneously cancel the functions $f$ and $g$. That is to say that in general we will have to find the pairs ($x_*,y_*$) that satisfy:

$f(x_*,y_*)=0$

$g(x_*,y_*)=0$

In the example we saw before we were able to calculate the fixed points.
But calculating the fixed points in more general cases will not be so simple. Let's introduce the concept of **nullclines** that will help us to determine the fixed points and give us some clues to the organization of the flow  and the type of stability of the fixed points.

The nullclines are the two curves determined in the plane by the conditions that one of the derivatives cancels separately.

That is to say:

- The first nullcline will be determined by all the points where $\dot{x}=0$ and at these points the vector field always points upwards or downwards (there is no horizontal flow just because $\dot{x}=0$), that is to say that **the trajectories always cut vertically to this nulcline**.
- The second nullcline will be determined by all the points where $\dot{y}=0$  and at these points the vector field always points to the right or left, i.e. **the trajectories always cut horizontally this second nulcline** (they have zero inclination in the vertical direction, hence the name).

The fixed points will be determined by the points where these two nullclines intersect, since at these points the two conditions $f(x_*,y_*)=0$ and $g(x_*,y_*)=0$ are satisfied simultaneously. In general, the nullclines will be curves in phase space. In the simplest case of a linear system they will be straight lines.

## Nullclines for the Duffing Oscillator 

For the case of the Duffing Oscillator, the first nullcline is the horizontal axis:

$y=0$

The condition for the second nullcline:

$-\gamma y + \beta x - x^3 = 0$

determines a cubic curve in the plane:

$y = \beta x/\gamma - x^3/\gamma$

In the graph below we plot the first nullcline in red and the second with a green line along with the vector fiels and one trajectory

"""

# ╔═╡ 1c1df8ca-382f-4887-80bf-0ae1dfc13e93
flow2d_nullclines(duffing!,[0.12;0.2],50.0,[0.15,0.5];vectorfield=true,regions=false,title="Duffing Oscillator. Vector Field and Nullclines")

# ╔═╡ d6dd2280-e6b5-41b4-a2d6-d08e71d833b2
md"""
The nullclines can give me additional information on how the flow is organized. Let's take the simplest case, the nullcline of $\dot{x}=0$
 which corresponds to the horizontal red line. That nullcline divides the phase space into two regions: above where all the flow goes to the right (because $\dot{x}>0$)  and below where all the flow goes to the left (because $\dot{x}<0$).
"""

# ╔═╡ af372f73-1ec6-4264-856b-444bf67bf76f
html"""
<div>
<img src="https://i.imgur.com/MIGLBCj.png" width="400px">
</div>
"""

# ╔═╡ eb478ae8-89f0-4eb2-ac7c-cd21db93e05e
md"""
In the case of the green curve in all the area below the inverted "N" curve the flow is pointing upwards $\dot{y}>0$ and at the other side of the curve the flow is pointing downwards $\dot{y}<0$.
"""

# ╔═╡ 3685dce3-5cba-4845-ac99-3c663711f12f
html"""
<div>
<img src="https://i.imgur.com/WecCj2e.png" width="400px">
</div>
"""

# ╔═╡ fe46ce5c-0add-4b9b-bb47-3bd64f962b27
md"""
If we combine both graphs we determine all the regions where the flow goes up-down-right-left and we have a very important clue about the global behavior of the flow.

To visualize this without having to calculate it by hand we will use the function flux2d_nullclines which in addition to plotting the nullclines curves $\dot{x}=0$
 in red and $\dot{y}=0$
 in green will differentiate the regions that divide these nullclines with transparent colors and will allow us to determine the direction of the flow in each of them.

We will use the following colors with transparency for the signs of each of the derivatives:

- green color for the flow to the right ($\dot{x}>0$)
- red color for the flow to the left ($\dot{x}<0$)
- yellow color for the upward flow ($\dot{y}>0$)
- blue color for downward flow ($\dot{y}<0$)

The combination of the colors will give four different types of regions:

- green + blue color = cyan where the flow goes to the right ($\dot{x}>0$) and downward ($\dot{y}<0$)
- color green + yellow = green where the flow goes to the right ($\dot{x}>0$) and up ($\dot{y}>0$)
- color red + yellow = orange where the flow goes to the left ($\dot{x}<0$) and up ($\dot{y}>0$)
- color red + blue = magenta where the flow goes to the left ($\dot{x}<0$) and down ($\dot{y}<0$)
"""

# ╔═╡ 8e14c209-4088-4570-b794-9df6c6e7e2a2
flow2d_nullclines(duffing!,[0.15,0.5];vectorfield=true,title="Duffing Oscillator. Nullclines")

# ╔═╡ 613fd8a4-f40b-4e9e-9e80-3c859f1fca79
md"""
# Stable and Unstable Manifolds of the Saddle Point

Another representation that is very useful for inferring the flow in the whole phase space is that of the stable and unstable manifolds of the saddle points. And they are defined by being the particular trajectories than converge to the saddle point going forward (stable) or backwards (unstable) in time and are invariant to the flow.

In particular the stable manifolds of the saddle points are relevant because they act as separatrices of the flow. To illustrate it better let's see it in the case of the flow above. We are only going to calculate the manifold of the saddle point for $\beta>0$.
"""

# ╔═╡ 30b42e19-b066-48aa-89d4-3ea79aabe330
md"""
γ : $(@bind γ2 Slider(0.0:0.01:1.0,default=0.7;show_value=true)) 
β : $(@bind β2 Slider(0.0:0.01:1.0,default=0.7;show_value=true)) 
"""

# ╔═╡ 8cd210e2-f014-4233-992f-d3565e9904a7
md"""
Note how the coiling and winding of the manifold is modified with the parameter $\gamma$. The two branches of the unstable manifolds spiral around each of the attractors and the branches of the stable manifold spiral outward around the three fixed points. 

For higher values of can be seen that the blue branches of the stable saddle manifold function as separatrices of the flow that will terminate at one or another fixed point. All the initial conditions that terminate in an attractor form the **basin of attraction** of the attractor. The stable manifolds are then forming the boundary between the two attractor basins.

"""

# ╔═╡ 7c99b3a9-5b92-4fa1-833b-1278c2e40cda
md"""
γ : $(@bind γ3 confirm(Slider(0.01:0.01:1.0,default=0.7;show_value=true))) \
β : $(@bind β3 confirm(Slider(0.01:0.01:1.0,default=0.7;show_value=true))) \
δ (resolution) $(@bind δ confirm(Slider(0.002:0.002:0.1,default=0.02;show_value=true))) \
"""

# ╔═╡ 96ff05e1-2d68-47fe-affc-9e17789f6305
attractor_basin(duffing!,[γ3,β3],[[sqrt(β3),0.0],[-sqrt(β3),0]],0.1;
    delta=δ,xlims=(-2.5,2.5),ylims=(-2.0,2.0),title="Duffing Oscillator. Attractor Basins")

# ╔═╡ 255dfea1-0594-43c1-9333-9056ddce3a06
function duffing_jac(u,p)
  J = Array{Float64, 2}(undef, 2, 2)
  J[1,1] = 0
  J[1,2] = 1.0
  J[2,1] = p[2]-2.0*u[1]*u[1]
  J[2,2] = -p[1]
  return J
end;

# ╔═╡ 7e365511-221e-4206-afd9-6dc99c302308
begin
	u0_array= [[-sqrt(β2);0],[0;0],[sqrt(β2);0]]	
	flow2d_manifolds(duffing!,duffing_jac,u0_array,[γ2,β2];
        tmax=100.0,delta=1e-5,xlims=[-2,2],ylims=[-2,2],title="Duffing Oscillator Stable and Unstable Manifolds")
end	

# ╔═╡ c7399360-fad1-466a-9324-f830b20b07a7
TableOfContents(title="📚 Table of Contents", indent=true, depth=4, aside=true)

# ╔═╡ 1f7f958e-b0c3-433f-bae9-3bf63da3de7a
html"""
<style>
input[type*="range"] {
	width: 30%;
}
</style>
"""

# ╔═╡ Cell order:
# ╟─c061212a-b7e5-420a-9fd2-350ff3ee807d
# ╠═8601d8d7-d4df-473f-b65d-0f03aeb8f5f4
# ╠═09b22c30-bb22-4633-9939-2e97bb0beb5e
# ╟─1f49c325-fc99-4b15-81ee-dc1c5bbe6f08
# ╟─2153a4e0-5c90-487d-b35d-6b113dc95b78
# ╟─c9f6e916-0c30-4ae6-b56a-cafda01db350
# ╟─a51cd197-9b85-4a17-91c2-44b99ec0d45b
# ╟─dfd8837d-d437-4d2b-91ea-97c56c3ab150
# ╟─f654e114-f3f7-4202-9451-eff4778d753d
# ╠═34617d33-99c3-4b66-8440-1257ca2c0ce7
# ╟─7cfeb364-59cd-4c66-aa76-6a27ace2e68e
# ╟─eeb50040-b5fe-475f-954e-af70bfb48770
# ╟─9b7ba5e5-0143-44c8-81a9-db10fa32d7e2
# ╟─970330b1-e37e-4a8d-91f5-a91d365a82b0
# ╟─fad1f89d-c3cd-4927-b3bb-c23814007665
# ╟─2def64ae-091e-4868-87ee-987cd6ab456f
# ╟─52c86644-3a24-4ea8-88ec-c259f6c1069a
# ╟─61908961-31f5-43e8-8330-bac3ce3fd274
# ╟─1c1df8ca-382f-4887-80bf-0ae1dfc13e93
# ╟─d6dd2280-e6b5-41b4-a2d6-d08e71d833b2
# ╟─af372f73-1ec6-4264-856b-444bf67bf76f
# ╟─eb478ae8-89f0-4eb2-ac7c-cd21db93e05e
# ╟─3685dce3-5cba-4845-ac99-3c663711f12f
# ╟─fe46ce5c-0add-4b9b-bb47-3bd64f962b27
# ╟─8e14c209-4088-4570-b794-9df6c6e7e2a2
# ╟─613fd8a4-f40b-4e9e-9e80-3c859f1fca79
# ╟─7e365511-221e-4206-afd9-6dc99c302308
# ╟─30b42e19-b066-48aa-89d4-3ea79aabe330
# ╟─8cd210e2-f014-4233-992f-d3565e9904a7
# ╟─96ff05e1-2d68-47fe-affc-9e17789f6305
# ╟─7c99b3a9-5b92-4fa1-833b-1278c2e40cda
# ╟─255dfea1-0594-43c1-9333-9056ddce3a06
# ╟─c7399360-fad1-466a-9324-f830b20b07a7
# ╟─1f7f958e-b0c3-433f-bae9-3bf63da3de7a
