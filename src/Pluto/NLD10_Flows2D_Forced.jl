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

# ╔═╡ fad1f89d-c3cd-4927-b3bb-c23814007665
md"""
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

# ╔═╡ 1c1df8ca-382f-4887-80bf-0ae1dfc13e93
flow2d_nullclines(duffing!,[0.12;0.2],50.0,[0.1,0.5];ylims=[-1.5,1.5],vectorfield=true,title="Duffing Oscillator")

# ╔═╡ 6b52504e-805e-417e-8fb2-6473b47b0ad0
md"""
# Forced Duffing Oscillator

Let us see now a periodic forced system where the response is much more irregular. We return to the Duffing oscillator (with linear friction). Recall that we arrived at this system by first writing the equation for the harmonic oscillator in its general form with a restoring force $K(x)$

$\dot{x}=y$ 

$\dot{y}=-\mu y + K(x)$

and choosing a restoring force with a linear and a cubic term:

$K(x) = \beta x - x^3$

As for large values of $x$ the cubic term will dominate, it is guaranteed that the system is globally attracting (if $x$ is positive $-K(x)$ is very negative and vice versa).

The Duffing oscillator is NOT a self oscillator because it has no negative friction (energy injection). In any case, what we are interested in studying here is the forced Duffing oscillator:

$\dot{x} = y$ 

$\dot{y} = -\mu y + \beta x -  x^3 + A cos(\phi)$ 

$\dot{\phi} = \omega$
"""

# ╔═╡ 6b14f6c6-536b-4de2-b669-e0f1f34dbbe2
function duffing_forced!(du,u,p,t)
    (γ,β,A,ω)=p
    du[1] = u[2]
    du[2] = -γ*u[2]+u[1]*(β-u[1]*u[1])+A*cos(u[3])
	du[3] = ω
    du
end;  

# ╔═╡ 4ee77c37-0347-45c0-b101-b96fcb7e004d
@bind pduff (
	PlutoUI.combine() do bind
		md"""
		γ: $(bind(Slider(0:0.01:1.0,default=0.8;show_value=true))) \
		β: $(bind(Slider(-2.0:0.02:2.0,default=1.0;show_value=true))) \
		A: $(bind(Slider(0.0:0.02:3.0,default=0.0;show_value=true))) \
		ω: $(bind(Slider(0:0.02:3.0,default=1.0;show_value=true))) \
		x0: $(bind(Slider(-2.0:0.02:2.0,default=1.0;show_value=true))) \
		ncycles: $(bind(Slider(1:10:200,default=1;show_value=true)))
		"""
	end
)	

# ╔═╡ 37ddd39f-4de1-4eb0-a230-46181c7e64db
flow2d_forced(duffing_forced!,[pduff[5],0.0,0.0],pduff,2*pi/pduff[4];tcycles=0,ncycles=pduff[6],xlims=(-2,2),ylims=(-1.5,1.5))

# ╔═╡ 8b29c214-4b49-4f44-9eb5-230899ec7321
md"""
# Attraction Basins

An important notion that will allow us to characterize the growing complexity of the behavior of this system is that of the basin of attraction.

The basin of attraction of a certain attractor (for example a limit cycle that in the Poincare section corresponds to a fixed point or a set of periodic points), is defined by all those initial conditions that converge to the attractor for long times.

To determine the basins it is necessary to evolve a grid of several thousand initial conditions over several cycles, therefore the graphs that follow can be very demanding. Run first with a low value of delta (the grid resolution) and make sure that Julia is taking in Threads.nthreads() the full number of processors in order to parallelize the problem.
"""

# ╔═╡ 0a270882-fde3-4dcc-8cef-291e86a4f8c2
md"""
Select one: $(@bind sel1 Select([([0.14,1.0,0.1,1.0],[[1.0,0.0],[-1.0,0.0],[-1,0.7],[0.3,0.3]])=>"γ=0.14,β=1,A=0.1,ω=1",([0.14,1.0,0.14,1.0],[[1.1,0.0],[-0.8,0.0]])=>"γ=0.14,β=1,A=0.14,ω=1",([0.14,1.0,0.2,1.0],[[1.2,0.0],[-0.9,0.0]])=>"γ=0.14,β=1,A=0.2,ω=1",([0.14,1.0,0.24,1.0],[[1.2,0.0],[-0.9,0.0],[0.1,1.1]])=>"γ=0.14,β=1,A=0.24,ω=1"]))
"""

# ╔═╡ dc6e3382-7239-48f1-a88d-bc12fda7cb73
attractor_basin(duffing_forced!,sel1[1],sel1[2],0.3;delta=0.01,tmax=30*pi,xlims=(-2.5,2.5),ylims=(-2.0,2.0))

# ╔═╡ 3078878b-a2f8-40f5-9398-401682ebb2f5
md"""
# Strange Attractor

For greater values of $A$, the trajectories no longer converge to limit cycles but istead they approach to a set with a fractal structure known as a **strange attractor**. On the Poincare section is a set of points that, unlike the one generated by a torus formed by quasiperiodic orbits, is not confined to a curve. Let's see the poincare section for a value $A=0.27$ that gives rise to a strange attractor:
"""

# ╔═╡ 59c6b45e-3959-45ff-be82-4c0409d17589
@bind par (
	PlutoUI.combine() do bind
		md"""
		μ: $(bind(Slider(0:0.01:0.3,default=0.14;show_value=true))) \
		β: $(bind(Slider(-2.0:0.02:2.0,default=1.0;show_value=true))) \
		A: $(bind(Slider(0.0:0.02:3.0,default=0.0;show_value=true))) \
		ω: $(bind(Slider(0:0.02:3.0,default=1.0;show_value=true))) \
		ncycles: $(bind(Slider(1000:1000:90000,default=10000;show_value=true)))
		"""
	end
)	

# ╔═╡ 64e04ff2-6662-4775-93cd-85a915487478
poincare_forced(duffing_forced!,[0.5,0.5,0],par,2*pi/par[4]; tcycles=30,ncycles=par[5],size=(900,600))

# ╔═╡ c7399360-fad1-466a-9324-f830b20b07a7
TableOfContents(title="📚 Table of Contents", indent=true, depth=4, aside=true)

# ╔═╡ 1f7f958e-b0c3-433f-bae9-3bf63da3de7a
html"""
<style>
input[type*="range"] {
	width: 50%;
}
</style>
"""

# ╔═╡ Cell order:
# ╟─c061212a-b7e5-420a-9fd2-350ff3ee807d
# ╠═8601d8d7-d4df-473f-b65d-0f03aeb8f5f4
# ╠═09b22c30-bb22-4633-9939-2e97bb0beb5e
# ╟─1f49c325-fc99-4b15-81ee-dc1c5bbe6f08
# ╟─c9f6e916-0c30-4ae6-b56a-cafda01db350
# ╟─a51cd197-9b85-4a17-91c2-44b99ec0d45b
# ╟─dfd8837d-d437-4d2b-91ea-97c56c3ab150
# ╟─f654e114-f3f7-4202-9451-eff4778d753d
# ╠═34617d33-99c3-4b66-8440-1257ca2c0ce7
# ╟─fad1f89d-c3cd-4927-b3bb-c23814007665
# ╟─2def64ae-091e-4868-87ee-987cd6ab456f
# ╟─52c86644-3a24-4ea8-88ec-c259f6c1069a
# ╠═1c1df8ca-382f-4887-80bf-0ae1dfc13e93
# ╟─6b52504e-805e-417e-8fb2-6473b47b0ad0
# ╠═6b14f6c6-536b-4de2-b669-e0f1f34dbbe2
# ╟─37ddd39f-4de1-4eb0-a230-46181c7e64db
# ╟─4ee77c37-0347-45c0-b101-b96fcb7e004d
# ╟─8b29c214-4b49-4f44-9eb5-230899ec7321
# ╟─0a270882-fde3-4dcc-8cef-291e86a4f8c2
# ╟─dc6e3382-7239-48f1-a88d-bc12fda7cb73
# ╟─3078878b-a2f8-40f5-9398-401682ebb2f5
# ╟─64e04ff2-6662-4775-93cd-85a915487478
# ╟─59c6b45e-3959-45ff-be82-4c0409d17589
# ╟─c7399360-fad1-466a-9324-f830b20b07a7
# ╟─1f7f958e-b0c3-433f-bae9-3bf63da3de7a
