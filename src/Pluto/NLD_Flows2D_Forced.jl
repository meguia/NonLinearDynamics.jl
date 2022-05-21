### A Pluto.jl notebook ###
# v0.19.5

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

# ╔═╡ 47e9a220-d944-11ec-1cb8-e7f659090d81
using Pkg; Pkg.activate("../..")

# ╔═╡ 8601d8d7-d4df-473f-b65d-0f03aeb8f5f4
using PlutoUI, Plots, DifferentialEquations, NonLinearDynamics

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

# ╔═╡ 4237b1cb-075a-49be-8f2c-a6eb9c5c9901
md"""
# Forced van der Pol

$\dot{x} = y$

$\dot{y} = \mu (1 -x^2)y - x + A cos(\phi)$

$\dot{\phi} = \omega$
"""

# ╔═╡ 906dd4e2-b975-4eb4-8469-a25ea89668dc
function fvdp!(du,u,p,t)
    du[1] = u[2]
    du[2] = p[1]*(1.0-u[1]*u[1])*u[2]-u[1]+p[2]*cos(u[3])
    du[3] = p[3]
    du
end    

# ╔═╡ 6282f7e6-4831-45ef-9b1d-3f077fdd74a4
@bind pvdp (
	PlutoUI.combine() do bind
		md"""
		μ: $(bind(Slider(0:0.02:4.0,default=0.1;show_value=true))) \
		A: $(bind(Slider(0.0:0.02:3.0,default=0.0;show_value=true))) \
		ω: $(bind(Slider(0:0.02:3.0,default=1.0;show_value=true))) \
		ncycles: $(bind(Slider(1:200,default=12;show_value=true)))
		"""
	end
)	

# ╔═╡ 2d29b217-fe25-4bb1-9282-51b4932356cf
flow2d_forced(fvdp!,[0.5,0.5,0],pvdp,2*pi/pvdp[3]; tcycles=10,ncycles=pvdp[4])

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

# ╔═╡ 4c3108b9-2328-42b7-9c69-bee0647e94b8
function duffing_forced!(du,u,p,t)
    (μ,β,A,ω)=p
    du[1] = u[2]
    du[2] = -μ*u[2]+u[1]*(β-u[1]*u[1])+A*cos(u[3])
    du[3] = ω
    du
end  

# ╔═╡ 4ee77c37-0347-45c0-b101-b96fcb7e004d
@bind pduff (
	PlutoUI.combine() do bind
		md"""
		μ: $(bind(Slider(0:0.02:4.0,default=0.8;show_value=true))) \
		β: $(bind(Slider(-2.0:0.02:2.0,default=1.0;show_value=true))) \
		A: $(bind(Slider(0.0:0.02:3.0,default=0.0;show_value=true))) \
		ω: $(bind(Slider(0:0.02:3.0,default=1.0;show_value=true))) \
		x0: $(bind(Slider(-2.0:0.02:2.0,default=1.0;show_value=true))) \
		ncycles: $(bind(Slider(1:200,default=1;show_value=true)))
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

# ╔═╡ 851d8073-e583-4c08-b53a-838dbe19eaf8
# Select One
#p = [0.14,1.0,0.1,1.0]; attractors=[[1.0,0.0],[-1.0,0.0],[-1,0.7],[0.3,0.3]];
p = [0.14,1.0,0.14,1.0]; attractors=[[1.1,0.0],[-0.8,0.0]];
# p = [0.14,1.0,0.2,1.0]; attractors=[[1.2,0.0],[-0.7,0.0]];
# p = [0.14,1.0,0.24,1.0]; attractors=[[1.2,0.0],[-0.7,0.0],[0.1,1.1]];

# ╔═╡ dc6e3382-7239-48f1-a88d-bc12fda7cb73
p1=attractor_basin(duffing_forced!,p,attractors,0.3;delta=0.01,tmax=30*pi,xlims=(-2.5,2.5),ylims=(-2.0,2.0))

# ╔═╡ 3078878b-a2f8-40f5-9398-401682ebb2f5
md"""
# Strange Attractor

For greater values of $A$, the trajectories no longer converge to limit cycles but istead they approach to a set with a fractal structure known as a **strange attractor**. On the Poincare section is a set of points that, unlike the one generated by a torus formed by quasiperiodic orbits, is not confined to a curve. Let's see the poincare section for a value $A=0.27$ that gives rise to a strange attractor:
"""

# ╔═╡ 59c6b45e-3959-45ff-be82-4c0409d17589
@bind par (
	PlutoUI.combine() do bind
		md"""
		μ: $(bind(Slider(0:0.02:0.5,default=0.14;show_value=true))) \
		β: $(bind(Slider(-2.0:0.02:2.0,default=1.0;show_value=true))) \
		A: $(bind(Slider(0.0:0.02:3.0,default=0.0;show_value=true))) \
		ω: $(bind(Slider(0:0.02:3.0,default=1.0;show_value=true))) \
		ncycles: $(bind(Slider(1000:1000:100000,default=10000;show_value=true)))
		"""
	end
)	

# ╔═╡ 64e04ff2-6662-4775-93cd-85a915487478
poincare_forced(duffing_forced!,[0.5,0.5,0],par,2*pi; tcycles=30,ncycles=par[5],size=(900,600))

# ╔═╡ Cell order:
# ╠═47e9a220-d944-11ec-1cb8-e7f659090d81
# ╠═8601d8d7-d4df-473f-b65d-0f03aeb8f5f4
# ╟─c7399360-fad1-466a-9324-f830b20b07a7
# ╟─1f7f958e-b0c3-433f-bae9-3bf63da3de7a
# ╟─4237b1cb-075a-49be-8f2c-a6eb9c5c9901
# ╠═906dd4e2-b975-4eb4-8469-a25ea89668dc
# ╠═2d29b217-fe25-4bb1-9282-51b4932356cf
# ╟─6282f7e6-4831-45ef-9b1d-3f077fdd74a4
# ╟─6b52504e-805e-417e-8fb2-6473b47b0ad0
# ╠═4c3108b9-2328-42b7-9c69-bee0647e94b8
# ╠═37ddd39f-4de1-4eb0-a230-46181c7e64db
# ╟─4ee77c37-0347-45c0-b101-b96fcb7e004d
# ╟─8b29c214-4b49-4f44-9eb5-230899ec7321
# ╠═851d8073-e583-4c08-b53a-838dbe19eaf8
# ╠═dc6e3382-7239-48f1-a88d-bc12fda7cb73
# ╟─3078878b-a2f8-40f5-9398-401682ebb2f5
# ╠═64e04ff2-6662-4775-93cd-85a915487478
# ╟─59c6b45e-3959-45ff-be82-4c0409d17589
