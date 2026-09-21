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

# ╔═╡ 8601d8d7-d4df-473f-b65d-0f03aeb8f5f4
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using NonLinearDynamics: attractor_basin, flow2d_forced, poincare_forced, poincare_forced!, poincare_forced_zoom, saddle_orbit2D, saddle_manifolds_forced
    using PlutoUI, Plots, DifferentialEquations, ForwardDiff, IntervalRootFinding, StaticArrays
end

# ╔═╡ 1f49c325-fc99-4b15-81ee-dc1c5bbe6f08
gr();

# ╔═╡ 6b52504e-805e-417e-8fb2-6473b47b0ad0
md"""
# Forced Duffing Oscillator

Let us see now a periodic forced system where the response is much more irregular. We return to the Duffing oscillator (with linear friction). Recall that we arrived at this system by first writing the equation for the harmonic oscillator in its general form with a restoring force $K(x)$

$\dot{x}=y$ 

$\dot{y}=-\mu y + K(x)$

and choosing a restoring force with a linear and a cubic term:

$K(x) = \beta x - x^3$

For large $|x|$ the cubic term dominates and the restoring force points toward the centre. With positive damping the unforced system dissipates energy; at zero damping this does not imply attraction.

The Duffing oscillator is NOT a self oscillator because it has no negative friction (energy injection). 

In the case of the Duffing oscillator we introduce the forcing as a periodic force with two control parameters: the angular frequency $\omega$ and the amplitude of the forcing $A$. Since it is a force, it is added as a term in the second equation, together with the restoring force and the nonlinear dissipation . 

However, before writing the equations it is worth remembering that when we defined a dynamical system we said that the evolution rules, given by the differential equations, were fixed and did not change in time. Yet, we can always write a time-dependent system by 'inventing' time as a new variable, or in the case of a periodic function, the phase of the forcing $\phi$. Thus we can write the forced Duffing oscillator as follows:

$\dot{x} = y$ 

$\dot{y} = -\mu y + \beta x -  x^3 + A cos(\phi)$ 

$\dot{\phi} = \omega$

at the cost of adding one more dimension, the flow is now 3D. The evolution of the third variable is trivial since as it corresponds to the phase it advances linearly in time as $\omega t$. 

Since this variable enters in the first two equations only through a sine function we can consider the third dimension to be periodic (formally we could consider the space formed by the product of a plane and a circle in the third variable). This sound strange but it can be represented as an infinite horizontal space between a ground and a ceil, and making the identity between the ground and the ceil, so every trayectory that goed through the ceil appears immediately at the ground in  the same point y viceversa. 

If we project in the plane we must take into account that now on the phase space (or rather on its projection) we are going to see trajectories that cross but that correspond to two different coordinates.
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
		γ: $(bind(Slider(0:0.01:1.0,default=0.8;show_value=true))) 
		β: $(bind(Slider(-2.0:0.02:2.0,default=1.0;show_value=true))) \
		A: $(bind(Slider(0.0:0.02:3.0,default=0.0;show_value=true))) 
		ω: $(bind(Slider(0.02:0.02:3.0,default=1.0;show_value=true))) \
		x0: $(bind(Slider(-2.0:0.02:2.0,default=1.0;show_value=true))) 
		ncycles: $(bind(Slider(1:10:200,default=1;show_value=true)))
		"""
	end
)	

# ╔═╡ 37ddd39f-4de1-4eb0-a230-46181c7e64db
flow2d_forced(duffing_forced!,[pduff[5],0.0,0.0],pduff,2*pi/pduff[4];tcycles=0,ncycles=pduff[6],xlims=(-2,2),ylims=(-1.5,1.5))

# ╔═╡ 4fe120b9-f5f6-4dcf-a560-535ec98b8f72
md"""
# Poincaré section

A very useful representation for quasi-periodic periodic flows (and as we will see later also chaotic) is the Poincaré section (or Poincaré map). For three-dimensional flows it consists of taking a plane that is transverse to the flow (i.e. no orbit is parallel to it) and taking the intersection of the trajectories with this plane as points $x_i$.

We can then study the dynamics of these points on the plane of the Poincaré section as a map (hence the name Poincaré map): $x_{i+1}=P(x_i)$
where the function $P$ is obtained by calculating the trajectory from the point given by the intersection of the flow over the plane $x_i$ 
 to the next intersection with the plane $x_{i+1}$
 (see figure)
"""

# ╔═╡ c30e6a34-bb1c-43a6-b5b4-87fecf3c4211
html"""
<div>
<img src="https://i.imgur.com/raoe46x.gif" width="250px">
</div>
"""

# ╔═╡ de2a21bb-e875-4dab-b3f2-c40dfa02dfb3
md"""
Then a periodic orbit in the original flow corresponds to a fixed point in the map, if the orbit has a period equal to the time of return to the plane, or to periodic points (which return after an integer number of iterations to the same point) if the closed orbit crosses the plane several times.

For the case of forced systems the choice of the Poincaré section is obvious because since the third dimension is periodic it is sufficient to choose a plane of constant phase. Without loss of generality we can choose the ground $\phi=0$. In this case the points on the Poincaré section correspond to the points on the previous graph.
"""

# ╔═╡ 7addbab2-b219-4a77-b7a1-56af59e028be
poincare_forced(duffing_forced!,[0.5,0.5,0],pduff,2*pi/pduff[4]; tcycles=10,ncycles=100,size=(900,600))

# ╔═╡ fbd1ace7-7d18-4e9b-aa32-4003fd5cd542
md"""
# Parameter exploration 

As you can see by playing a little with the parameters, the repertoire of behaviors is highly varied. A good way to explore it is to start from the unforced system $A=0$ with and a fixed value of the other parameters (remember the effect that $\gamma$ and $\beta$ had on the original system). For example, with we have the double well potential. In that case remember that we had a saddle type fixed point at the origin and symmetric attractor foci at $x^*_{2,3}=\pm\sqrt{\beta}$. 

For small $A$ values the attractors are transformed into *stable limit cycles* with the period of the forcing. For intermediate $A$ values other attractor cycles appear (through a saddle-node bifurcation of limit cycles) surrounding the smaller limit cycle. Stable limit cycles may also appear surrounding both wells and with different period values as integer multiplies of the forcing.

On the other hand, the saddle point of the unforced system becomes a saddle-type limit cycle, i.e. it is not an attractor (except in a certain direction), but as in 2D systems it organizes the flow. Later we will see that the stable and unstable manifolds of the saddle orbit (which are actually like folded sheets in 3D) are important to understand the organization of the flow and the occurrence of chaos.

The coexistence of these attractor limit cycles already makes the behavior more complex and unpredictable (try varying the initial condition along the horizontal axis for example for the values ($\gamma=0.1,\beta=1,A=0.1,\omega=1$) to see the different possible destinations depending on the starting point or $A=0.24$ for a more complex behavior).

However for certain values of parameters and initial conditions the orbit never seems to converge to an attractor and is left spinning following an irregular path (e.g. for $\gamma=0.1,\beta=1,A=0.1,\omega=1$ ). We will see that this is a new type of attractor (called "strange" and having fractal structure) and that it is a signature of chaos.

But first let's explore the coexistence of attractors and define the basin of attraction of an attractor. If we explore different initial conditions in the case suggested above with values of $A=0.1$ and $A=0.24$:
"""

# ╔═╡ 8b29c214-4b49-4f44-9eb5-230899ec7321
md"""
# Attraction Basins

An important notion that will allow us to characterize the growing complexity of the behavior of this system is that of the basin of attraction.

The basin of attraction of a certain attractor (for example a limit cycle that in the Poincare section corresponds to a fixed point or a set of periodic points), is defined by all those initial conditions that converge to the attractor for long times.

To determine the basins it is necessary to evolve a grid of several thousand initial conditions over several cycles, therefore the graphs that follow can be very demanding. Run first with a large value of delta (the grid resolution) and make sure that Julia is taking in Threads.nthreads() the full number of processors in order to parallelize the problem.
"""

# ╔═╡ 0a270882-fde3-4dcc-8cef-291e86a4f8c2
# ╠═╡ disabled = true
#=╠═╡
md"""
Select one: $(@bind sel1 Select([([0.14,1.0,0.1,1.0],[[1.0,0.0],[-1.0,0.0],[-1,0.7],[0.3,0.3]])=>"γ=0.14,β=1,A=0.1,ω=1",([0.14,1.0,0.14,1.0],[[1.1,0.0],[-0.8,0.0]])=>"γ=0.14,β=1,A=0.14,ω=1",([0.14,1.0,0.2,1.0],[[1.2,0.0],[-0.9,0.0]])=>"γ=0.14,β=1,A=0.2,ω=1",([0.14,1.0,0.24,1.0],[[1.2,0.0],[-0.9,0.0],[0.1,1.1]])=>"γ=0.14,β=1,A=0.24,ω=1"]))
"""
  ╠═╡ =#

# ╔═╡ dc6e3382-7239-48f1-a88d-bc12fda7cb73
# ╠═╡ disabled = true
#=╠═╡
attractor_basin(duffing_forced!,sel1[1],sel1[2],0.3;state_dimension=3,delta=0.01,tmax=30*pi,xlims=(-2.5,2.5),ylims=(-2.0,2.0))
  ╠═╡ =#

# ╔═╡ 3078878b-a2f8-40f5-9398-401682ebb2f5
md"""
# Strange Attractor

For greater values of $A$, the trajectories no longer converge to limit cycles but istead they approach to a set with a fractal structure known as a **strange attractor**. On the Poincare section is a set of points that, unlike the one generated by a torus formed by quasiperiodic orbits, is not confined to a curve. Let's see the poincare section for a value $A=0.27$ that gives rise to a strange attractor:
"""

# ╔═╡ 64e04ff2-6662-4775-93cd-85a915487478
# ╠═╡ disabled = true
#=╠═╡
poincare_forced(duffing_forced!,[0.5,0.5,0],[0.14,1.0,0.27,1.0],2*pi; tcycles=30,ncycles=50000,size=(900,600))
  ╠═╡ =#

# ╔═╡ 19d2efe0-f5bc-4353-af5e-e16e424edd38
md"""
Do not confuse this structure with that of the basins of attraction. These points are the limit set (which in the case of cycles correspond to a single point) and are formed by infinite points that form a fractal structure. In this case the basin of attraction occupies the whole plane but for other values of can coexist with periodical orbits.

If we zoom in to a detail of the attractor we see that it has a structure with even more detail:
"""

# ╔═╡ bb2042f9-84bf-452e-be4b-b60b8550e787
# ╠═╡ disabled = true
#=╠═╡
poincare_forced_zoom(duffing_forced!,[0.5,0.5,0],[0.14,1,0.27,1],2*pi;npts=30000,maxiter=1000,size=(900,600),xlims=[-0.75,-0.5],ylims=[0,0.2])
  ╠═╡ =#

# ╔═╡ 6e2672fa-670c-4a6a-ac1d-3aa897a6210f
md"""
# Unstable and Stable Manifolds of the Saddle in the Poincare Section

To understand how chaos develops (we will give a more rigorous definition later) it is useful to study the stable and unstable manifolds of the saddle orbit (or the saddle point of the map). 
"""

# ╔═╡ 7b41aa57-598c-4f87-ba26-ec3c54a12024
# ╠═╡ disabled = true
#=╠═╡
md"""
Select one: $(@bind sel2 Select([([0.14,1.0,0.1,1.0],[-0.05,0.0],1500)=>"γ=0.14,β=1,A=0.1,ω=1",([0.14,1.0,0.14,1.0],[-0.07,0.0],2500)=>"γ=0.14,β=1,A=0.14,ω=1",([0.14,1.0,0.24,1.0],[-0.1,0.0],5500)=>"γ=0.14,β=1,A=0.24,ω=1",([0.14,1.0,0.27,1.0],[-0.1,0.0],5500)=>"γ=0.14,β=1,A=0.27,ω=1"]))
"""
  ╠═╡ =#

# ╔═╡ 788140f8-9705-4784-8089-4bb4c942c3a8
# ╠═╡ disabled = true
#=╠═╡
begin
	p = sel2[1]
	period = 2*pi
	npts = 1500
	us,conv = saddle_orbit2D(duffing_forced!,sel2[2],p,period)
	if conv
	    #println(us)
	    p1=saddle_manifolds_forced(duffing_forced!,duffing_jac,us,p,period;ncycles=[8,2],npts=npts,delta=1e-4)
		if (p[3]<0.27)
	    	poincare_forced!(p1,duffing_forced!,[1.2,0.0,0.0],p, period; tcycles=40,ncycles=41,msize=3.0,col=:black)
	    	poincare_forced!(p1,duffing_forced!,[-0.8,0.0,0.0],p, period; tcycles=40,ncycles=41,msize=3.0,col=:black)
			if (p[3]<0.14)
	    		poincare_forced!(p1,duffing_forced!,[-1,0.7,0.0],p, period; tcycles=40,ncycles=41,msize=3.0,col=:black)
	    		poincare_forced!(p1,duffing_forced!,[0.3,0.3,0.0],p, period; tcycles=40,ncycles=41,msize=3.0,col=:black)
			end
		end	
	    xlims!(p1,(-2.5,2.5)); ylims!(p1,(-1.5,1.5))
		title!(p1,"Forced Duffing Oscillator. γ=0.14, β=1, A=0.14, ω=1")
	end	
end	
  ╠═╡ =#

# ╔═╡ 5e7d4da9-1ba7-44fd-b5ca-8931fa6dc00c
# ╠═╡ disabled = true
#=╠═╡
function duffing_jac(u,p)
  J = Array{Float64, 2}(undef, 3, 3)
  J[1,1] = 0
  J[1,2] = 1.0
  J[1,3] = 0.0  
  J[2,1] = p[2]-3.0*u[1]*u[1]
  J[2,2] = -p[1]  
  J[2,3] = -p[3]*sin(u[3]) 
  J[3,1] = 0.0
  J[3,2] = 0.0
  J[3,3] = 0.0  
  return J
end;
  ╠═╡ =#

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

# ╔═╡ d4847686-80fe-539b-9bcc-bee2dec3d50a
md"""
## Compare a time trace, a stroboscopic section, and recurrence

Sample at $t_n=nT$, with $T=2\pi/\omega$, after discarding a transient.
One point suggests forcing-period locking; a finite set suggests a multiple
of that period. A smooth closed curve is consistent with quasiperiodicity.
A more complicated set calls for additional checks before claiming chaos.

The same equations above are used for the comparison below, with
$\gamma=0.14$, $\beta=1$, $\omega=1$, and adjustable forcing amplitude.
Try $A=0.1$ and $A=0.27$. The transient and the displayed segment are solved
on one continuous time interval, preserving the forcing phase.
"""

# ╔═╡ b444e370-6014-524f-987c-9b7fe8294396
md"""
forcing amplitude A: $(@bind course_forcing_A Slider(0.0:0.01:0.35; default=0.27, show_value=true))
"""

# ╔═╡ b59600a9-9720-52b2-acd9-9ed33517a22f
course_forcing_u0 = [0.5,0.5,0.0]

# ╔═╡ feada8e0-3663-5d83-9f50-1c25360ca23f
course_forcing_p = [0.14,1.0,course_forcing_A,1.0]

# ╔═╡ 33cb060e-5b6a-5ed1-9045-38910522837f
course_forcing_period = 2pi/course_forcing_p[4]

# ╔═╡ 0e81c277-2397-5bb1-9b18-07427a025f92
course_forcing_tspan = (0.0,160course_forcing_period)

# ╔═╡ 9940cae4-90a6-599b-8f69-6e0dceb402b7
course_forcing_prob = ODEProblem(duffing_forced!, course_forcing_u0, course_forcing_tspan, course_forcing_p)

# ╔═╡ 8ee3834c-fb55-5b5e-a1ee-c693cd9ed70b
course_forcing_sol = solve(course_forcing_prob, Tsit5(); abstol=1e-9, reltol=1e-8, dtmax=course_forcing_period/40);

# ╔═╡ 767a781c-4cb6-5756-8dff-5311b956cd98
let period=course_forcing_period
    section_times=(60:160).*period
    points=Array(course_forcing_sol(section_times;idxs=[1,2]))
    trace=plot(course_forcing_sol;idxs=(0,1),tspan=(140period,160period),
        xlabel="t",ylabel="x",legend=false)
    section=scatter(points[1,:],points[2,:];xlabel="x(nT)",ylabel="v(nT)",
        markersize=2,markerstrokewidth=0,legend=false)
    plot(trace,section;layout=(1,2),size=(950,380),margin=5*Plots.mm)
end

# ╔═╡ 116e0262-55b9-5d5c-b9ba-fec9a998321e
md"""
A recurrence plot asks when the system returns near a previous **full state**:

```math
R_{ij}=\begin{cases}1&\|\mathbf z_i-\mathbf z_j\|\leq\varepsilon,\\0&\text{otherwise},
\end{cases}
\qquad
\mathbf z(t)=(x(t),v(t),\cos\omega t,\sin\omega t).
```

The two phase coordinates identify the forcing state without a discontinuity
at $2\pi$. Comparing only $(x,v)$ could mark a return at a different forcing
phase. All coordinates here are dimensionless with unit scales; changing
scales or the threshold changes the notion of “near”.

Repeated diagonal bands indicate repeated evolution; broken bands indicate
less regular recurrence. Neither a recurrence image nor an irregular time
trace alone establishes chaos.
"""

# ╔═╡ 475ac4f6-70bb-5462-817e-c732e4caf754
md"""
recurrence threshold ε: $(@bind course_epsilon Slider(0.05:0.05:0.5; default=0.2, show_value=true))
"""

# ╔═╡ b3f24b1e-485e-5f9b-ab7a-927f7f9a7a6e
course_recurrence_times = collect(range(120course_forcing_period,140course_forcing_period;length=401))

# ╔═╡ 71d510c8-598e-5129-8e4c-be2869b510c0
course_recurrence_states = vcat(
    Array(course_forcing_sol(course_recurrence_times;idxs=[1,2])),
    permutedims(cos.(course_forcing_p[4].*course_recurrence_times)),
    permutedims(sin.(course_forcing_p[4].*course_recurrence_times)))

# ╔═╡ 3413d3f9-8a8d-5d73-90bc-126805aa5627
let states=course_recurrence_states, times=course_recurrence_times
    # Each column is [x, v, cos(phase), sin(phase)] at one sample time.
    n = size(states, 2)
    R = [sum(abs2, states[:,i] - states[:,j]) <= course_epsilon^2
         for i in 1:n, j in 1:n]
    heatmap(times, times, Int.(R); color=cgrad([:white, :black]), clims=(0,1),
        colorbar=false, xlabel="tᵢ", ylabel="tⱼ", title="Recurrence: distance ≤ ε",
        aspect_ratio=1, size=(550,500), margin=5*Plots.mm)
end

# ╔═╡ d8ccbb3c-be01-5cf5-9609-6ce0898074ca
md"""
**Try:** compare the two amplitudes with the same threshold, sampling window,
and discarded transient. Then change only $\varepsilon$. Which features
persist? Return to NLD01 and explain why a periodic orbit of a flow can
be studied as a fixed or periodic point of a map.
"""

# ╔═╡ Cell order:
# ╠═8601d8d7-d4df-473f-b65d-0f03aeb8f5f4
# ╟─1f49c325-fc99-4b15-81ee-dc1c5bbe6f08
# ╟─6b52504e-805e-417e-8fb2-6473b47b0ad0
# ╠═6b14f6c6-536b-4de2-b669-e0f1f34dbbe2
# ╠═37ddd39f-4de1-4eb0-a230-46181c7e64db
# ╟─4ee77c37-0347-45c0-b101-b96fcb7e004d
# ╟─4fe120b9-f5f6-4dcf-a560-535ec98b8f72
# ╟─c30e6a34-bb1c-43a6-b5b4-87fecf3c4211
# ╟─de2a21bb-e875-4dab-b3f2-c40dfa02dfb3
# ╟─7addbab2-b219-4a77-b7a1-56af59e028be
# ╟─fbd1ace7-7d18-4e9b-aa32-4003fd5cd542
# ╟─8b29c214-4b49-4f44-9eb5-230899ec7321
# ╟─0a270882-fde3-4dcc-8cef-291e86a4f8c2
# ╟─dc6e3382-7239-48f1-a88d-bc12fda7cb73
# ╟─3078878b-a2f8-40f5-9398-401682ebb2f5
# ╟─64e04ff2-6662-4775-93cd-85a915487478
# ╟─19d2efe0-f5bc-4353-af5e-e16e424edd38
# ╠═bb2042f9-84bf-452e-be4b-b60b8550e787
# ╟─6e2672fa-670c-4a6a-ac1d-3aa897a6210f
# ╟─7b41aa57-598c-4f87-ba26-ec3c54a12024
# ╟─788140f8-9705-4784-8089-4bb4c942c3a8
# ╟─5e7d4da9-1ba7-44fd-b5ca-8931fa6dc00c
# ╟─c7399360-fad1-466a-9324-f830b20b07a7
# ╟─1f7f958e-b0c3-433f-bae9-3bf63da3de7a
# ╟─d4847686-80fe-539b-9bcc-bee2dec3d50a
# ╟─b444e370-6014-524f-987c-9b7fe8294396
# ╠═b59600a9-9720-52b2-acd9-9ed33517a22f
# ╠═feada8e0-3663-5d83-9f50-1c25360ca23f
# ╠═33cb060e-5b6a-5ed1-9045-38910522837f
# ╠═0e81c277-2397-5bb1-9b18-07427a025f92
# ╠═9940cae4-90a6-599b-8f69-6e0dceb402b7
# ╠═8ee3834c-fb55-5b5e-a1ee-c693cd9ed70b
# ╠═767a781c-4cb6-5756-8dff-5311b956cd98
# ╟─116e0262-55b9-5d5c-b9ba-fec9a998321e
# ╟─475ac4f6-70bb-5462-817e-c732e4caf754
# ╠═b3f24b1e-485e-5f9b-ab7a-927f7f9a7a6e
# ╠═71d510c8-598e-5129-8e4c-be2869b510c0
# ╠═3413d3f9-8a8d-5d73-90bc-126805aa5627
# ╟─d8ccbb3c-be01-5cf5-9609-6ce0898074ca
