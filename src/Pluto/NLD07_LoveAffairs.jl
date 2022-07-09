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

# ╔═╡ 673fa0ec-754a-41a7-884f-d153f411e41c
using Pkg;Pkg.add("LinearAlgebra")

# ╔═╡ f527ceca-f2d1-11ec-3cc5-bd33fdb53d6e
using DifferentialEquations, Plots, ForwardDiff,IntervalRootFinding, StaticArrays, PlutoUI, LinearAlgebra

# ╔═╡ cf2a4cec-ea32-4bd1-a56d-8fde23d3e31c
include("../NLD_utils.jl")

# ╔═╡ 54dd97fe-4ae1-4f29-a30f-b2d11f8cbf62
TableOfContents()

# ╔═╡ d777dcf0-600c-4bce-b6cf-da7e57b1c1ae
md"""
# Nonlinear Dynamics of Love Affairs 💘

# Classic Linear model of Romantic Love 👫

Based on the one-page influential work of Strogatz (1988) of the Shaskespearean love affair of Romeo and Juliet. Here $R$ is Romeo’s love (or hate if negative) for Juliet and $J$ is Juliet’s love for Romeo. The simplest model is linear with:

$\dot{R}= aR+bJ$

$\dot{J}= cR+dJ$

where $a$ and $b$ specify Romeo’s “romantic style,” and $c$ and $d$ specify Juliet’s style. The parameter $a$($c$) describes the extent to which Romeo (Juliet) is encouraged by his(her) own feelings, and $b$($d$) is the extent to which he is encouraged by his(her) partner feelings.

The resulting dynamics is
two-dimensional, governed by the initial conditions and the four
parameters, which may be positive or negative.
"""


# ╔═╡ 42638d65-5778-4776-b7f6-20550c515289
function lovelinear!(du,u,p,t)
	(a,b,c,d) = p
	du[1] = a*u[1]+b*u[2]
	du[2] = c*u[1]+d*u[2]
end

# ╔═╡ 0f57818a-9608-4f24-84a0-b5579e72fc2d
md"""
Romeo can exhibit one of four romantic styles depending on the signs of a and b, with names adapted from those suggested by Strogatz and his students:

1. Eager beaver: $a > 0$, $b > 0$ (Romeo is encouraged by his own feelings as well as Juliet’s.)
2. Narcissistic nerd: $a > 0$, $b < 0$ (Romeo wants more of what he feels but retreats from Juliet’s feelings.)
3. Cautious (or secure) lover: $a < 0$, $b > 0$ (Romeo retreats from his own feelings but is encouraged by Juliet’s.)
4. Hermit: $a < 0$, $b < 0$ (Romeo retreats from his own feelings as well as Juliet’s.) 

The same four styles applies for Juliet following the sign of the parametes $c$ and $d$
"""	

# ╔═╡ b35cce8b-796c-4f5c-b312-68d7084d5e90
md"""
Partner Name 1 $(@bind name1 TextField((20,1);default="Juliet"))
$(@bind heshe1 Select(["she", "he"])) \
Partner Name 1 $(@bind name2 TextField((20,1);default="Romeo")) 
$(@bind heshe2 Select(["he", "she"])) 
"""

# ╔═╡ dc301dcf-2b27-431d-8604-1c94f6c6d353
@bind pars1 (
	PlutoUI.combine() do bind
		md"""
		a: $(bind(Slider(-1.0:0.02:1.0,default=0.1;show_value=true))) 
		b: $(bind(Slider(-1.0:0.02:1.0,default=0.1;show_value=true))) \
		c: $(bind(Slider(-1.0:0.02:1.0,default=0.1;show_value=true))) 
		d: $(bind(Slider(-1.0:0.02:1.0,default=0.1;show_value=true))) \
		"""
	end
)

# ╔═╡ 59d4e977-7d25-4ef1-98de-9a2a712dfa6c
@bind u1 (
	PlutoUI.combine() do bind
		md"""
		x$0$: $(bind(Slider(-1.0:0.02:1.0,default=0.1;show_value=true))) 
		y0: $(bind(Slider(-1.0:0.02:1.0,default=0.1;show_value=true))) \
		"""
	end
)

# ╔═╡ ad546f0a-fdfa-422c-8ed5-3e2b59f9cbfb
md"""

Although Strogatz’s model was originally intended more to motivate students than
as a serious description of love affairs, it makes several interesting and
plausible predictions and suggests extensions that produce an even wider
range of behavior.


"""

# ╔═╡ d2988970-f026-4eb3-8ec7-b92efc4e0b4a
md"""
# Classification of the combination of loving styles (or the stability of the fixed point)

Let us now see how we can classify this pair of interacting behaviors, and in general of this linear system. It is a linear system because it only depends on the variables in the form of a constant per variable and there are no quadratic or cubic terms or any other functional form.

As can be seen also in the 2D flows, linear systems can only have one fixed point.

To study the stability of this fixed point we are going to evolve a set of initial conditions, as if we were dealing with two fixed loving styles but which combine under different circumstances. 


When studying the stability of the fixed point, we will see that the description is a little more subtle than in the 1D case where we could only have an attractor or a repeller (or at most a neutral point as a marginal case).

"""

# ╔═╡ 706f8d27-9357-4e17-a288-867c08969fa2
@bind pars2 (
	PlutoUI.combine() do bind
		md"""
		a: $(bind(Slider(-1.0:0.02:1.0,default=0.1;show_value=true))) 
		b: $(bind(Slider(-1.0:0.02:1.0,default=0.1;show_value=true))) \
		c: $(bind(Slider(-1.0:0.02:1.0,default=0.1;show_value=true))) 
		d: $(bind(Slider(-1.0:0.02:1.0,default=0.1;show_value=true))) \
		"""
	end
)

# ╔═╡ fc617669-c8ed-47fa-a430-db6299b5cfde
begin 
	#u0_array = vec([[cos(i*pi/10);sin(i*pi/10)] for i=1:20, j=1:20]);
	u0_array = vec([[-1+i*0.2;-1+j*0.2] for i=1:10, j=1:10]);
    condition(u,t,integrator) = (u[1]*u[1]+u[2]*u[2]) > 1
    affect!(integrator) = terminate!(integrator)
	prob = ODEProblem(lovelinear!,u0_array[1],(0.0,2.0),pars2)
    ensamble_prob = EnsembleProblem(prob,prob_func=(prob,i,repeat;u0=u0_array)->(remake(prob,u0=u0[i])))
    sol = solve(ensamble_prob,EnsembleThreads(),trajectories=length(u0_array),
		callback=DiscreteCallback(condition,affect!))
    plot(sol,vars=(1,2),c=:black,arrow=true,xlims=(-1,1),ylims=(-1,1),linewidth=0.5,labels="",xlabel="x",ylabel="y",size = (400,400))
end

# ╔═╡ c4da9d4d-6035-4792-9e58-5006bd566441
begin
	function node1!(du,u,p,t)
	    du[1] = -u[1]
	    du[2] = -2.0*u[2]
	end   
	u0_arr = vec([[cos(i*pi/10);sin(i*pi/10)] for i=1:20])
	prob1 = ODEProblem(node1!, u0_arr[1], (0,3.0))
	function prob_func(prob,i,repeat)
	  remake(prob,u0=u0_arr[i])
	end
	sol1 = solve(EnsembleProblem(prob1,prob_func=prob_func),
	    EnsembleThreads(),trajectories=length(u0_arr))
	pn1 = plot(sol1,vars=(1,2),arrows=(:head,1.0),legend=false,xlabel="x",ylabel="y",
	    title="Attractor Node")
	pn2 = plot(sol1,vars=(1,2),arrows=(:tail,1.0),legend=false,xlabel="x",ylabel="y",
	    title="Repeller Node")
	plot(pn1,pn2,layout=(1,2),size=(900,450))
end	

# ╔═╡ 21b1ac8b-e2a0-416a-a84c-227152546c3e
begin
	function spiral1!(du,u,p,t)
	    du[1] = 2*u[2]
	    du[2] = -1*u[1]-2*u[2]
	    return 
	end    
	
	#u0_arr = vec([[cos(i*pi/10);sin(i*pi/10)] for i=1:20, j=1:20])
	prob3 = ODEProblem(spiral1!, u0_arr[1], (0,5.0))
	sol3 = solve(EnsembleProblem(prob3,prob_func=prob_func),
	    EnsembleThreads(),trajectories=length(u0_arr))
	pf1 = plot(sol3,vars=(1,2),arrows=(:head,1.0),legend=false,xlabel="x",ylabel="y",
	    title="Attractor Focus")
	pf2 = plot(sol3,vars=(1,2),arrows=(:tail,1.0),legend=false,xlabel="x",ylabel="y",
	    title="Repeller Focus")
	plot(pf1,pf2,layout=(1,2),size=(900,450))
end	

# ╔═╡ a40cee7a-994e-4860-9c97-394d7dbfec77
begin
	function saddle1!(du,u,p,t)
	    du[1] = u[1]
	    du[2] = -2.0*u[2]
	end    
	prob2 = ODEProblem(saddle1!, u0_arr[1], (0,1.0))
	sol2 = solve(EnsembleProblem(prob2,prob_func=prob_func),
	    EnsembleThreads(),trajectories=length(u0_arr))
	plot(sol2,vars=(1,2),arrows=true,legend=false,xlabel="x",ylabel="y",
	    size=(400,400),fmt=:png,title="Saddle Point")
end	

# ╔═╡ 4aa62401-a1c6-47f8-807f-7444d542fbde
@bind parscl (
	PlutoUI.combine() do bind
		md"""
		Trace: $(bind(Slider(-1.0:0.02:1.0,default=0.1;show_value=true))) 
		Determinant: $(bind(Slider(-1.0:0.02:1.0,default=0.1;show_value=true))) \
		d: $(bind(Slider(-1.0:0.02:1.0,default=0.1;show_value=true))) 
		"""
	end
)

# ╔═╡ b4bf8b41-dff9-4e87-94c8-bda9365bd78e
begin 
	(trace,determ,d1)=parscl
    b1=1.0;
    a1=trace-d1;
    c1=(a1*d1-determ)/b1
    A = round.([a1 b1; c1 d1],digits=2)
    classification_linear(A)
end	

# ╔═╡ 17da3943-629d-49cf-88b4-febc97bb26bb
md"""
# Nonlinear (Non-romantic?) Love

"""

# ╔═╡ d2008fc9-dd22-4725-9355-dacb7ff5821b
function lovenlinear1!(du,u,p,t)
	(a,b,c,d,ϵ) = p
	du[1] = a*u[1]+b*u[2]*(1-ϵ*u[2]*u[2])
	du[2] = d*u[2]+c*u[1]*(1-ϵ*u[1]*u[1])
end

# ╔═╡ 40ae473c-3362-468d-a988-0e0c1d85b5ea
@bind pars3 (
	PlutoUI.combine() do bind
		md"""
		a: $(bind(Slider(-1.0:0.02:1.0,default=0.1;show_value=true))) 
		b: $(bind(Slider(-1.0:0.02:1.0,default=0.1;show_value=true))) \
		c: $(bind(Slider(-1.0:0.02:1.0,default=0.1;show_value=true))) 
		d: $(bind(Slider(-1.0:0.02:1.0,default=0.1;show_value=true))) \
		ϵ: $(bind(Slider(-0.5:0.001:0.5,default=0.0;show_value=true))) 
		tmax : $(bind(Slider(0:10.0:1000.0,default=100.0;show_value=true)))
		"""
	end
)

# ╔═╡ d622eafe-f2d5-4992-ac3f-fabb392bb12a
# Classification of Romantic Styles for R&J
begin
	(a2,b2,c2,d2,ϵ,_)=pars3
	if (a2>0)
		if (b2>0)
			rtitle2 = "Romeo: Eager Beaver"
		else
			rtitle2 = "Romeo: Narcissistic"
		end
	else
		if (b2>0)
			rtitle2 = "Romeo: Cautious"
		else
			rtitle2 = "Romeo: Hermit"
		end
	end
	if (c2>0)
		if (d2>0)
			jtitle2 = "Juliet: Eager Beaver"
		else
			jtitle2 = "Juliet: Narcissistic"
		end
	else
		if (d2>0)
			jtitle2 = "Juliet: Cautious"
		else
			jtitle2 = "Juliet: Hermit"
		end
	end
	if (ϵ>0)
		title2 = string(rtitle2," | ",jtitle2," | Positive Nonlinearity")
	elseif (ϵ<0)
		title2 = string(rtitle2," | ",jtitle2," | Negative Nonlinearity")
	else
		title2 = string(rtitle2," | ",jtitle2," | Zero Nonlinearity")
	end	
		
end

# ╔═╡ 7862b57f-93c1-4293-8e83-57b7c18dc3fc
@bind u2 (
	PlutoUI.combine() do bind
		md"""
		x$0$: $(bind(Slider(-1.0:0.02:1.0,default=0.1;show_value=true))) 
		y0: $(bind(Slider(-1.0:0.02:1.0,default=0.1;show_value=true))) \
		"""
	end
)

# ╔═╡ 86eb471e-cedf-4e5e-bcac-84476eb21274
function flow2d(f::Function,u0::Vector{Float64},tmax::Float64,p;
    xlims=[-1.0,1.0],ylims=[-1.0,1.0],size=(700,500),plotops...)
    xrange = xlims[2]-xlims[1]
    yrange = ylims[2]-ylims[1]
    condition(u,t,integrator) = (u[1]*u[1]+u[2]*u[2]) > max(xrange*xrange,yrange*yrange)
    affect!(integrator) = terminate!(integrator)
    sol = solve(ODEProblem(f,u0,(0.0,tmax),p),callback=DiscreteCallback(condition,affect!))
    p1 = plot(sol,vars=(0,1:2))
    p2 = plot(sol,vars=(1,2),c=:black,arrow=true,xlims=xlims,ylims=ylims,labels="",xlabel="x",ylabel="y")
    plot(p1,p2,layout=(1,2),size = size)
end;

# ╔═╡ 5db59da7-984a-4f90-9a90-6a35cb61778e
flow2d(lovelinear!,collect(u1),200.0,pars1;size=(800,400))

# ╔═╡ a0662db2-84d5-473c-88a7-19066f8f9de6
flow2d(lovenlinear1!,collect(u2),pars3[6],pars3;size=(800,400),xlims=[-5,5],ylims=[-5,5])

# ╔═╡ 211a4a3a-640a-492c-a7b8-e1f0eae4d7f4
# Classification of Romantic Styles 
function romanticlass(pars,name1,name2,heshe1,heshe2)
	(a,b,c,d)=pars
	if heshe1=="he"
		hisher1="his"
	else
		hisher1="her"
	end
	if heshe2=="he"
		hisher2="his"
	else
		hisher2="her"
	end
	if (a>0)
		if (b>0)
			rtitle = "$name1 (EAGER) is encouraged by $hisher1 own feelings as well as $(name2)’s"
		else
			rtitle = "$name1 (NARCISSISTIC) wants more of what $heshe1 feels but retreats from $(name2)’s feelings"
		end
	else
		if (b>0)
			rtitle = "$name1 (CAUTIOUS) Retreats from $hisher1 own feelings but is encouraged by $(name2)'s"
		else
			rtitle = "$name1 (HERMIT) retreats from $hisher1 own feelings as well as $(name2)’s"
		end
	end
	if (c>0)
		if (d>0)
			jtitle = "$name2 (EAGER) is encouraged by $hisher2 own feelings as well as $(name1)’s"
		else
			jtitle = "$name2 (NARCISSISTIC) wants more of what $heshe2 feels but retreats from $(name1)’s feelings"
		end
	else
		if (d>0)
			jtitle = "$name2 (CAUTIOUS) Retreats from $hisher2 own feelings but is encouraged by $(name1)'s"
		else
			jtitle = "$name1 (HERMIT) retreats from $hisher1 own feelings as well as $(name2)’s"
		end
	end
	return rtitle , jtitle
end;

# ╔═╡ ca9f485c-585a-4e53-899e-b908d9dd7963
begin
	tit1, tit2 = romanticlass(pars1,name1,name2,heshe1,heshe2)
	println(tit1)
	println(tit2)
end	

# ╔═╡ e9a8a3f7-065f-4a51-a688-433765742a09
begin
	tit1b, tit2b = romanticlass(pars2,name1,name2,heshe1,heshe2)
	println(tit1b)
	println(tit2b)
end	

# ╔═╡ Cell order:
# ╠═673fa0ec-754a-41a7-884f-d153f411e41c
# ╠═f527ceca-f2d1-11ec-3cc5-bd33fdb53d6e
# ╠═cf2a4cec-ea32-4bd1-a56d-8fde23d3e31c
# ╟─54dd97fe-4ae1-4f29-a30f-b2d11f8cbf62
# ╟─d777dcf0-600c-4bce-b6cf-da7e57b1c1ae
# ╠═42638d65-5778-4776-b7f6-20550c515289
# ╟─0f57818a-9608-4f24-84a0-b5579e72fc2d
# ╟─b35cce8b-796c-4f5c-b312-68d7084d5e90
# ╟─ca9f485c-585a-4e53-899e-b908d9dd7963
# ╟─5db59da7-984a-4f90-9a90-6a35cb61778e
# ╟─dc301dcf-2b27-431d-8604-1c94f6c6d353
# ╟─59d4e977-7d25-4ef1-98de-9a2a712dfa6c
# ╟─ad546f0a-fdfa-422c-8ed5-3e2b59f9cbfb
# ╟─d2988970-f026-4eb3-8ec7-b92efc4e0b4a
# ╟─e9a8a3f7-065f-4a51-a688-433765742a09
# ╟─fc617669-c8ed-47fa-a430-db6299b5cfde
# ╟─706f8d27-9357-4e17-a288-867c08969fa2
# ╟─c4da9d4d-6035-4792-9e58-5006bd566441
# ╟─21b1ac8b-e2a0-416a-a84c-227152546c3e
# ╟─a40cee7a-994e-4860-9c97-394d7dbfec77
# ╟─b4bf8b41-dff9-4e87-94c8-bda9365bd78e
# ╟─4aa62401-a1c6-47f8-807f-7444d542fbde
# ╟─17da3943-629d-49cf-88b4-febc97bb26bb
# ╠═d2008fc9-dd22-4725-9355-dacb7ff5821b
# ╟─d622eafe-f2d5-4992-ac3f-fabb392bb12a
# ╠═a0662db2-84d5-473c-88a7-19066f8f9de6
# ╟─40ae473c-3362-468d-a988-0e0c1d85b5ea
# ╟─7862b57f-93c1-4293-8e83-57b7c18dc3fc
# ╟─86eb471e-cedf-4e5e-bcac-84476eb21274
# ╟─211a4a3a-640a-492c-a7b8-e1f0eae4d7f4
