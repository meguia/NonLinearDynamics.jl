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

# ╔═╡ 38c4e220-d708-11ec-3968-fbbd41c26155
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using NonLinearDynamics: flow1D
    using PlutoUI, Plots, DifferentialEquations, ForwardDiff, IntervalRootFinding, StaticArrays
end

# ╔═╡ e2271efe-9a67-4496-891b-46f076b367f3
theme(:default)

# ╔═╡ 62cba86d-4406-4d16-8715-b29bee7d3178
md"""
# Logistic Equation

In what follows we will work with simple population dynamics models where the continuous variable can represent the density of a population (i.e. both the variable and time evolve continuously). In this context it is usual to denote the growth rate as $R$ (rate). Therefore the equation:


$\dot{x} = Rx$

corresponds to the unlimited exponential growth of the population density at rate $R$. 

But exponential growth cannot be maintained forever, so if we want to model in a more realistic way a magnitude that grows with limited resources we must limit the growth rate $R$. One possible way (the simplest) is to limit $R$ with a linear function becoming zero when the population reaches the maximum capacity $K$. That is, let our variable growth rate be: $R(x)=R(1-x/K)$. Replacing this rate in the original equation we obtain the logistic equation 

$\dot{x}=Rx\left(1-\frac{x}{K}\right)$

wuth $K>0$ and $R>0$.

This is the simplest system that models the growth of a population with K capacity.
"""

# ╔═╡ 6d64cc27-1059-4dd7-8f17-1c6fd9a2eca0
html"""
<div style="position: relative; right: 0; top: 0; z-index: 300;"><iframe src="https://www.youtube.com/embed/JwYhhnuuINk" width=500 height=250  frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe></div>
"""

# ╔═╡ 4dc32b5a-5794-4d2d-844e-cdb80bec2e7b
logistic(x,p,t)=p[1]*x*(1.0-x/p[2]);

# ╔═╡ aa22c728-662c-4fdb-aae9-b3b411e311e2
md"""
Now we are going to plot, in addition to the time evolution, the function $f(x)$ in order to visualize the fixed points and the flow on the line.

In order to do this, we are going to use the function of the course package `flow1D` whose arguments are (for the basic method)

`flow1D(f,x0,tmax,p;xlims)`
- f is a function f(x,p,t) that defines the flow (returns the time derivative).
- x0 is a scalar which gives the initial condition
- tmax is a scalar that gives the final time of integration (always starts from t=0)
- p is a vector (1D array) with the parameters of the system. Even if we have only one parameter we have to pass it as an array of only one element
- the tuple xlims=(xmin,xmax) limits the graph of f(x) between the interval (xmin,xmax). This parameter is optional so it is after the semicolon.


We can see that $f(x)$ is an inverted parabola that cuts the horizontal axis always in the fixed points $x=0$ and $x=K$
and always with the same stability. That is, there are no topological changes (no bifurcations).
Why inverted? How would it be the other way around in the positives?
"""

# ╔═╡ 77766c99-4142-4acb-a855-ab133e67aa66
@bind pars2 (
	PlutoUI.combine() do bind
		md"""
		R: $(bind(Slider(0:0.02:2.0,default=0.5;show_value=true))) \
		K: $(bind(Slider(0.02:0.02:2.0,default=1.0;show_value=true))) \
		x0: $(bind(Slider(0:0.01:2.0,default=0.01;show_value=true)))
		"""
	end
)

# ╔═╡ 16a156b5-92ab-4e1c-ab70-c3b348186c57
flow1D(logistic,pars2[3],100.0,pars2;xlims=[-0.1,2.0],title="Logistic")

# ╔═╡ ad67a3cb-a42f-451d-92a5-301559322c30
md"""
## Math addendum 

This equation has two fixed points. One at $x=0$ and the other at $x=K$. To evaluate its stability we calculate the derivative of the vector field $f(x)=Rx(1-x/K)$, which is equal to $f'(x)=R-2Rx/K$, and evaluate it at the fixed points.

- Fixed point $x_*=0$ : $f'(0)=R$ is always positive $\rightarrow$ the fixed point is unstable (repulsive). In its neighborhood, for a small positive population there is an exponential growth with rate R since 1 is much larger than $x/K$.

- Fixed point $x_*=K$: $f'(K)=-R$ is always negative $\rightarrow$ the fixed point is stable (attractor). It represents the maximum population that reaches growth.  If $x$ is slightly less than $K$ then $(1-x/K)$ is positive, the derivative is positive and the population grows. If instead $x$ is slightly higher than $K$ then $(1-x/K)$ is negative, the derivative is negative and the population decreases until reaching equilibrium at $x=K$.
"""

# ╔═╡ ce4da9bf-316b-45d2-8504-f2a63cdda247
md"""
# Logistic Equation with Harvest

Now, we will take a look at a simple population model that traverses an attractor-repeller (or Saddle-Node in 1D) bifurcation. 

All the population models that we will see below are based on the logistic model and add an additional term that accounts for predation, either by another species or by exploitation of the species as a resource. In its simplest version, a constant decreasing term is added to represent harvesting (harvesting, hunting, fishing)

$\dot{x} = Rx\left(1-\frac{x}{K}\right) - H$ 

Warning: this is not a realistic model for a population because it may give negative x values. Because of this we must introduce a cut-off condition when the variable becomes negative: $(u)->(u<0)$.

"""

# ╔═╡ 29474ca0-839b-45de-a7bf-7d91c5ea59aa
logharvest1(x,p,t)=p[1]*x*(1.0-x/p[2])-p[3];

# ╔═╡ 8fc407f4-a242-4cdc-b7e0-062c8f5782d1
@bind pars_harvest (
	PlutoUI.combine() do bind
		md"""
		R: $(bind(Slider(0:0.02:2.0,default=0.5;show_value=true))) \
		K: $(bind(Slider(0.02:0.02:2.0,default=1.0;show_value=true))) \
		H: $(bind(Slider(0.0:0.02:2.0,default=0.1;show_value=true))) \
		x0: $(bind(Slider(0:0.02:2.0,default=0.5;show_value=true)))
		"""
	end
)

# ╔═╡ 794936ae-d9ae-4319-994f-4282c99c4238
flow1D(logharvest1,pars_harvest[4],300.0,pars_harvest,(u)->(u<0);xlims=[0.0,2.0],title="Logistic with Harvest")

# ╔═╡ b13d261e-1c1b-41b0-ae86-e20cd2304727
md"""
## Critical Slowing Down
"""

# ╔═╡ c3757398-060c-412c-997b-58331b2dee04
@bind pars_csd (
	PlutoUI.combine() do bind
		md"""
		H: $(bind(Slider(0.2:0.002:0.25,default=0.2;show_value=true))) \
		S: $(bind(Slider(0:0.001:0.1,default=0.0;show_value=true)))
		"""
	end
)

# ╔═╡ f1f04dd4-0447-4e48-9bab-f34058a6a0f1
flow1D(logharvest1,0.5+sqrt(0.25-pars_csd[1]),200.0,[1.0,1.0,pars_csd[1]],10.0,pars_csd[2],(u)->(u<0);xlims=[0.0,1.0],title="Log whith Harvest perturbed")

# ╔═╡ c803ed90-958c-4d0e-84d5-3e8cc9c8fbfe
md"""
# Consumer Equation

The next population model is a bit more realistic and is known as the **consumer equation**. It is a model that appears in macroeconomics as a simplification of the dynamics of consumption of a renewable resource and is made by the logistic equation with a harvest or consumption term that is directly proportional to the abundance of the resource $Px$. The parameter $P$ corresponds to the rate of consumption of the resource and $R$ can be reinterpreted as the rate of production or generation of the renewable resource.  


$\dot{x} = Rx\left(1-\frac{x}{K}\right) - Px$ 
"""

# ╔═╡ b4ea8364-a582-443c-bdb1-46e75a5c4d4e
# Consumer Equation
consumer(x,p,t)=p[1]*x*(1.0-x/p[2])-p[3]*x;

# ╔═╡ 5afa3aff-af15-429b-a0fb-9a0dfcad740d
@bind pars_consumer (
	PlutoUI.combine() do bind
		md"""
		R: $(bind(Slider(0:0.02:1.0,default=0.5;show_value=true))) \
		K: $(bind(Slider(0.02:0.02:2.0,default=1.0;show_value=true))) \
		P: $(bind(Slider(0.0:0.01:0.5,default=0.0;show_value=true))) \
		x0: $(bind(Slider(0:0.02:2.0,default=0.5;show_value=true)))
		"""
	end
)

# ╔═╡ 5018d342-0a51-4275-8c0f-1b5d253c2ec0
flow1D(consumer,pars_consumer[4],300.0,pars_consumer,xlims=[0.0,2.0],title="Consumer Equation")

# ╔═╡ 84a596c1-3a0a-4edd-a7ea-f3fcd8839c2a
md"""
# Logistic Equation with Outbreak

$\dot{x} = Rx\left(1-\displaystyle\frac{x}{K}\right)-P\displaystyle\frac{x^2}{1+x^2}$
"""

# ╔═╡ 91b2718a-a921-4bae-99c0-91ec9c7e6479
html"""
<div style="position: relative; right: 0; top: 0; z-index: 300;"><iframe src="https://www.youtube.com/embed/1CSKTCS6st8" width=500 height=250  frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe></div>
"""

# ╔═╡ 5a0b156e-b0b2-43b8-b15e-bcdb65c154bf
logoutbreak(x,p,t)=p[1]*x*(1.0-x/p[2])-p[3]*x*x/(1+x*x);

# ╔═╡ 6528501f-3bb2-48a7-b24e-d32400033538
@bind pars_outbreak (
	PlutoUI.combine() do bind
		md"""
		R: $(bind(Slider(0:0.01:2.0,default=0.25;show_value=true))) \
		K: $(bind(Slider(0.01:0.02:10.0,default=8.71;show_value=true))) \
		P: $(bind(Slider(0.0:0.01:1.0,default=0.48;show_value=true))) \
		x0: $(bind(Slider(0:0.02:8.0,default=0.1;show_value=true)))
		"""
	end
)

# ╔═╡ a3270d8b-f167-4fb9-bcc5-c62c91b39649
flow1D(logoutbreak,pars_outbreak[4],300.0,pars_outbreak;xlims=[-0.2,8.0],title="Log with Outbreak")

# ╔═╡ 0d0d6d2f-56b2-4e01-b9a4-6e469eac90ad
html"""
<div style="position: relative; right: 0; top: 0; z-index: 300;"><iframe src="https://www.youtube.com/embed/MlwAI3BlDsU" width=500 height=250  frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe></div>
"""

# ╔═╡ 76cfd5c8-6298-4ae2-9d2d-9e60deb12f49
html"""
<div style="position: relative; right: 0; top: 0; z-index: 300;"><iframe src="https://www.youtube.com/embed/yUEXpUyi404" width=500 height=250  frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe></div>
"""

# ╔═╡ 0bf40552-b627-4fa7-b3ba-1d8c7ed4a568
html"""
<div style="position: relative; right: 0; top: 0; z-index: 250;"><iframe src="https://www.youtube.com/embed/Vu9oNWXv4Uk" width=500 height=250  frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe></div>
"""

# ╔═╡ 2ac364c2-cbdd-49b4-9f26-9fe89382be5e
TableOfContents(title="📚 Table of Contents", indent=true, depth=4, aside=true)

# ╔═╡ 6b6a5e19-5298-4969-b1c9-699aa1cb2996
html"""
<style>
input[type*="range"] {
	width: 50%;
}
</style>
"""

# ╔═╡ d38f2b14-9068-591d-b95d-0776b4cc9410
md"""
## Harvest threshold and recovery time

For $r,K>0$, constant harvesting has a saddle-node at
$H_c=rK/4$. Below it, the larger equilibrium attracts and the smaller one is
a threshold for extinction:

```math
N_\pm=\frac K2\left(1\pm\sqrt{1-\frac{4H}{rK}}\right),\qquad
\tau_{\rm recovery}=\frac{1}{r\sqrt{1-H/H_c}}.
```

The recovery time is obtained by linearizing at $N_+$, so it describes small
perturbations. Its divergence explains **critical slowing down** before collapse.
"""

# ╔═╡ a42dc65a-896a-566d-9717-c3cad70ad775
let r=pars_harvest[1], K=pars_harvest[2]
    if r==0
        md"With zero growth there is no positive sustainable constant harvest."
    else
        Hc=r*K/4
        hs=range(0,Hc;length=201)
        root=sqrt.(max.(0,1 .-hs./Hc))
        a=plot(hs,K/2 .* (1 .+root);label="stable",xlabel="H",ylabel="N*")
        plot!(a,hs,K/2 .* (1 .-root);label="unstable",linestyle=:dash)
        scatter!(a,[Hc],[K/2];label="saddle-node",color=:black)
        fractions=range(0,0.99;length=201)
        b=plot(fractions,1 ./ (r.*sqrt.(1 .-fractions));
            xlabel="H/Hc",ylabel="recovery time",legend=false)
        plot(a,b;layout=(1,2),size=(900,370),margin=5*Plots.mm)
    end
end

# ╔═╡ 87089c60-83a2-59b7-834f-a3e448867065
md"""
## Outbreaks, two potential wells, and hysteresis

Saturating predation can allow both a small and a large population to persist:

```math
\dot N=rN(1-N/K)-P\frac{N^2}{1+N^2},\qquad
V(N)=-\frac{rN^2}{2}+\frac{rN^3}{3K}+P(N-\arctan N).
```

Here $-\partial_N V=\dot N$. Two minima separated by a maximum mean two
attractors separated by an unstable population threshold. This differs from
constant harvesting: $N=0$ is an equilibrium and positive solutions stay nonnegative.
"""

# ╔═╡ 65428e0c-65b5-56d1-9ecf-4521f4292a12
let r=pars_outbreak[1], K=pars_outbreak[2], P=pars_outbreak[3]
    ns=range(0,max(K,1);length=501)
    f=r.*ns.*(1 .-ns./K) .-P.*ns.^2 ./ (1 .+ns.^2)
    V=-r.*ns.^2 ./2 .+r.*ns.^3 ./(3K) .+P.*(ns .-atan.(ns))
    a=plot(ns,f;xlabel="N",ylabel="dN/dt",legend=false)
    hline!(a,[0];color=:gray)
    b=plot(ns,V;xlabel="N",ylabel="V(N)",legend=false)
    plot(a,b;layout=(1,2),size=(900,370),margin=5*Plots.mm)
end

# ╔═╡ aae8b77f-8b71-5c3c-b0fb-bff2c33caacc
md"""
For positive equilibria, solving for the growth parameter gives a branch
without root-finding:

```math
r(N)=\frac{PN}{(1+N^2)(1-N/K)},\quad 0<N<K.
```

The next plot uses $K=8$, $P=1$ to show both folds. Stability is evaluated
from $f'(N)=r(1-2N/K)-2PN/(1+N^2)^2$. A slow sweep follows an attracting
branch until it ends, then jumps to the other one; the jump back occurs at a
different fold. Arrows would describe this sweep, not the direction of time along
the equilibrium curve itself.
"""

# ╔═╡ 2ff1e9d4-36e3-51cb-959c-c3b992448d3d
let K=8.0, P=1.0, ns=range(0.001,7.5;length=1200)
    rates=P.*ns ./((1 .+ns.^2).*(1 .-ns./K))
    slopes=rates.*(1 .-2 .*ns./K) .-2P.*ns./(1 .+ns.^2).^2
    stable=ifelse.(slopes.<0,ns,NaN)
    unstable=ifelse.(slopes.>=0,ns,NaN)
    f=plot(rates,stable;label="stable",xlabel="r",ylabel="N*",xlims=(0.3,0.8))
    plot!(f,rates,unstable;label="unstable",linestyle=:dash)
end

# ╔═╡ 46253e5c-cd04-54e8-bf52-4acea1053f62
md"""
**Try:** use $r=0.5$, $K=8$, $P=1$ and compare initial populations $0.1$ and $7$.
Locate their destinations in the potential. Then vary $r$ slowly in both
directions, using the previous final population as the next initial condition.
Does the same $r$ always give the same final population?
"""

# ╔═╡ Cell order:
# ╠═38c4e220-d708-11ec-3968-fbbd41c26155
# ╟─e2271efe-9a67-4496-891b-46f076b367f3
# ╟─62cba86d-4406-4d16-8715-b29bee7d3178
# ╟─6d64cc27-1059-4dd7-8f17-1c6fd9a2eca0
# ╠═4dc32b5a-5794-4d2d-844e-cdb80bec2e7b
# ╟─aa22c728-662c-4fdb-aae9-b3b411e311e2
# ╠═16a156b5-92ab-4e1c-ab70-c3b348186c57
# ╟─77766c99-4142-4acb-a855-ab133e67aa66
# ╟─ad67a3cb-a42f-451d-92a5-301559322c30
# ╟─ce4da9bf-316b-45d2-8504-f2a63cdda247
# ╠═29474ca0-839b-45de-a7bf-7d91c5ea59aa
# ╟─794936ae-d9ae-4319-994f-4282c99c4238
# ╟─8fc407f4-a242-4cdc-b7e0-062c8f5782d1
# ╟─b13d261e-1c1b-41b0-ae86-e20cd2304727
# ╟─f1f04dd4-0447-4e48-9bab-f34058a6a0f1
# ╟─c3757398-060c-412c-997b-58331b2dee04
# ╟─c803ed90-958c-4d0e-84d5-3e8cc9c8fbfe
# ╠═b4ea8364-a582-443c-bdb1-46e75a5c4d4e
# ╟─5018d342-0a51-4275-8c0f-1b5d253c2ec0
# ╟─5afa3aff-af15-429b-a0fb-9a0dfcad740d
# ╟─84a596c1-3a0a-4edd-a7ea-f3fcd8839c2a
# ╟─91b2718a-a921-4bae-99c0-91ec9c7e6479
# ╠═5a0b156e-b0b2-43b8-b15e-bcdb65c154bf
# ╟─a3270d8b-f167-4fb9-bcc5-c62c91b39649
# ╟─6528501f-3bb2-48a7-b24e-d32400033538
# ╟─0d0d6d2f-56b2-4e01-b9a4-6e469eac90ad
# ╟─76cfd5c8-6298-4ae2-9d2d-9e60deb12f49
# ╟─0bf40552-b627-4fa7-b3ba-1d8c7ed4a568
# ╟─2ac364c2-cbdd-49b4-9f26-9fe89382be5e
# ╟─6b6a5e19-5298-4969-b1c9-699aa1cb2996
# ╟─d38f2b14-9068-591d-b95d-0776b4cc9410
# ╠═a42dc65a-896a-566d-9717-c3cad70ad775
# ╟─87089c60-83a2-59b7-834f-a3e448867065
# ╠═65428e0c-65b5-56d1-9ecf-4521f4292a12
# ╟─aae8b77f-8b71-5c3c-b0fb-bff2c33caacc
# ╠═2ff1e9d4-36e3-51cb-959c-c3b992448d3d
# ╟─46253e5c-cd04-54e8-bf52-4acea1053f62
