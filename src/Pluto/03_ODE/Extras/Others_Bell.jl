### A Pluto.jl notebook ###
# v0.19.26

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

# ╔═╡ 8b3b5b86-3474-11ee-08d0-452e0ff17997
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", "..", ".."))
    Pkg.instantiate()
    using Plots, DifferentialEquations, PlutoUI
end

# ╔═╡ f8ecadc1-ea69-4945-a89e-07a5fb5f9380
md"""
# Bell and Clapper dynamics

This is based on "The Dynamics of a Ringing Church Bell"
J. Woodhouse et al.

In Equation (11) they display the equations of motion for the angle of the bell relative to the vertical $\theta$ and the angle of the clapper relative to the bell axis $\phi$ in normalized units:

$\dot{\theta} = v_{\theta}$

$\dot{v_{\theta}} = - sin(\theta)$

$\dot{\phi} = v_{\phi}$

$\dot{v_{\phi}} = sin(\theta)\left(1+\frac{r}{L_c}cos(\phi)\right) -\frac{r}{L_c}sin(\phi) v_{\theta}^2 - \frac{L_b}{L_c}sin(\theta+\phi)-F_c$

$F_c = \frac{k}{2}\operatorname{sign}(\phi)(|\phi|-\phi_{max})^2+\gamma v_\phi \quad : \quad|\phi|>\phi_{max}$

$F_c=0$ otherwise. This compliant contact and its damping are additions to the free-motion equations.

where $L_b=I_b/(Ma)$ and $L_c=I_c/(mb)$ and the nondimensional time is $\tau=t\sqrt{\frac{g}{L_b}}$.

The bell (clapper) has mass $M$ ($m$) moment of inerta $I_b$ ($I_c$) and the center of mass is at distance $a$ ($b$) from the pivot.
"""

# ╔═╡ 593d3f09-9a1f-4014-9282-3c86ec160428
function bell_clapper!(du,u,p,t)
	(R,L,k,ϕₘ,γ) = p
	if u[3] > ϕₘ
		Fₖ = 0.5*k*(u[3]-ϕₘ)^2
		γₖ = γ 
	elseif u[3] < -ϕₘ
		Fₖ = -0.5*k*(u[3]+ϕₘ)^2
		γₖ = γ 
	else
		Fₖ = 0 
		γₖ = 0
	end	
	du[1] = u[2]
	du[2] = -sin(u[1])
	du[3] = u[4]
	du[4] = sin(u[1])*(1+R*cos(u[3]))-R*sin(u[3])*u[2]^2-L*sin(u[1]+u[3])-Fₖ-γₖ*u[4]
end

# ╔═╡ 8974c1c7-97e6-41db-8f1d-950ed6fbdd5d
md"""In the hard-contact example, the clapper rebounds with restitution coefficient 0.7. As impacts accumulate, it rests against the bell until the free acceleration points back into the allowed region."""

# ╔═╡ b198d4d5-9e2f-4dda-93a7-145c8b37dc83


# ╔═╡ 1191acf0-4e5d-5065-a6e4-a1fa7637f3fe
md"""
Hard contact uses the free-flight equations below, with ``R=r/L_c`` and ``L=L_b/L_c``.

```math
\begin{aligned}
\dot{\theta} &= v_\theta, & \dot{v}_\theta &= -\sin\theta,\\
\dot{\phi} &= v_\phi, & \dot{v}_\phi &= \sin\theta(1+R\cos\phi)-R\sin\phi\,v_\theta^2-L\sin(\theta+\phi),\\
v_\phi^+ &= -0.7v_\phi^- \quad\text{at an outward impact on } |\phi|=26\pi/180.
\end{aligned}
```

At resting contact, the implementation sets both clapper derivatives to zero while the free acceleration points into the wall. Impacts with speed below 10⁻⁷ are treated as rest; motion resumes when acceleration points inward.
"""

# ╔═╡ 6499063c-819f-4fe8-97a9-aeb2e2a3d36d
function bell_clapper_hard!(du,u,p,t)
	(R,L) = p
	du[1] = u[2]
	du[2] = -sin(u[1])
	du[3] = u[4]
	du[4] = sin(u[1])*(1+R*cos(u[3]))-R*sin(u[3])*u[2]^2-L*sin(u[1]+u[3])
	if abs(u[3]) >= 26*pi/180 - 1e-12 && abs(u[4]) < 1e-7 && u[3]*du[4] > 0
		du[3] = 0
		du[4] = 0
	end
end

# ╔═╡ ec45044a-7bb4-4289-b22e-f5ad85cd9f8f
function collision(u,t,integrator)
  abs(u[3])-26*pi/180
end

# ╔═╡ 853145f7-29cb-4f6e-afe6-347f187fd93b
function bounce!(integrator)
	integrator.u[3] = sign(integrator.u[3])*26*pi/180
	integrator.u[4] = abs(integrator.u[4]) < 1e-7 ? 0.0 : -0.7*integrator.u[4]
	set_proposed_dt!(integrator, min(0.01, abs(integrator.dt)))
end	

# ╔═╡ 1ab24db2-221b-44a9-86a4-386ef00a2302


# ╔═╡ ea5401b7-09dd-48ad-890a-773d9a22b75e
ϕₘ = 26*pi/180

# ╔═╡ ae98a306-9c10-4fee-84f4-db3e94b9e753
p=[0.44,1.685,1000000,ϕₘ,10]

# ╔═╡ 5506e9e7-1374-4db3-84db-6bdf3608aacc
probh = ODEProblem(bell_clapper_hard!,[-3.0,0.0,0.41,0.0],(0,50.0),p[1:2]);

# ╔═╡ 924ba207-bcba-44d7-a10e-64cf626df977
solh = solve(probh, Rosenbrock23(); callback=ContinuousCallback(collision,bounce!,nothing), abstol=1e-9, reltol=1e-8, dtmax=0.05);

# ╔═╡ 6a18af4a-7033-4107-9bc7-fa0a9224efcf
md"""
θ₀ $(@bind θ₀ Slider(-3.14:0.001:3.14,default=3.0;show_value=true)) 
ϕ₀ $(@bind ϕ₀ Slider(-0.453:0.001:0.453,default=0.4;show_value=true)) \
r $(@bind r Slider(0.001:0.001:2.0,default=0.445;show_value=true)) 
Lb $(@bind Lb Slider(0.001:0.001:3.0,default=1.68;show_value=true)) \
k $(@bind k Slider(1000:1000:100000,default=10000;show_value=true)) 
μ $(@bind μ Slider(0:0.001:20.0,default=0.0;show_value=true)) \
tmax $(@bind tmax Slider(5.0:1.0:100.0,default=5.0;show_value=true))
"""

# ╔═╡ 31ca30e2-7c80-442a-afe3-7c0ad0b859ec
begin
	prob = ODEProblem(bell_clapper!,[θ₀,0.0,ϕ₀,0],(0,tmax),[r,Lb,k,ϕₘ,μ]);
	sol = solve(prob,Rosenbrock23());
	plot(sol,idxs=(0,1))
	plot!(sol,idxs=(0,3))
	plot!(sol,idxs=(0,4))
end	

# ╔═╡ Cell order:
# ╠═8b3b5b86-3474-11ee-08d0-452e0ff17997
# ╟─f8ecadc1-ea69-4945-a89e-07a5fb5f9380
# ╠═593d3f09-9a1f-4014-9282-3c86ec160428
# ╟─8974c1c7-97e6-41db-8f1d-950ed6fbdd5d
# ╠═b198d4d5-9e2f-4dda-93a7-145c8b37dc83
# ╟─1191acf0-4e5d-5065-a6e4-a1fa7637f3fe
# ╠═6499063c-819f-4fe8-97a9-aeb2e2a3d36d
# ╠═ec45044a-7bb4-4289-b22e-f5ad85cd9f8f
# ╠═853145f7-29cb-4f6e-afe6-347f187fd93b
# ╠═1ab24db2-221b-44a9-86a4-386ef00a2302
# ╠═ea5401b7-09dd-48ad-890a-773d9a22b75e
# ╠═ae98a306-9c10-4fee-84f4-db3e94b9e753
# ╠═5506e9e7-1374-4db3-84db-6bdf3608aacc
# ╠═924ba207-bcba-44d7-a10e-64cf626df977
# ╟─6a18af4a-7033-4107-9bc7-fa0a9224efcf
# ╠═31ca30e2-7c80-442a-afe3-7c0ad0b859ec
