### A Pluto.jl notebook ###
# v0.19.5

using Markdown
using InteractiveUtils

# ╔═╡ 94d7d7a2-d9ee-11ec-344c-0f555ffdf705
using Pkg; Pkg.activate("../..")

# ╔═╡ 6a3d6236-a28d-4040-aa15-d5ef8ddbe3ed
using PlutoUI, Plots, DifferentialEquations, NonLinearDynamics

# ╔═╡ 5b2eef3e-cd79-4989-8c94-e2dcc215bc2a
md"""
# Lorenz System

This is undoubtedly the most famous dynamical system featuring chaos. The origin goes back to an original work done by Edward Lorenz in 1963, proposing a reduction to three differential equations of an already simplified model of atmospheric convection (the original model had 12 equations).

In the model, the lower layer of the atmosphere is at a higher temperature than the upper layer. In the case where the temperature difference is small there is a temperature gradient (a linear variation of temperature with height) but above a certain critical value the hot air rises and the cold air falls forming convection rolls as seen in the figure.
"""

# ╔═╡ d7c7c3dd-8eab-4a9d-b372-804a8b2c9477
html"""
<div>
<img src="https://i.imgur.com/ukyDvpM.png" width="700px">
</div>
"""

# ╔═╡ c74fb54e-b49b-4687-8c0d-5ada8769ddbf
md"""
For a history of how Lorenz discovered the sensitivity to the initial conditions of this system from a truncation of the numerical simulation printed on paper and the true origin of the term "butterfly effect" see

$\dot{x} = \sigma(y-x)$

$\dot{y} = \rho x - y - xz$

$\dot{z} = xy - \beta z$

In this system, the variable $x$ corresponds to the intensity of convection (how fast the coils rotate), the variable $y$ to the temperature difference between the rising hot air stream and the falling cold one, and $z$ to the deviation of the linear temperature variation with height. The parameters also have physical significance, although there are not ver common outside the domain of fluid dynamics: $\sigma$ is the Prandtl number, $\rho$ the Rayleigh number and $\beta$ a geometric factor.

In order to obtain chaos, the values $\sigma=10$, $\rho=28$ and $\beta=8/3$ are traditionally used, although these values do not correspond to any particular physical system. For high Rayleigh number values $\rho>1$ the truncation made by Lorenz to only three modes is no longer valid. Therefore, although the model can explain the origin of the convection rolls (for $\rho=1$) the chaotic regime does not correspond to the physical model. However, being a relatively simple system to analyze and mainly for historical reasons, the Lorenz model became the "model of models" displaying chaotic behavior.

"""

# ╔═╡ 80630ff0-64d2-49c0-b879-cf5de0d8fad0
function lorenz!(du,u,p,t)
    (σ,ρ,β)=p
    du[1]=σ*(u[2]-u[1])
    du[2]=ρ*u[1]-u[2]-u[1]*u[3]
    du[3]=u[1]*u[2]-β*u[3]
    du
end    

# ╔═╡ c2338e97-1c4c-4b61-bb83-794a54db455f
sol = solve(ODEProblem(lorenz!,[0.1,0.1,5.0],(0.0,100),[10.0,24.5,8/3]));

# ╔═╡ 4de8b3eb-4a77-47d1-b7cc-7a2d5c5bfd8d
plotly()

# ╔═╡ 40597d23-a13c-4954-8fca-d6423e6b36a0
plot(sol,vars=(1,2,3),label="lorenz")

# ╔═╡ 7ab51dfc-d044-47f3-a725-695ece51653b
gr()

# ╔═╡ 3592613a-538c-49cd-a27e-65d41b8b9f5d
begin 
	p1 = plot(sol,vars=(1,2,3),legend=false,title="Lorenz")
	p2 = plot(sol,vars=(0,1),label="x")
	p3 = plot(sol,vars=(0,2),label="y")
	p4 = plot(sol,vars=(0,3),label="z")
	plot(p1,p2,p3,p4,layout=@layout [a{0.5w} grid(3,1)])
end	

# ╔═╡ e82494b8-3ace-4bb9-851a-0be8336233ad


# ╔═╡ Cell order:
# ╠═94d7d7a2-d9ee-11ec-344c-0f555ffdf705
# ╠═6a3d6236-a28d-4040-aa15-d5ef8ddbe3ed
# ╟─5b2eef3e-cd79-4989-8c94-e2dcc215bc2a
# ╟─d7c7c3dd-8eab-4a9d-b372-804a8b2c9477
# ╟─c74fb54e-b49b-4687-8c0d-5ada8769ddbf
# ╠═80630ff0-64d2-49c0-b879-cf5de0d8fad0
# ╠═c2338e97-1c4c-4b61-bb83-794a54db455f
# ╠═4de8b3eb-4a77-47d1-b7cc-7a2d5c5bfd8d
# ╠═40597d23-a13c-4954-8fca-d6423e6b36a0
# ╠═7ab51dfc-d044-47f3-a725-695ece51653b
# ╠═3592613a-538c-49cd-a27e-65d41b8b9f5d
# ╠═e82494b8-3ace-4bb9-851a-0be8336233ad
