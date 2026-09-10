### A Pluto.jl notebook ###
# v0.19.39

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

# ╔═╡ a860f150-f213-11ee-183b-c5168dc31a64
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, IntervalRootFinding, LinearAlgebra, Plots, PlutoUI
end

# ╔═╡ f572cb9a-a288-47c2-8967-325707299580
md"""
# LFF Excitable model with Noise

$dx = y\,dt + \eta\,dW_t$

$dy = [\epsilon_1-y+x(1+\epsilon_2x+y-x^2)]\,dt$

#![Bifuraction Diagram](https://i.imgur.com/jNWr2Od.png)
"""

# ╔═╡ 5c859d3f-b048-40e3-bfa3-a0eade5da658
function lff(du, u, p, t)
	du[1] = u[2]
	du[2] = p[1]-u[2]+u[1]*(1.0+p[2]*u[1]+u[2]-u[1]^2)
end

# ╔═╡ ccc11b0a-1686-4579-99e1-f50010d78a80
function noise(du,u,p,t)
	du[1] = p[3]
	du[2] = 0.0
end

# ╔═╡ 25208988-7b14-44d4-889a-4d166424a5ab
function root_from_e2(e1,e2)
	f(x) = e1 + x + e2*x^2 - x^3
	roots(f,-3..3)
end

# ╔═╡ e222dbce-506f-4e7f-b53a-704ec6463f3d
#saddle node
begin
	e2 = 0:0.01:2;
	x2 = (e2-sqrt.(e2.^2 .+3))/3
	sn = x2.*(x2.^2 .- e2 .* x2 .- 1.0)
end;	

# ╔═╡ f4d76d8e-f4de-430e-9424-32ec7a8239ab
md"""Compute the global bifurcation curve (slow): $(@bind track_global CheckBox(default=false))"""

# ╔═╡ 70cda775-ebe8-4b4d-af4b-1c5f357b0b4a
# trackeo manual de la bif global
#for e1=0:0.01:0.28
begin
	e2g = Float64[]
	e1g = track_global ? (0.1:0.0005:0.2845) : Float64[]
	for e1 = e1g
		mm = 1.8
		if e1>0.1
			e2v = e2g[end]-0.001
		else
			e2v = 0.34
		end	
		for attempt in 1:2000
			pp = root_from_e2(e1,e2v)
			x0 = maximum(mid(p.interval) for p in pp) # repulsor
			solt = solve(ODEProblem(lff, [x0+0.1, 0.0], (0.0, 300.0),[e1,e2v,0]),saveat=1.0)
			mm = maximum([x[1] for x in solt.u[end-100:end]])
			mm <= 0 && break
			e2v += 0.0001
		end	
		mm <= 0 || error("Global bifurcation search did not converge for ϵ1=$e1")
		push!(e2g,e2v)
	end	
end	

# ╔═╡ 0201ef70-5adc-4568-a364-12881d5b9163
md"""
ϵ1 : $(@bind ϵ1 Slider(0.2:0.002:0.3,default=0.26;show_value=true))\
ϵ2 : $(@bind ϵ2 Slider(0.2:0.002:1.0,default=0.44;show_value=true)) \
η : $(@bind η Slider(0.0:0.01:0.3,default=0.0;show_value=true)) \
"""

# ╔═╡ 8708d3ec-657f-4b91-905f-74ec563d2b58
begin
	plot(e2,sn)
	plot!(e2g,e1g)
	scatter!([ϵ2],[ϵ1])
end	

# ╔═╡ 0ef06404-e9fc-4688-8f24-861ea384a42d
begin
	sol = solve(SDEProblem(lff, noise, [1, 0.0], (0.0, 3000.0),[ϵ1,ϵ2,η]))
	plot(sol, idxs = (0,1),size=(1200,300),label="")
end	

# ╔═╡ 6ff35c81-91eb-466a-ae1a-e5c8f2172e26
function lff_cycle(du, u, p, t)
	(ϵ10,ϵ20,ϵ1A,ϵ2A,α,ω) = p 
	ϵ1 = ϵ10 + ϵ1A*sin(ω*t)*cos(α) + ϵ2A*cos(ω*t)*sin(α)
	ϵ2 = ϵ20 + ϵ2A*cos(ω*t)*cos(α) - ϵ1A*sin(ω*t)*sin(α)
	du[1] = u[2]
	du[2] = ϵ1-u[2]+u[1]*(1.0+ϵ2*u[1]+u[2]-u[1]^2)
end

# ╔═╡ 4eb71e59-395c-49d4-acac-c1a2682dd676
md"""
ϵ10 : $(@bind ϵ10 Slider(0.2:0.002:0.3,default=0.284;show_value=true))\
ϵ20 : $(@bind ϵ20 Slider(0.2:0.001:1.0,default=0.362;show_value=true)) \
ϵ1A : $(@bind ϵ1A Slider(0.0:0.01:0.2,default=0.0;show_value=true))\
ϵ2A : $(@bind ϵ2A Slider(0.0:0.01:1.0,default=0.0;show_value=true)) \
α : $(@bind α Slider(-0.2:0.001:0.2,default=0.0;show_value=true)) \
ω : $(@bind ω Slider(0.0:0.001:0.1,default=0.0;show_value=true)) \
"""

# ╔═╡ a441d4af-c8b8-4c75-95d5-14ba0152d223
begin
	ϕ = 0:0.01:2*pi
	#α = -pi/24
	plot(e2,sn)
	plot!(e2g,e1g)
	plot!(ϵ20 .+ ϵ2A*cos.(ϕ).*cos(α) .- ϵ1A*sin.(ϕ).*sin(α),ϵ10 .+ ϵ1A*sin.(ϕ).*cos(α) .+ ϵ2A*cos.(ϕ).*sin(α))
end	

# ╔═╡ 50b807b7-4b9b-4f4c-aa9b-034483157ec2
begin
	sol1 = solve(ODEProblem(lff_cycle, [-0.1, 0.0], (0.0, 3000.0),[ϵ10,ϵ20,ϵ1A,ϵ2A,α,ω]))
	plot(sol1, idxs = (0,1),size=(1200,300),label="")
end	

# ╔═╡ de0d0f5d-dd4f-45da-86b7-5027d82d5ef8
html"""
<style>
main {
    max-width: 1200px;
}
input[type*="range"] {
	width: 90%;
}
</style>
"""

# ╔═╡ 79c4f5ef-6b20-49bb-8e3e-df390a67a88d
md"""
# Double Hopf

"Easy Case"
"""

# ╔═╡ 3d8d06c3-87c5-4e75-91b7-8fcf94d554dc
function dhopf(du, u, p, t)
	(μ1,μ2,ω1,ω2,θ,δ) = p
	ξ1 = u[1]^2+u[2]^2
	ξ2 = u[3]^2+u[4]^2
	du[1] = μ1*u[1] - ω1*u[2] - u[1]*ξ1 - θ*u[1]*ξ2
	du[2] = ω1*u[1] + μ1*u[2] - u[2]*ξ1 - θ*u[2]*ξ2
	du[3] = μ2*u[3] - ω2*u[4] - u[3]*ξ2 - δ*u[3]*ξ1
	du[4] = ω2*u[3] + μ2*u[4] - u[4]*ξ2 - δ*u[4]*ξ1
end

# ╔═╡ 4ae3e029-dc5f-4c37-ad89-b8c5c40ba236
md"""
# Double Hopf with Noise
"""

# ╔═╡ 073062ad-2f9d-4b1c-bf8e-a50b4eb460f3
function noise2(du,u,p,t)
	du[1] = p[end-1]
	du[2] = 0.0
	du[3] = p[end]
	du[4] = 0.0
end

# ╔═╡ 71047ce7-71a6-4dea-b363-1111879ae2c0
md"""
μ1 : $(@bind μ1 Slider(-1.0:0.002:1.0,default=0.3;show_value=true))\
μ2 : $(@bind μ2 Slider(-1.0:0.001:1.0,default=0.4;show_value=true)) \
ω1 : $(@bind ω1 Slider(0.0:0.002:1.0,default=0.3;show_value=true))\
ω2 : $(@bind ω2 Slider(0.0:0.002:1.0,default=0.13;show_value=true))\
θ : $(@bind θ Slider(-1:0.1:3.0,default=1.1;show_value=true))\
δ : $(@bind δ Slider(-1:0.1:3.0,default=1.1;show_value=true))\
η1 : $(@bind η1 Slider(0.0:0.001:0.3,default=0.0;show_value=true))\
η2 : $(@bind η2 Slider(0.0:0.001:0.3,default=0.0;show_value=true))\
tmax = $(@bind tmax Slider(100.0:100.0:5000.0,default=500.0;show_value=true))\
"""

# ╔═╡ 848c35f4-0ba2-4ccf-b3fb-a56eaebcb25d
# ╠═╡ disabled = true
#=╠═╡
begin
	sol2 = solve(ODEProblem(dhopf, [0.035, -0.01,0.01,0.01], (0.0, 1000.0),[μ1,μ2,ω1,ω2,θ,δ]))
	plot(sol2, idxs = (0,1),size=(1200,300),label="1")
	plot!(sol2, idxs = (0,3),size=(1200,300),label="2")
end	
  ╠═╡ =#

# ╔═╡ 1558b198-afef-4f69-bb43-6fb367561f86
begin
	sol3 = solve(SDEProblem(dhopf, noise2, [0.035, -0.01,0.01,0.01], (0.0, tmax),[μ1,μ2,ω1,ω2,θ,δ,η1,η2]))
	plot(sol3, idxs = (0,2),size=(1200,300),label="1")
	plot!(sol3, idxs = (0,4),size=(1200,300),label="2")
end	

# ╔═╡ 66f9c93b-7d88-46c8-a365-1027878c50ff
md"""
## Double Hopf standard form
"""

# ╔═╡ 9bffe395-0595-4b46-981f-6c313e1e5512
function dhopfs(du, u, p, t)
	(μ1,μ2,k1,k2,θ,δ) = p
	ξ1 = u[1]^2+u[2]^2
	ξ2 = u[3]^2+u[4]^2
	du[1] = u[2]
	du[2] = -k1*u[1] + μ1*u[2] - u[2]*ξ1 - θ*u[2]*ξ2
	du[3] = u[4]
	du[4] = -k2*u[3] + μ2*u[4] - u[4]*ξ2 - δ*u[4]*ξ1
end

# ╔═╡ 4d63ba5d-9d53-4a07-bfee-658c6cd9e9eb
md"""
Standard Form \
μ1 : $(@bind μ1s Slider(-1.0:0.002:1.0,default=0.3;show_value=true))\
μ2 : $(@bind μ2s Slider(-1.0:0.001:1.0,default=0.4;show_value=true)) \
k1 : $(@bind k1 Slider(0.0:0.002:1.0,default=0.3;show_value=true))\
k2 : $(@bind k2 Slider(0.0:0.002:1.0,default=0.13;show_value=true))\
θ : $(@bind θs Slider(-1:0.1:3.0,default=1.1;show_value=true))\
δ : $(@bind δs Slider(-1:0.1:3.0,default=1.1;show_value=true))\
tmax = $(@bind tmaxs Slider(100.0:100.0:5000.0,default=500.0;show_value=true))\
"""

# ╔═╡ 1f40e730-0957-47cb-9bb7-175c09db8fc7
begin
	sol4 = solve(ODEProblem(dhopfs, [0.035, -0.01,0.01,0.01], (0.0, tmaxs),[μ1s,μ2s,k1,k2,θs,δs]))
	plot(sol4, idxs = (0,1),size=(1200,300),label="1")
	plot!(sol4, idxs = (0,3),size=(1200,300),label="2")
end	

# ╔═╡ Cell order:
# ╠═a860f150-f213-11ee-183b-c5168dc31a64
# ╠═f572cb9a-a288-47c2-8967-325707299580
# ╠═5c859d3f-b048-40e3-bfa3-a0eade5da658
# ╠═ccc11b0a-1686-4579-99e1-f50010d78a80
# ╠═25208988-7b14-44d4-889a-4d166424a5ab
# ╠═e222dbce-506f-4e7f-b53a-704ec6463f3d
# ╟─f4d76d8e-f4de-430e-9424-32ec7a8239ab
# ╠═70cda775-ebe8-4b4d-af4b-1c5f357b0b4a
# ╠═8708d3ec-657f-4b91-905f-74ec563d2b58
# ╟─0201ef70-5adc-4568-a364-12881d5b9163
# ╠═0ef06404-e9fc-4688-8f24-861ea384a42d
# ╠═6ff35c81-91eb-466a-ae1a-e5c8f2172e26
# ╠═a441d4af-c8b8-4c75-95d5-14ba0152d223
# ╟─4eb71e59-395c-49d4-acac-c1a2682dd676
# ╠═50b807b7-4b9b-4f4c-aa9b-034483157ec2
# ╟─de0d0f5d-dd4f-45da-86b7-5027d82d5ef8
# ╟─79c4f5ef-6b20-49bb-8e3e-df390a67a88d
# ╠═3d8d06c3-87c5-4e75-91b7-8fcf94d554dc
# ╠═848c35f4-0ba2-4ccf-b3fb-a56eaebcb25d
# ╟─4ae3e029-dc5f-4c37-ad89-b8c5c40ba236
# ╠═073062ad-2f9d-4b1c-bf8e-a50b4eb460f3
# ╟─71047ce7-71a6-4dea-b363-1111879ae2c0
# ╠═1558b198-afef-4f69-bb43-6fb367561f86
# ╟─66f9c93b-7d88-46c8-a365-1027878c50ff
# ╠═9bffe395-0595-4b46-981f-6c313e1e5512
# ╟─4d63ba5d-9d53-4a07-bfee-658c6cd9e9eb
# ╠═1f40e730-0957-47cb-9bb7-175c09db8fc7
