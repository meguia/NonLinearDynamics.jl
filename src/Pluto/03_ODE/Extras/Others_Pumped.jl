### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 8ef08220-20ff-11ee-2e03-df7546b1b262
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots
end

# ╔═╡ a4e679e0-7436-4ddc-a104-90afd6a2b71e
step(x,ϵ) = (x/sqrt(ϵ+x*x)+1)/2

# ╔═╡ e64a1f29-e152-41e9-ace1-6d6466caa64e
function facc(du,u,p,t)
	(m,ϵ,K,c) = p
	x = sin(c*t)
	du[1] = u[2]
	du[2] = step(x-u[1],ϵ)*(x-u[1])*K/m
	du[3] = c*cos(c*t)
end

# ╔═╡ 7e1110f2-8b7e-4b58-a46f-5acd7a7adcd6
begin
	u0 = [0.0,0.0,0.0]
	tspan = (0.0, 20.0)
	p = [50.0,1e-8,1.0,1.0]
	prob = ODEProblem(facc, u0, tspan,p)
	sol = solve(prob, Tsit5(), reltol = 1e-8, abstol = 1e-8)
end;

# ╔═╡ 98991644-16f7-4f41-a437-2a767c3ef797
force(t,x0,x) = (t,step(x-x0,p[2])*(x-x0)*p[3]/p[1])

# ╔═╡ 1a154b1c-0b2c-4c66-bbe8-1fe0d60b8d0c
ω = sqrt(p[3]/p[1])

# ╔═╡ 1c325481-010b-4883-bbea-bb16c19b0e92
c = p[4]

# ╔═╡ 53e10af6-e61c-4323-ab1f-e9769d6c2528
begin
	t = 0:0.01:20
	x0 = @. ω/(ω^2-c^2)*(-c*sin(ω*t)+ω*sin(c*t))
	v0 = @. ω^2*c/(ω^2-c^2)*(cos(c*t)-cos(ω*t))
	a0 = @. ω^2*c/(ω^2-c^2)*(ω*sin(ω*t)-c*sin(c*t))
end		

# ╔═╡ 5a2fc8c9-077e-4bf4-9201-9be2d5bb4da0
begin
	plot(sol,idxs=(0,1),label="x₀")
	plot!(t,sin.(c*t),label="x")
end	

# ╔═╡ 60040ce3-cc66-4844-8185-0684d235b0ec
begin
	plot(sol,idxs=(0,2),label="v₀")
	#plot!(sol,idxs=(0,3),label="x")
	plot!(sol,idxs=(force,0,1,3),label="force")
	plot!(t,x0,label="X0")
	plot!(t,a0,label="a0")
	plot!(t,v0,label="V0")
end

# ╔═╡ Cell order:
# ╠═8ef08220-20ff-11ee-2e03-df7546b1b262
# ╠═a4e679e0-7436-4ddc-a104-90afd6a2b71e
# ╠═e64a1f29-e152-41e9-ace1-6d6466caa64e
# ╠═7e1110f2-8b7e-4b58-a46f-5acd7a7adcd6
# ╠═98991644-16f7-4f41-a437-2a767c3ef797
# ╠═1a154b1c-0b2c-4c66-bbe8-1fe0d60b8d0c
# ╠═1c325481-010b-4883-bbea-bb16c19b0e92
# ╠═53e10af6-e61c-4323-ab1f-e9769d6c2528
# ╠═5a2fc8c9-077e-4bf4-9201-9be2d5bb4da0
# ╠═60040ce3-cc66-4844-8185-0684d235b0ec
