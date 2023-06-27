### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 8984e3bc-39c9-47bf-886e-48e9bc9626a8
import Pkg; Pkg.add("ForwardDiff")

# ╔═╡ 96e6e962-71f3-4c7a-8c2d-ab7a6cc7baf4
using PlutoUI, DifferentialEquations

# ╔═╡ 3c1b2921-4fd2-4481-937b-a47b597b8f46
using Plots; plotly()

# ╔═╡ bc3b5438-df2a-481d-9cfb-c86eae6ebe9b
include("../NLD_utils.jl")

# ╔═╡ b15c712e-61bd-4568-8ba8-579212bb6e6c
html"""
<style>
	main {
		margin: 0 auto;
		max-width: 2000px;
    	padding-left: max(160px, 10%);
    	padding-right: max(160px, 10%);
	}
</style>
"""

# ╔═╡ 57f14c67-8e05-4c8a-90e5-6a4f3377506e
TableOfContents()

# ╔═╡ 92bbe726-f115-40fe-8a15-0c5bea8589df
md"""

# Dynamical Systems on the Circle S^1

## Adler's Equation
"""

# ╔═╡ b2ea050c-7f55-48ba-a165-41709ae5f78e
adler(u,p,t) = p[1]-p[2]*sin(u)

# ╔═╡ 99f2f5e8-f91d-4de9-b0b9-3e0d1c7b2ec9
begin
	ω = 1.1
	sol = solve(ODEProblem(adler, 0.1, (0,50.0), [ω,1]))
	plot(sol,idxs=((t,x)->(t,mod(x,2*pi)/(2*pi)),0,1),label="phase")
	plot!(sol,idxs=((t,x)->(t,sin(x)),0,1),label="sin",size=(1400,300))
end	

# ╔═╡ 540a0af5-321f-41a4-8a90-756b798c84f2
plot(sol,idxs=(0,1),label="unwrapped phase",size=(1400,300))

# ╔═╡ 12f373d5-bc45-480b-93ba-d83b0c09f873
function condition(u, t, integrator) 
    (u-pi)*(u+pi)
end

# ╔═╡ 416e0aed-a7bc-42e0-8a32-512471bc56f3
function affect!(integrator)
    integrator.u -= sign(integrator.u)*2*pi
end

# ╔═╡ 5ec39c31-afe9-41a0-8b4b-f6e68928144c
md"""
# Osciladores de Adler acoplados
"""

# ╔═╡ 3169f185-9d52-4a80-95a0-dde68804ac84
function adler_coupled!(du,u,p,t)
    (ω1,ω2,k) = p
    du[1] = ω1-sin(u[1])+k*sin(u[2]-u[1])
    du[2] = ω2-sin(u[2])+k*sin(u[1]-u[2])
    return
end

# ╔═╡ 9e9dbfe7-668e-494a-8159-710ca914f2de
begin
	(ω1,ω2,k) = [1.2,1.3,0.1]
	tmax = 200.0
	t = 0:0.1:tmax
	sol2=solve(ODEProblem(adler_coupled!, [0.0,0.0], (0,tmax), [ω1,ω2,k]),Rosenbrock23())
    p1 = plot(sol2,idxs=((t,x)->(t,mod(x,2*pi)),0,1),size=(1400,300))
    plot!(p1,sol2,idxs=((t,x)->(t,mod(x,2*pi)),0,2))
end	

# ╔═╡ 2912533c-f785-4ca2-be22-5efeb5947567
function affect!(integrator, idx)
	integrator.u[idx] -= sign(integrator.u[idx])*2*pi
end

# ╔═╡ f1d13184-8e0c-41b0-af66-7e1925bf6b91
function condition(out,u, t, integrator) 
    out[1] = (u[1]-pi)*(u[1]+pi)
	out[2] = (u[2]-pi)*(u[2]+pi)
end

# ╔═╡ 55ed709b-bc02-4752-b026-8e38e6f5d17f
begin
	sol1 = solve(ODEProblem(adler, 0.1, (0,50.0), [ω,1]),callback = ContinuousCallback(condition, affect!))
	plot(sol1,idxs=(0,1),label="unwrapped phase",size=(1400,300))
end	


# ╔═╡ da35c5c5-90ee-4335-82b7-49ae7bde66c7
begin
	cb = VectorContinuousCallback(condition, affect!, 2)
	sol3 = solve(ODEProblem(adler_coupled!, [0.0,0.0], (0,tmax), [ω1,ω2,k]),Rosenbrock23(),callback = cb)
	x1 = [sol3(ti)[1] for ti in t]
	y1 = [sol3(ti)[2] for ti in t];
	sol3b = solve(ODEProblem(adler_coupled!, [0.0,2.0], (0,tmax), [ω1,ω2,k]),Rosenbrock23(),callback = cb)
	x2 = [sol3b(ti)[1] for ti in t]
	y2 = [sol3b(ti)[2] for ti in t];
    #plot(x1,x2,size=(500,500))
end	

# ╔═╡ f520b8f5-cdfa-44ce-8dca-dcd8bf982111
function traj_torus(θ1,θ2,R)
	X = (R .+ cos.(θ2)).*cos.(θ1)
	Y = (R .+ cos.(θ2)).*sin.(θ1)
	Z = sin.(θ2)
	return (X,Y,Z)
end

# ╔═╡ 66fef05b-4acb-4f13-83a1-2c46e9ead2f8
begin
	Θ₁ = -π:0.01:π
	Θ₂ = -π:0.01:π
	X_torus = [(2.5 + cos(Θ₂ᵢ)) * cos(Θ₁ᵢ) for Θ₁ᵢ in Θ₁, Θ₂ᵢ in Θ₂]
	Y_torus = [(2.5 + cos(Θ₂ᵢ)) * sin(Θ₁ᵢ) for Θ₁ᵢ in Θ₁, Θ₂ᵢ in Θ₂]
	Z_torus = [sin(Θ₂ᵢ) for Θ₁ᵢ in Θ₁, Θ₂ᵢ in Θ₂]
	surface(X_torus, Y_torus, Z_torus, alpha=0.8, colorbar=:none, legend = false);
	X1,Y1,Z1 = traj_torus(x1,y1,2.5)
	X2,Y2,Z2 = traj_torus(x2,y2,2.5)
	plot!(X1,Y1,Z1,linecolor=:blue,linewidth=2,size=(1200,600))
	plot!(X2,Y2,Z2,linecolor=:red,linewidth=2,size=(1200,600))
end	

# ╔═╡ f67da2ec-7816-420a-96fe-aa6d72ba0aae
md"""

# Ecuaciones lineales en el toro

"""

# ╔═╡ a4fb3bc3-db31-4516-91b7-2c04abb169f5
function lin_torus!(du,u,p,t)
    (a,b,c,d) = p
    du[1] = a*sin(u[1])+b*sin(u[2])
    du[2] = c*sin(u[1])+d*sin(u[2])
    return
end

# ╔═╡ 16968a95-d19f-4aa1-b3f4-7932971d9796
begin
	(a,b,c,d) = [-0.0,-1,1,-0.0]
	t2 = 100.0
	sol4 = solve(ODEProblem(lin_torus!, [3.1,0.1], (0,t2), [a,b,c,d]),Tsit5(),callback = cb)
	θ_1 = [sol4(ti)[1] for ti in 0:0.01:t2]
	θ_2 = [sol4(ti)[2] for ti in 0:0.01:t2]
    plot(θ_1,θ_2)
end	

# ╔═╡ a867d3fc-3373-4937-af30-fea2aab10f79
begin
	surface(X_torus, Y_torus, Z_torus, alpha=0.8, colorbar=:none, legend = false);
	X3,Y3,Z3 = traj_torus(θ_1,θ_2,2.5)
	plot!(X3,Y3,Z3,linecolor=:black,linewidth=2,size=(1200,600))
end	

# ╔═╡ Cell order:
# ╠═8984e3bc-39c9-47bf-886e-48e9bc9626a8
# ╠═96e6e962-71f3-4c7a-8c2d-ab7a6cc7baf4
# ╠═3c1b2921-4fd2-4481-937b-a47b597b8f46
# ╠═b15c712e-61bd-4568-8ba8-579212bb6e6c
# ╠═bc3b5438-df2a-481d-9cfb-c86eae6ebe9b
# ╟─57f14c67-8e05-4c8a-90e5-6a4f3377506e
# ╟─92bbe726-f115-40fe-8a15-0c5bea8589df
# ╠═b2ea050c-7f55-48ba-a165-41709ae5f78e
# ╠═99f2f5e8-f91d-4de9-b0b9-3e0d1c7b2ec9
# ╠═540a0af5-321f-41a4-8a90-756b798c84f2
# ╠═12f373d5-bc45-480b-93ba-d83b0c09f873
# ╠═416e0aed-a7bc-42e0-8a32-512471bc56f3
# ╠═55ed709b-bc02-4752-b026-8e38e6f5d17f
# ╟─5ec39c31-afe9-41a0-8b4b-f6e68928144c
# ╠═3169f185-9d52-4a80-95a0-dde68804ac84
# ╠═9e9dbfe7-668e-494a-8159-710ca914f2de
# ╟─2912533c-f785-4ca2-be22-5efeb5947567
# ╟─f1d13184-8e0c-41b0-af66-7e1925bf6b91
# ╠═da35c5c5-90ee-4335-82b7-49ae7bde66c7
# ╟─f520b8f5-cdfa-44ce-8dca-dcd8bf982111
# ╠═66fef05b-4acb-4f13-83a1-2c46e9ead2f8
# ╟─f67da2ec-7816-420a-96fe-aa6d72ba0aae
# ╠═a4fb3bc3-db31-4516-91b7-2c04abb169f5
# ╠═16968a95-d19f-4aa1-b3f4-7932971d9796
# ╠═a867d3fc-3373-4937-af30-fea2aab10f79
