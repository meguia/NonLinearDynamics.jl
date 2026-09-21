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

# ╔═╡ 96e6e962-71f3-4c7a-8c2d-ab7a6cc7baf4
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using PlutoUI, DifferentialEquations
    using Plots; plotly()
end

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

# ╔═╡ e298f3c0-ecc4-53e2-a709-310b68dca883
md"""
Adler phase equation; the examples below use a = 1.

```math
\begin{aligned}
\dot{\theta}=\omega-a\sin\theta.
\end{aligned}
```
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

# ╔═╡ 796874fc-ebdd-508e-859a-8810320684e7
md"""
## From a fixed phase to a rotating phase

Adler's equation is a flow on a circle:

```math
\dot\theta=\omega-a\sin\theta,\qquad
T=\int_0^{2\pi}\frac{d\theta}{\omega-a\sin\theta}
  =\frac{2\pi}{\sqrt{\omega^2-a^2}}\quad(\omega>a\geq0).
```

For $0<\omega<a$ there are two fixed phases, one stable and one unstable.
At $\omega=a$ they collide on the circle; just above it the trajectory spends
most of each turn near the vanished pair. This is a **saddle-node on an
invariant circle** (SNIC, also called SNIPER). Its period diverges, unlike the
finite onset period of a generic Hopf bifurcation.

We integrate the unwrapped phase and wrap it modulo $2\pi$ only for display.
"""

# ╔═╡ 9f854251-4be9-5469-ad55-34d89e30906c
function course_adler!(du,u,p,t)
    ω,a=p
    du[1]=ω-a*sin(u[1])
    nothing
end

# ╔═╡ 6a34a69d-9b35-5d79-b0a3-0b85a3c12e56
md"""
ω (a = 1): $(@bind course_omega Slider(0.5:0.01:1.5; default=1.05, show_value=true))
"""

# ╔═╡ 147c066b-0ca4-5ef5-b793-08fc6561943e
course_adler_u0 = [0.0]

# ╔═╡ ca40de0e-805a-5b47-bd5c-edb18d4b7e05
course_adler_tspan = (0.0,80.0)

# ╔═╡ 3a787114-4f7e-5cf8-b455-ad75e2e35c4a
course_adler_p = [course_omega,1.0]

# ╔═╡ 853c2184-d82c-5058-a2d2-15ed202205c2
course_adler_prob = ODEProblem(course_adler!, course_adler_u0, course_adler_tspan, course_adler_p)

# ╔═╡ 989c3787-fda2-5d26-a22d-ad1148fa5b86
course_adler_sol = solve(course_adler_prob, Tsit5(); abstol=1e-9, reltol=1e-8);

# ╔═╡ 6cdfa56c-bd82-5726-98ad-2bac14a20221
let ts=range(course_adler_tspan...;length=1601), θ=range(0,2pi;length=401)
    phases=course_adler_sol(ts;idxs=1).u
    a=plot(θ,course_omega .-sin.(θ);xlabel="θ",ylabel="dθ/dt",legend=false)
    hline!(a,[0];color=:gray)
    b=plot(ts,mod.(phases,2pi);xlabel="t",ylabel="θ mod 2π",legend=false)
    c=plot(ts,sin.(phases);xlabel="t",ylabel="sin θ",legend=false)
    plot(a,b,c;layout=(1,3),size=(1000,340),margin=5*Plots.mm)
end

# ╔═╡ 9cb434bb-7639-5753-a544-a0d2f2d500f3
let frequencies=range(1.001,1.5;length=400)
    f=plot(frequencies,2pi ./sqrt.(frequencies.^2 .-1);
        xlabel="ω (a = 1)",ylabel="period T",label="exact rotating solution")
    if course_omega>1
        scatter!(f,[course_omega],[2pi/sqrt(course_omega^2-1)];label="selected ω")
    end
    f
end

# ╔═╡ 5ec39c31-afe9-41a0-8b4b-f6e68928144c
md"""
# Osciladores de Adler acoplados
"""

# ╔═╡ 85398bab-e2c9-565a-848f-ff022caec73d
md"""
Two coupled phases, each defined modulo 2π.

```math
\begin{aligned}
\dot{\theta}_1 &= \omega_1-\sin\theta_1+k\sin(\theta_2-\theta_1),\\
\dot{\theta}_2 &= \omega_2-\sin\theta_2+k\sin(\theta_1-\theta_2).
\end{aligned}
```
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
	plot(sol1,idxs=(0,1),label="wrapped phase",size=(1400,300))
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

# ╔═╡ 5451b0b3-d3f1-51fd-88ee-cf2e7047a0c9
md"""
A sine vector field on the torus.

```math
\begin{aligned}
\dot{\theta}_1 &= a\sin\theta_1+b\sin\theta_2,\\
\dot{\theta}_2 &= c\sin\theta_1+d\sin\theta_2.
\end{aligned}
```
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

# ╔═╡ 3f4bd706-7d31-5397-8e8f-99456df26f3a
md"""
**Try:** compare $\omega=0.95$, $1$, and $1.05$. Does an almost constant
phase necessarily mean a stable equilibrium? For the coupled phases, compare
the **unwrapped phase difference**: bounded variation is consistent with
1:1 locking; repeated $2\pi$ slips indicate drift. Two similar sine traces
over a short interval are not enough to establish locking.
"""

# ╔═╡ cee51649-f248-597d-8b0c-4918f6de5718
let s=sol2, ts=range(first(sol2.t),last(sol2.t);length=1201)
    delta=[s(t)[2]-s(t)[1] for t in ts]
    plot(ts,delta;xlabel="t",ylabel="θ₂ − θ₁ (unwrapped)",legend=false)
end

# ╔═╡ Cell order:
# ╠═96e6e962-71f3-4c7a-8c2d-ab7a6cc7baf4
# ╠═b15c712e-61bd-4568-8ba8-579212bb6e6c
# ╟─57f14c67-8e05-4c8a-90e5-6a4f3377506e
# ╟─92bbe726-f115-40fe-8a15-0c5bea8589df
# ╟─e298f3c0-ecc4-53e2-a709-310b68dca883
# ╠═b2ea050c-7f55-48ba-a165-41709ae5f78e
# ╠═99f2f5e8-f91d-4de9-b0b9-3e0d1c7b2ec9
# ╠═540a0af5-321f-41a4-8a90-756b798c84f2
# ╠═12f373d5-bc45-480b-93ba-d83b0c09f873
# ╠═416e0aed-a7bc-42e0-8a32-512471bc56f3
# ╠═55ed709b-bc02-4752-b026-8e38e6f5d17f
# ╟─796874fc-ebdd-508e-859a-8810320684e7
# ╠═9f854251-4be9-5469-ad55-34d89e30906c
# ╟─6a34a69d-9b35-5d79-b0a3-0b85a3c12e56
# ╠═147c066b-0ca4-5ef5-b793-08fc6561943e
# ╠═ca40de0e-805a-5b47-bd5c-edb18d4b7e05
# ╠═3a787114-4f7e-5cf8-b455-ad75e2e35c4a
# ╠═853c2184-d82c-5058-a2d2-15ed202205c2
# ╠═989c3787-fda2-5d26-a22d-ad1148fa5b86
# ╠═6cdfa56c-bd82-5726-98ad-2bac14a20221
# ╠═9cb434bb-7639-5753-a544-a0d2f2d500f3
# ╟─5ec39c31-afe9-41a0-8b4b-f6e68928144c
# ╟─85398bab-e2c9-565a-848f-ff022caec73d
# ╠═3169f185-9d52-4a80-95a0-dde68804ac84
# ╠═9e9dbfe7-668e-494a-8159-710ca914f2de
# ╟─2912533c-f785-4ca2-be22-5efeb5947567
# ╟─f1d13184-8e0c-41b0-af66-7e1925bf6b91
# ╠═da35c5c5-90ee-4335-82b7-49ae7bde66c7
# ╟─f520b8f5-cdfa-44ce-8dca-dcd8bf982111
# ╠═66fef05b-4acb-4f13-83a1-2c46e9ead2f8
# ╟─f67da2ec-7816-420a-96fe-aa6d72ba0aae
# ╟─5451b0b3-d3f1-51fd-88ee-cf2e7047a0c9
# ╠═a4fb3bc3-db31-4516-91b7-2c04abb169f5
# ╠═16968a95-d19f-4aa1-b3f4-7932971d9796
# ╠═a867d3fc-3373-4937-af30-fea2aab10f79
# ╟─3f4bd706-7d31-5397-8e8f-99456df26f3a
# ╠═cee51649-f248-597d-8b0c-4918f6de5718
