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

# ╔═╡ 45047118-156c-11ee-2ef1-f3e884a0906f
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots, PlutoUI, ForwardDiff, StaticArrays, IntervalRootFinding
    using NonLinearDynamics
end

# ╔═╡ 0e01b152-bc20-5deb-ae42-3a6d6d6f4199
md"""
```math
\begin{aligned}
\dot{x} &= y-x[\eta_1+\gamma x(1+x)],\\
\dot{y} &= -x+y(\eta_2-2\gamma y^2).
\end{aligned}
```
"""

# ╔═╡ fdbb3f29-8de6-49b3-99e6-6b8aede86a8d
function ghopf!(du,u,p,t)
	(η1,η2,γ) = p
	du[1] = u[2]-u[1]*(η1+γ*u[1]*(1+u[1]))
	du[2] = -u[1]+u[2]*(η2-2*γ*u[2]*u[2])
end	

# ╔═╡ 095bd3fb-03ec-4f40-bd56-76c08e802a2a
u0_arr=[[xi,0] for xi=-1.01:0.1:-0.1]

# ╔═╡ 13d05f5e-76ee-4a77-a876-92cf1e708312
md"""
η1 $(@bind η1 Slider(0:0.001:2.0,default=0.85;show_value=true)) 
η2 $(@bind η2 Slider(0:0.001:2.0,default=0.55;show_value=true)) \
γ $(@bind γ Slider(0.002:0.001:0.1,default=0.01;show_value=true)) 
tmax $(@bind tmax Slider(10:10:1000,default=100;show_value=true)) \
"""

# ╔═╡ 69c13d15-4ee4-4961-9ec0-72825fbeef5c
prob = ODEProblem(ghopf!, [0.001,0.001], (0,tmax),[η1,η2,γ])

# ╔═╡ 16e6ce75-c342-492e-ba9d-2ae210453d42
ensamble_prob = EnsembleProblem(prob,prob_func=(prob,i,repeat;u0=u0_arr)->(remake(prob,u0=u0[i])))

# ╔═╡ a5087e7a-80e2-425a-b1c5-3348e90261da
sol = solve(ensamble_prob,EnsembleThreads(),trajectories=length(u0_arr));

# ╔═╡ a045f8a8-ab78-4aa3-b08f-4788df7c5b29
plot(sol,idxs=(1,2), xlim=(-5,5),ylim=(-5,5), linecolor=:black,linealpha=0.05)

# ╔═╡ Cell order:
# ╠═45047118-156c-11ee-2ef1-f3e884a0906f
# ╟─0e01b152-bc20-5deb-ae42-3a6d6d6f4199
# ╠═fdbb3f29-8de6-49b3-99e6-6b8aede86a8d
# ╠═69c13d15-4ee4-4961-9ec0-72825fbeef5c
# ╠═095bd3fb-03ec-4f40-bd56-76c08e802a2a
# ╠═16e6ce75-c342-492e-ba9d-2ae210453d42
# ╠═a5087e7a-80e2-425a-b1c5-3348e90261da
# ╠═a045f8a8-ab78-4aa3-b08f-4788df7c5b29
# ╠═13d05f5e-76ee-4a77-a876-92cf1e708312
