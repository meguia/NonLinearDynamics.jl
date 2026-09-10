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

# ╔═╡ 66270ae2-1483-11ee-2f94-9fa61a5bff52
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", "..", ".."))
    Pkg.instantiate()
    using ModelingToolkit, DifferentialEquations, Plots, PlutoUI
end

# ╔═╡ dd568e45-8431-4032-b535-2b440b6a8fc2
@variables t x1(t)=1 x2(t)=1 x3(t)=1 x4(t)=1

# ╔═╡ b9b0ee31-f055-40a4-85c3-109695e48f92
@parameters η1=1.85 η2=1.55 η3=0.77

# ╔═╡ fd5bb01b-6def-434a-a777-824c2e300e67
D = Differential(t)

# ╔═╡ 0ae113bf-b572-4bd9-9d69-96260ff880c8
function dhopf!(du,u,p,t)
	(η1,η2,η3) = p
	du[1] = η1*(0.5*u[1]+u[2]-u[4]-u[1]*(0.6*u[1]+u[1]*u[1]))
	du[2] = -η3*u[1]
	du[3] = (1+sqrt(2))*u[4]
	du[4] = (2-sqrt(2))*(u[1]-u[3]-η2*u[4])
end

# ╔═╡ 2001e25f-a361-40ab-aec6-b5a8bf4c8be0
eqs = [D(x1) ~ η1*(0.5*x1+x2-x4-x1*(0.6*x1+x1*x1))
	D(x2) ~ -η3*x1
	D(x3) ~ (1+sqrt(2))*x4
	D(x4) ~ (2-sqrt(2))*(x1-x3-η2*x4)]

# ╔═╡ 63c97c2a-d0c4-441b-940a-91011e7f1caa
@named sys = ODESystem(eqs, t)

# ╔═╡ 9794ee2f-8411-4679-8ac4-05770eb860b8
simpsys = structural_simplify(sys)

# ╔═╡ 9da4331a-7c0e-4a67-8bae-3594e798c07a
prob = ODEProblem(simpsys, [], (0,100.0))

# ╔═╡ c0051142-7135-4677-962f-26f6148739f7
sol_mtk = solve(prob, Tsit5());

# ╔═╡ 39a9a118-39da-45bb-baa9-6f767ad5e50a
let figure = plot(sol_mtk; idxs=(x1,x2), label="x1, x2")
    plot!(figure, sol_mtk; idxs=(x3,x4), label="x3, x4")
end

# ╔═╡ e76f9756-9804-4b7c-8fba-b6811a78084e
md"""
MTK builds an ODE from symbolic variables, parameters, and equations; the direct ODE example below retains its own initial conditions and parameter controls.

p1 $(@bind p1 Slider(1.7:0.001:2.5,default=1.85;show_value=true)) 
p2 $(@bind p2 Slider(1.4:0.001:2.0,default=1.55;show_value=true)) \
p3 $(@bind p3 Slider(0.5:0.001:1.0,default=0.75;show_value=true)) 
tmax $(@bind tmax Slider(100:10:1000,default=100;show_value=true)) \
"""

# ╔═╡ 39ccc3cb-7119-4302-b332-481bf9624ae0
prob1 = ODEProblem(dhopf!, [1.0,1.0,0.0,0.0], (0,tmax),[p1,p2,p3])

# ╔═╡ ebb5d99f-ef6d-49bc-9cca-e078cb331d95
sol = solve(prob1);

# ╔═╡ ba873076-94c8-48ca-8cb2-c6acae09195f
begin
	plot(sol,idxs=(1,2))
	plot!(sol,idxs=(3,4))
end	

# ╔═╡ Cell order:
# ╠═66270ae2-1483-11ee-2f94-9fa61a5bff52
# ╠═dd568e45-8431-4032-b535-2b440b6a8fc2
# ╠═b9b0ee31-f055-40a4-85c3-109695e48f92
# ╠═fd5bb01b-6def-434a-a777-824c2e300e67
# ╠═0ae113bf-b572-4bd9-9d69-96260ff880c8
# ╠═2001e25f-a361-40ab-aec6-b5a8bf4c8be0
# ╠═63c97c2a-d0c4-441b-940a-91011e7f1caa
# ╠═9794ee2f-8411-4679-8ac4-05770eb860b8
# ╠═9da4331a-7c0e-4a67-8bae-3594e798c07a
# ╠═c0051142-7135-4677-962f-26f6148739f7
# ╠═39a9a118-39da-45bb-baa9-6f767ad5e50a
# ╟─e76f9756-9804-4b7c-8fba-b6811a78084e
# ╠═39ccc3cb-7119-4302-b332-481bf9624ae0
# ╠═ebb5d99f-ef6d-49bc-9cca-e078cb331d95
# ╠═ba873076-94c8-48ca-8cb2-c6acae09195f
