### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 4df18ed1-628e-546a-a9f0-44a7912b76d0
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots
end

# ╔═╡ f1b40b34-ca27-5718-9166-0958242f6280
md"""
# Population growth with harvesting

A constant harvest lowers the logistic equilibrium; this example stops at extinction instead of continuing to negative population.
"""

# ╔═╡ e148d1f4-92be-5e85-a78a-3896e4825940
function model!(du, u, p, t)
    N = only(u)
    r, K, H = p
    du[1] = r*N*(1-N/K) - H
    nothing
end

# ╔═╡ a409ae6b-93fe-5b85-b2f3-28f6fa8c1aa3
u0 = [4.0]

# ╔═╡ 40136758-4be9-5dae-b930-250f5463c65a
tspan = (0.0, 20.0)

# ╔═╡ 053562d4-fdc6-5ee2-a144-e2ec5e3c6743
p = [1.0, 10.0, 1.5]  # r, K, H

# ╔═╡ 68f2e8f6-dd9d-5f89-a80e-19b8f46a8213
prob = ODEProblem(model!, u0, tspan, p)

# ╔═╡ 5778074b-0427-5c24-9266-728bd6b3c481
extinction = ContinuousCallback((u,t,integrator) -> u[1], nothing, terminate!)

# ╔═╡ 2cbba37a-6d8d-584e-b5f4-2b5b1411f8f0
sol = solve(prob, Tsit5(); abstol=1e-9, reltol=1e-7, callback=extinction, saveat=0.02);

# ╔═╡ eb8a3047-1699-554c-85fa-0d56c56e3595
plot(sol; idxs=[1], xlabel="t", ylabel="state")

# ╔═╡ Cell order:
# ╠═4df18ed1-628e-546a-a9f0-44a7912b76d0
# ╟─f1b40b34-ca27-5718-9166-0958242f6280
# ╠═e148d1f4-92be-5e85-a78a-3896e4825940
# ╠═a409ae6b-93fe-5b85-b2f3-28f6fa8c1aa3
# ╠═40136758-4be9-5dae-b930-250f5463c65a
# ╠═053562d4-fdc6-5ee2-a144-e2ec5e3c6743
# ╠═68f2e8f6-dd9d-5f89-a80e-19b8f46a8213
# ╠═5778074b-0427-5c24-9266-728bd6b3c481
# ╠═2cbba37a-6d8d-584e-b5f4-2b5b1411f8f0
# ╠═eb8a3047-1699-554c-85fa-0d56c56e3595
