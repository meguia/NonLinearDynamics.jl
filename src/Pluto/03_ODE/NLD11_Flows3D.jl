### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ c01ad9e6-5574-520e-8444-ba5c873d4990
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots
end

# ╔═╡ da649440-5e10-546f-ad57-4f9025ad9247
md"""
# Flows in 3D: Lorenz system

Three coupled equations produce the Lorenz attractor; short-time agreement is meaningful even when long chaotic trajectories separate.
"""

# ╔═╡ 59a9ae19-20f8-5ad4-8777-09f0be32063a
md"""
```math
\begin{aligned}
\dot{x} &= \sigma(y-x),\\
\dot{y} &= x(\rho-z)-y,\\
\dot{z} &= xy-\beta z.
\end{aligned}
```
"""

# ╔═╡ 122c348a-4bf4-5c4f-8f33-caabe781eb3d
function model!(du, u, p, t)
    x, y, z = u
    σ, ρ, β = p
    du[1] = σ*(y-x)
    du[2] = x*(ρ-z)-y
    du[3] = x*y-β*z
    nothing
end

# ╔═╡ cef048c5-9987-56a9-aa15-7ebe72f2698c
u0 = [1.0, 0.0, 0.0]

# ╔═╡ 89f85229-c35d-5675-bbfd-0d476ad033c2
tspan = (0.0, 100.0)

# ╔═╡ a9c202d5-e0c4-5806-bb88-a9fd5d607d79
p = [10.0, 28.0, 8/3]  # σ, ρ, β

# ╔═╡ a526e38f-7d08-5b1e-b6f2-54c3cf416225
prob = ODEProblem(model!, u0, tspan, p)

# ╔═╡ 519c061b-820b-5669-b5ae-b2beb5b86f3b
sol = solve(prob, Tsit5(); abstol=1e-9, reltol=1e-7, saveat=0.01);

# ╔═╡ 10397b16-3994-5868-bded-88280657ddec
plot(sol; idxs=[1, 2, 3], xlabel="t", ylabel="state")

# ╔═╡ a6011ee8-5b75-503b-9d9e-2df15082e739
plot(sol; idxs=(1, 2, 3), legend=false)

# ╔═╡ Cell order:
# ╠═c01ad9e6-5574-520e-8444-ba5c873d4990
# ╟─da649440-5e10-546f-ad57-4f9025ad9247
# ╟─59a9ae19-20f8-5ad4-8777-09f0be32063a
# ╠═122c348a-4bf4-5c4f-8f33-caabe781eb3d
# ╠═cef048c5-9987-56a9-aa15-7ebe72f2698c
# ╠═89f85229-c35d-5675-bbfd-0d476ad033c2
# ╠═a9c202d5-e0c4-5806-bb88-a9fd5d607d79
# ╠═a526e38f-7d08-5b1e-b6f2-54c3cf416225
# ╠═519c061b-820b-5669-b5ae-b2beb5b86f3b
# ╠═10397b16-3994-5868-bded-88280657ddec
# ╠═a6011ee8-5b75-503b-9d9e-2df15082e739
