### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 21d57702-f5af-54c3-8a03-3a125964037c
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots
end

# ╔═╡ fa128836-5ed9-57bd-861f-4410bb72451c
md"""
# Love affairs: a linear system

The coefficients a and d describe each person’s own response, while b and c describe their response to the partner.
"""

# ╔═╡ b562a1ef-d4d2-5cee-bb93-53051819bdbb
md"""
```math
\begin{aligned}
\dot{R} &= aR+bJ, & \dot{J} &= cR+dJ.
\end{aligned}
```
"""

# ╔═╡ cb1ff608-36c5-5617-8a31-c1cba8c61420
function model!(du, u, p, t)
    R, J = u
    a, b, c, d = p
    du[1] = a*R + b*J
    du[2] = c*R + d*J
    nothing
end

# ╔═╡ 27c5e722-eb3f-5ebe-b287-30a62bf04bb9
u0 = [1.0, 0.0]

# ╔═╡ 73d6777e-36ca-517e-acc6-454aa565a23d
tspan = (0.0, 30.0)

# ╔═╡ 9e65dea9-58b6-5d8c-83ac-9eca424cc319
p = [-0.1, 1.0, -1.0, -0.1]  # a, b, c, d

# ╔═╡ 8a8d9145-160a-531f-9879-51e8c55709c2
prob = ODEProblem(model!, u0, tspan, p)

# ╔═╡ c4931c16-b627-5872-aaa1-ee01641c2753
sol = solve(prob, Tsit5(); abstol=1e-9, reltol=1e-7, saveat=0.02);

# ╔═╡ d7a5ebb9-1a1a-537b-9fe8-535098a4b603
plot(sol; idxs=[1, 2], xlabel="t", ylabel="state")

# ╔═╡ 6f716367-a116-5483-9b0c-830c07bc295b
plot(sol; idxs=(1, 2), legend=false)

# ╔═╡ Cell order:
# ╠═21d57702-f5af-54c3-8a03-3a125964037c
# ╟─fa128836-5ed9-57bd-861f-4410bb72451c
# ╟─b562a1ef-d4d2-5cee-bb93-53051819bdbb
# ╠═cb1ff608-36c5-5617-8a31-c1cba8c61420
# ╠═27c5e722-eb3f-5ebe-b287-30a62bf04bb9
# ╠═73d6777e-36ca-517e-acc6-454aa565a23d
# ╠═9e65dea9-58b6-5d8c-83ac-9eca424cc319
# ╠═8a8d9145-160a-531f-9879-51e8c55709c2
# ╠═c4931c16-b627-5872-aaa1-ee01641c2753
# ╠═d7a5ebb9-1a1a-537b-9fe8-535098a4b603
# ╠═6f716367-a116-5483-9b0c-830c07bc295b
