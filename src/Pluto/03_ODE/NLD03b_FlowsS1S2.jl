### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 03c69f7b-1eeb-585d-a787-5f968b7eaa08
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots
end

# ╔═╡ 0e026a8c-edbd-521c-938f-421f82bb3ec1
md"""
# Coupled phases on a torus

Integrate the unwrapped Adler phases and plot their sines to see whether the oscillators lock.
"""

# ╔═╡ 1ffd3f2d-1e4f-5b37-8536-84b15e1277db
function model!(du, u, p, t)
    θ1, θ2 = u
    ω1, ω2, κ = p
    du[1] = ω1 - sin(θ1) + κ*sin(θ2-θ1)
    du[2] = ω2 - sin(θ2) + κ*sin(θ1-θ2)
    nothing
end

# ╔═╡ e45a6c49-f2f3-5598-a4d8-70c4504da728
u0 = [0.0, 0.2]

# ╔═╡ 1294d097-6879-58b3-98e8-20c5d7588479
tspan = (0.0, 80.0)

# ╔═╡ d155dc43-216e-59e4-9a14-6ab09f5ccc31
p = [1.2, 1.4, 0.3]  # ω1, ω2, κ

# ╔═╡ 2935661f-c3b9-5b21-b878-3ceca93d1491
prob = ODEProblem(model!, u0, tspan, p)

# ╔═╡ 42aba460-4bc1-508b-be55-db624b15da0d
sol = solve(prob, Tsit5(); abstol=1e-9, reltol=1e-7, saveat=0.02);

# ╔═╡ 0bbc70dd-3183-52b2-b5d5-d0438d473b56
plot(sol.t, [sin.(sol[1,:]) sin.(sol[2,:])];
    xlabel="t", ylabel="sin(θ)", label=["oscillator 1" "oscillator 2"])

# ╔═╡ Cell order:
# ╠═03c69f7b-1eeb-585d-a787-5f968b7eaa08
# ╟─0e026a8c-edbd-521c-938f-421f82bb3ec1
# ╠═1ffd3f2d-1e4f-5b37-8536-84b15e1277db
# ╠═e45a6c49-f2f3-5598-a4d8-70c4504da728
# ╠═1294d097-6879-58b3-98e8-20c5d7588479
# ╠═d155dc43-216e-59e4-9a14-6ab09f5ccc31
# ╠═2935661f-c3b9-5b21-b878-3ceca93d1491
# ╠═42aba460-4bc1-508b-be55-db624b15da0d
# ╠═0bbc70dd-3183-52b2-b5d5-d0438d473b56
