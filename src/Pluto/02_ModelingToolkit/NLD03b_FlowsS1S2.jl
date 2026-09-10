### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 2c28dbfe-830e-57b5-969f-ff37be88b17c
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots, ModelingToolkit
end

# ╔═╡ a976c390-96a7-5701-b457-1c15921426e2
md"""
# Coupled phases on a torus

Integrate the unwrapped Adler phases and plot their sines to see whether the oscillators lock.
"""

# ╔═╡ 14430aa5-bc16-5fe2-b91a-311abc18571a
# Variables
begin
    @independent_variables t
    @variables θ1(t) θ2(t)
end

# ╔═╡ 11141439-32b1-53cc-b1a8-b3a2cfa7d196
# Parameters
@parameters ω1 ω2 κ

# ╔═╡ 8023c427-6266-564c-9707-9ec7c78f8f9c
# Time derivative
D = Differential(t)

# ╔═╡ 938720fa-7889-51be-9863-64c022e3678d
equations = [
    D(θ1) ~ ω1 - sin(θ1) + κ*sin(θ2-θ1),
    D(θ2) ~ ω2 - sin(θ2) + κ*sin(θ1-θ2)
]

# ╔═╡ ba7f989f-d982-58e0-848d-f41745fb4c49
@named system = ODESystem(equations, t)

# ╔═╡ 2454f6f6-6c54-512d-9036-f32f02acdfd2
simplified = structural_simplify(system)

# ╔═╡ b1a123d1-75c3-56dc-be75-027fd57c872f
u0 = [θ1 => 0.0, θ2 => 0.2]

# ╔═╡ ab53c853-adc0-5bd2-a249-0687b2865261
tspan = (0.0, 80.0)

# ╔═╡ 3ca978e9-429e-598f-af0b-78891ad17cdc
p = [ω1 => 1.2, ω2 => 1.4, κ => 0.3]

# ╔═╡ feb46389-932b-5127-874d-f1c76473b830
prob = ODEProblem(simplified, u0, tspan, p)

# ╔═╡ c14d10a0-a05a-55ca-ab7c-20f69e405e59
sol = solve(prob, Tsit5(); abstol=1e-9, reltol=1e-7, saveat=0.02);

# ╔═╡ a8cae562-a347-5554-9e19-bf695f6a5ffb
plot(sol.t, [sin.(sol[θ1]) sin.(sol[θ2])];
    xlabel="t", ylabel="sin(θ)", label=["oscillator 1" "oscillator 2"])

# ╔═╡ Cell order:
# ╠═2c28dbfe-830e-57b5-969f-ff37be88b17c
# ╟─a976c390-96a7-5701-b457-1c15921426e2
# ╠═14430aa5-bc16-5fe2-b91a-311abc18571a
# ╠═11141439-32b1-53cc-b1a8-b3a2cfa7d196
# ╠═8023c427-6266-564c-9707-9ec7c78f8f9c
# ╠═938720fa-7889-51be-9863-64c022e3678d
# ╠═ba7f989f-d982-58e0-848d-f41745fb4c49
# ╠═2454f6f6-6c54-512d-9036-f32f02acdfd2
# ╠═b1a123d1-75c3-56dc-be75-027fd57c872f
# ╠═ab53c853-adc0-5bd2-a249-0687b2865261
# ╠═3ca978e9-429e-598f-af0b-78891ad17cdc
# ╠═feb46389-932b-5127-874d-f1c76473b830
# ╠═c14d10a0-a05a-55ca-ab7c-20f69e405e59
# ╠═a8cae562-a347-5554-9e19-bf695f6a5ffb
