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
    using LinearAlgebra
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

# ╔═╡ 4d1ecd5b-6435-54bd-9a97-901b88d07a5f
md"""
## Equilibria organize the three-dimensional flow

In the convection interpretation, $x$ measures circulation, $y$ a horizontal
temperature contrast, and $z$ a vertical temperature-profile correction.
The three-mode truncation is an instructive dynamical system, not a complete
weather model.

For $\sigma,\beta>0$, the origin is stable for $\rho<1$. For $\rho>1$,
two equilibria appear:

```math
C_\pm=(\pm\sqrt{\beta(\rho-1)},\ \pm\sqrt{\beta(\rho-1)},\ \rho-1),\qquad
J=\begin{pmatrix}
-\sigma&\sigma&0\\ \rho-z&-1&-x\\ y&x&-\beta
\end{pmatrix}.
```

For $\sigma>\beta+1$, their linear Hopf threshold is
$\rho_H=\sigma(\sigma+\beta+3)/(\sigma-\beta-1)$, about $24.74$ for
$\sigma=10$, $\beta=8/3$. This local threshold does **not** locate all
global transitions or guarantee chaos at every larger $\rho$.
"""

# ╔═╡ b03b4e9a-db1a-5658-8fa7-99f2e5f4bee0
let
    σ,ρ,β=p
    points=[[0.0,0.0,0.0]]
    if ρ>1 && β>0
        q=sqrt(β*(ρ-1))
        append!(points,[[q,q,ρ-1],[-q,-q,ρ-1]])
    end
    [(point=u,eigenvalues=eigvals([-σ σ 0;ρ-u[3] -1 -u[1];u[2] u[1] -β])) for u in points]
end

# ╔═╡ f6449519-09b6-5c2e-80f8-f6d3b92e7b6a
md"""
## Nearby initial states and the prediction horizon

Integrate the same equations from $\mathbf u_0$ and
$\mathbf u_0+(10^{-7},0,0)$ with the same tight tolerances. Compare the
full-state separation $d(t)=\|\mathbf u_1(t)-\mathbf u_2(t)\|$ on a log scale.
In a chaotic regime, an initial growth interval is followed by saturation at
the attractor's scale. A finite-time separation curve is not by itself a
Lyapunov-exponent calculation; also check tolerance sensitivity.
"""

# ╔═╡ 203fd202-777e-5d9c-aa9c-47d5744ad198
course_lorenz_u0 = [1.0,0.0,0.0]

# ╔═╡ 9830d2a4-1d52-5086-b5ab-6e8aae69787c
course_lorenz_tspan = (0.0,60.0)

# ╔═╡ 184b94b9-e3c9-5fc1-abb1-c2cf65894a51
course_lorenz_p = p

# ╔═╡ 447cdaa6-aa3c-55fc-a206-d602720af57a
course_lorenz_prob = ODEProblem(model!, course_lorenz_u0, course_lorenz_tspan, course_lorenz_p)

# ╔═╡ 56d9eae6-4545-53d0-9ca7-7738a5f0e181
course_lorenz_sol = solve(course_lorenz_prob, Tsit5(); abstol=1e-9, reltol=1e-8);

# ╔═╡ 8df0c9cb-c516-5565-a977-912fe6032c4a
course_lorenz_near = solve(remake(course_lorenz_prob;u0=course_lorenz_u0+[1e-7,0,0]),
    Tsit5();abstol=1e-9,reltol=1e-8);

# ╔═╡ ef1193aa-0828-5e38-8e72-6cb52f2cca0e
let ts=range(course_lorenz_tspan...;length=3001)
    distance=[norm(course_lorenz_sol(t)-course_lorenz_near(t)) for t in ts]
    a=plot(course_lorenz_sol;idxs=(0,1),label="reference",xlabel="t",ylabel="x")
    plot!(a,course_lorenz_near;idxs=(0,1),label="perturbed",alpha=0.7)
    b=plot(ts,max.(distance,eps());yscale=:log10,xlabel="t",ylabel="‖δu(t)‖",legend=false)
    plot(a,b;layout=(2,1),size=(900,500),margin=5*Plots.mm)
end

# ╔═╡ 8c8dc7e9-eb70-527e-8af0-ac6df18a539e
md"""
## Reduce a flow to a return map

Record successive maxima of $z(t)$ after the transient. A maximum occurs
when $\dot z=xy-\beta z$ crosses zero **from positive to negative**.
The callback locates that event between solver steps; sampling a coarse
time grid can miss it. Plot $z_{n+1}$ against $z_n$.

At the classical chaotic parameters this projection gives an almost
one-dimensional folded relation. It illustrates stretching and folding,
but is a projection of a return map rather than an exact scalar evolution law.
"""

# ╔═╡ 419228a4-01b9-5dfd-8c23-b532a73dc741
course_lorenz_maxima = let times=Float64[], maxima=Float64[]
    callback=ContinuousCallback(
        (u,t,i)->u[1]*u[2]-i.p[3]*u[3],nothing,
        i->begin
            if i.t>30
                push!(times,i.t)
                push!(maxima,i.u[3])
            end
        end;save_positions=(false,false))
    solution=solve(remake(course_lorenz_prob;tspan=(0.0,200.0)),Tsit5();
        callback,abstol=1e-9,reltol=1e-8,dtmax=0.05,save_everystep=false)
    (times=times,values=maxima,retcode=solution.retcode)
end

# ╔═╡ 85242f18-97b8-5e90-948b-4e426831aef8
let z=course_lorenz_maxima.values
    if length(z)>=3 && maximum(z)-minimum(z)>1e-5
        f=scatter(z[1:end-1],z[2:end];xlabel="zₙ (maximum)",ylabel="zₙ₊₁",
            markersize=2,markerstrokewidth=0,legend=false,size=(550,450),margin=5*Plots.mm)
        plot!(f,collect(extrema(z)),collect(extrema(z));color=:gray,linestyle=:dash)
    else
        md"No resolved sequence of distinct maxima at these parameters; try σ = 10, ρ = 28, β = 8/3."
    end
end

# ╔═╡ 054125f3-6639-57bf-a49e-adaaef332854
md"""
**Try:** compare $\rho=0.5$, $10$, and $28$ with $\sigma=10$, $\beta=8/3$.
Explain the final behavior using equilibria, the separation curve, and the
maxima map together. Near bifurcations, extend the transient before deciding
that irregular motion is a persistent attractor.
"""

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
# ╟─4d1ecd5b-6435-54bd-9a97-901b88d07a5f
# ╠═b03b4e9a-db1a-5658-8fa7-99f2e5f4bee0
# ╟─f6449519-09b6-5c2e-80f8-f6d3b92e7b6a
# ╠═203fd202-777e-5d9c-aa9c-47d5744ad198
# ╠═9830d2a4-1d52-5086-b5ab-6e8aae69787c
# ╠═184b94b9-e3c9-5fc1-abb1-c2cf65894a51
# ╠═447cdaa6-aa3c-55fc-a206-d602720af57a
# ╠═56d9eae6-4545-53d0-9ca7-7738a5f0e181
# ╠═8df0c9cb-c516-5565-a977-912fe6032c4a
# ╠═ef1193aa-0828-5e38-8e72-6cb52f2cca0e
# ╟─8c8dc7e9-eb70-527e-8af0-ac6df18a539e
# ╠═419228a4-01b9-5dfd-8c23-b532a73dc741
# ╠═85242f18-97b8-5e90-948b-4e426831aef8
# ╟─054125f3-6639-57bf-a49e-adaaef332854
