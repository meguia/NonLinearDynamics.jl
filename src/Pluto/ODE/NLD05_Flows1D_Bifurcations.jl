### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ db4f14e7-32f9-5c87-93ee-f988cc5aacca
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots
end

# ╔═╡ 62518ca0-5495-5907-93cd-6535924c39f1
md"""
# One-dimensional bifurcations: pitchfork

Change μ through zero: the stable origin gives way to two stable equilibria in this supercritical pitchfork.
"""

# ╔═╡ 800d9c5c-bae0-5ed5-91bb-8fdc7fcade41
md"""
```math
\begin{aligned}
\dot{x} &= \mu x-x^3.
\end{aligned}
```
"""

# ╔═╡ 9bd171c2-0dbd-5768-91f4-a5023265fbb0
function model!(du, u, p, t)
    x = only(u)
    μ = only(p)
    du[1] = μ*x - x^3
    nothing
end

# ╔═╡ 7bb14eaf-4b73-5036-a893-afbbcbcc932e
u0 = [0.1]

# ╔═╡ d13525f4-c075-5f9b-a75d-e1a55ca1e303
tspan = (0.0, 30.0)

# ╔═╡ f00db1e4-6298-5fb1-b3fe-df8bc73d8365
p = [0.5]  # μ

# ╔═╡ ba0270a1-ef90-586c-8682-fd104d5b00d1
prob = ODEProblem(model!, u0, tspan, p)

# ╔═╡ 23ede884-9339-547d-8047-550719420738
sol = solve(prob, Tsit5(); abstol=1e-9, reltol=1e-7, saveat=0.02);

# ╔═╡ 50fca1fa-bcab-5504-af22-0302e55ea7df
plot(sol; idxs=[1], xlabel="t", ylabel="state")

# ╔═╡ 92cd25ba-520d-5fba-936b-5b37aa16bcc3
let negative=range(-1,0;length=201), positive=range(0,1;length=201)
    f=plot(negative,zero.(negative);label="stable origin",xlabel="μ",ylabel="equilibrium x",color=:steelblue)
    plot!(f,positive,zero.(positive);label="unstable origin",linestyle=:dash,color=:darkorange)
    plot!(f,positive,sqrt.(positive);label="stable branches",color=:steelblue)
    plot!(f,positive,-sqrt.(positive);label=false,color=:steelblue)
end

# ╔═╡ 5b21731b-e879-5abd-a1c2-ea34401fd517
md"""
## Three ways to change the number or stability of equilibria

A **bifurcation diagram** plots equilibria against a parameter, not against
time. Solid branches below are attracting; dashed branches are repelling.
The local tests are $f(x_*)=0$ and $f'(x_*)\lessgtr0$.

| Normal form | Equilibria | What happens at $\mu=0$? |
| --- | --- | --- |
| $\dot x=\mu-x^2$ | $x_*=\pm\sqrt\mu$ for $\mu>0$ | Saddle-node: a pair is created |
| $\dot x=\mu x-x^2$ | $x_*=0,\mu$ | Transcritical: branches exchange stability |
| $\dot x=\mu x-x^3$ | $0$; also $\pm\sqrt\mu$ for $\mu>0$ | Pitchfork: symmetry produces two attracting branches |

These normal forms describe a neighborhood of the bifurcation; their remote
trajectories need not remain bounded. At $\mu=0$ the linear test is inconclusive.
"""

# ╔═╡ dd90d3e7-ef6c-5183-a74a-cd458b7dede1
function course_normal_form!(du,u,p,t)
    μ,kind=p
    x=u[1]
    du[1] = kind==1 ? μ-x^2 : kind==2 ? μ*x-x^2 : μ*x-x^3
    nothing
end

# ╔═╡ 0309fd38-15e0-55cf-9b37-03abbf5994a0
course_kind = 3  # 1: saddle-node; 2: transcritical; 3: pitchfork

# ╔═╡ 4c0724a8-8ba2-5b14-8428-71f1743eb66d
course_mu = 0.2  # μ

# ╔═╡ fd5b4121-8d0a-5be7-99da-d8deb0ec7e4e
course_normal_u0 = [0.15]

# ╔═╡ a27e6e38-1be3-5f7b-8382-cae5c24ecdc4
course_normal_tspan = (0.0,25.0)

# ╔═╡ 295f9e97-c69b-5ec8-a629-4e2ef9b4fd34
course_normal_p = [course_mu,course_kind]

# ╔═╡ bcec42e1-c96e-59af-a663-7c8cee43d9f6
course_normal_prob = ODEProblem(course_normal_form!, course_normal_u0, course_normal_tspan, course_normal_p)

# ╔═╡ 2edfe101-589e-57e7-aea3-92d1352ee0af
course_normal_sol = solve(course_normal_prob, Tsit5(); abstol=1e-9, reltol=1e-8, callback=ContinuousCallback((u,t,i)->4-abs(u[1]),nothing,terminate!));

# ╔═╡ cb80123b-8743-563a-b36d-24fd5e571b49
let xs=range(-1.2,1.2;length=401), dn=zeros(1)
    f=[(course_normal_form!(dn,[x],course_normal_p,0.0);dn[1]) for x in xs]
    field=plot(xs,f;label="f(x)",xlabel="x",ylabel="dx/dt")
    hline!(field,[0];color=:gray,label=false)
    trace=plot(course_normal_sol;label="x(t)",xlabel="t")
    plot(field,trace;layout=(1,2),size=(900,370),margin=5*Plots.mm)
end

# ╔═╡ 726a770d-f00a-590a-9fde-fad0eb8789e3
let positive=range(0,0.5;length=201), negative=range(-0.5,0;length=201)
    figures=Plots.Plot[]
    for (kind,title) in enumerate(("Saddle-node","Transcritical","Pitchfork"))
        f=plot(;xlabel="μ",ylabel="x*",title,legend=false)
        if kind==1
            plot!(f,positive,sqrt.(positive);color=:steelblue)
            plot!(f,positive,-sqrt.(positive);color=:darkorange,linestyle=:dash)
        else
            plot!(f,negative,zero.(negative);color=:steelblue)
            plot!(f,positive,zero.(positive);color=:darkorange,linestyle=:dash)
            if kind==2
                plot!(f,positive,positive;color=:steelblue)
                plot!(f,negative,negative;color=:darkorange,linestyle=:dash)
            else
                for sign in (-1,1)
                    plot!(f,positive,sign.*sqrt.(positive);color=:steelblue)
                end
            end
        end
        scatter!(f,[0],[0];color=:black,markersize=3)
        push!(figures,f)
    end
    plot(figures...;layout=(1,3),size=(1000,340),margin=5*Plots.mm)
end

# ╔═╡ 60a15fb2-50d3-5761-b73f-316a7fc8d19d
md"""
**Try:** move $\mu$ through zero, then start on the other side of an equilibrium.
Which branch attracts, and which separates different destinations? For the
saddle-node just below zero, look for slow passage through the **ghost** of the
vanished pair. The example stops at $|x|=4$ to keep this local experiment finite.
"""

# ╔═╡ Cell order:
# ╠═db4f14e7-32f9-5c87-93ee-f988cc5aacca
# ╟─62518ca0-5495-5907-93cd-6535924c39f1
# ╟─800d9c5c-bae0-5ed5-91bb-8fdc7fcade41
# ╠═9bd171c2-0dbd-5768-91f4-a5023265fbb0
# ╠═7bb14eaf-4b73-5036-a893-afbbcbcc932e
# ╠═d13525f4-c075-5f9b-a75d-e1a55ca1e303
# ╠═f00db1e4-6298-5fb1-b3fe-df8bc73d8365
# ╠═ba0270a1-ef90-586c-8682-fd104d5b00d1
# ╠═23ede884-9339-547d-8047-550719420738
# ╠═50fca1fa-bcab-5504-af22-0302e55ea7df
# ╠═92cd25ba-520d-5fba-936b-5b37aa16bcc3
# ╟─5b21731b-e879-5abd-a1c2-ea34401fd517
# ╠═dd90d3e7-ef6c-5183-a74a-cd458b7dede1
# ╠═0309fd38-15e0-55cf-9b37-03abbf5994a0
# ╠═4c0724a8-8ba2-5b14-8428-71f1743eb66d
# ╠═fd5b4121-8d0a-5be7-99da-d8deb0ec7e4e
# ╠═a27e6e38-1be3-5f7b-8382-cae5c24ecdc4
# ╠═295f9e97-c69b-5ec8-a629-4e2ef9b4fd34
# ╠═bcec42e1-c96e-59af-a663-7c8cee43d9f6
# ╠═2edfe101-589e-57e7-aea3-92d1352ee0af
# ╠═cb80123b-8743-563a-b36d-24fd5e571b49
# ╠═726a770d-f00a-590a-9fde-fad0eb8789e3
# ╟─60a15fb2-50d3-5761-b73f-316a7fc8d19d
