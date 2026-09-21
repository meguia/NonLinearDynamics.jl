### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ ff4ccfc6-c478-5eb8-ad55-6a263f056929
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using NonLinearDynamics: flow2d_nullclines
    using DifferentialEquations, Plots
    using LinearAlgebra
end

# ╔═╡ cadba734-5b4b-510b-abd1-bddc2792cbec
md"""
# The Duffing oscillator

A double-well potential gives two stable resting positions; change the initial state to explore their basins.
"""

# ╔═╡ 68460e98-3998-537c-8a5f-be8c7b6b1331
md"""
```math
\begin{aligned}
\dot{x} &= v, & \dot{v} &= -\gamma v+\beta x-x^3.
\end{aligned}
```
"""

# ╔═╡ 4bc05396-18f8-59fb-be29-1f988fa376a9
function model!(du, u, p, t)
    x, v = u
    γ, β = p
    du[1] = v
    du[2] = -γ*v + β*x - x^3
    nothing
end

# ╔═╡ b9eb2d07-4791-520a-a3c8-6c2b1fbf55a4
u0 = [0.1, 0.7]

# ╔═╡ a1df769d-1131-5b35-86cc-2c969cb2d19c
tspan = (0.0, 60.0)

# ╔═╡ d1ba8bea-1d62-51f2-8044-4cf10095fae1
p = [0.15, 1.0]  # γ, β

# ╔═╡ c67fa65d-32d0-5cef-a854-33de4d66d906
prob = ODEProblem(model!, u0, tspan, p)

# ╔═╡ 9a436f19-7c6a-5fe4-a884-4d0e268b0f8b
sol = solve(prob, Tsit5(); abstol=1e-9, reltol=1e-7, saveat=0.02);

# ╔═╡ 3125de08-bb5b-5dee-b608-cda49b83d138
plot(sol; idxs=[1, 2], xlabel="t", ylabel="state")

# ╔═╡ a8bb9eb9-8562-5c7f-8656-ed6ab53c0dd3
plot(sol; idxs=(1, 2), legend=false)

# ╔═╡ 158b8eb6-a0a3-5827-a0e4-9222abe430a6
md"""
## From nullclines to local stability and energy

Solve both nullcline equations, then evaluate the Jacobian at each intersection:

```math
v=0,\quad x(\beta-x^2)=0,\qquad
J(x_*,0)=\begin{pmatrix}0&1\\\beta-3x_*^2&-\gamma\end{pmatrix}.
```

For $\beta>0$ the origin is a saddle and the outer equilibria are attracting
when $\gamma>0$. The saddle's stable manifold separates their basins.
The energy gives complementary global information:

```math
V(x)=\frac{x^4}{4}-\frac{\beta x^2}{2},\qquad
E=\frac{v^2}{2}+V(x),\qquad \dot E=-\gamma v^2.
```

At $\gamma=0$ there are centers and conserved-energy orbits, rather than
attractors. A double-well potential alone does not imply attraction.
"""

# ╔═╡ b0e810f4-a0ea-5a31-9752-bd313f4d2a34
let parameters=p, xs=range(-1.6,1.6;length=401)
    gamma,beta=parameters
    a=plot(xs,xs.^4 ./4 .-beta.*xs.^2 ./2;xlabel="x",ylabel="V(x)",legend=false)
    b=flow2d_nullclines(model!,parameters;xlims=[-1.6,1.6],ylims=[-1.5,1.5],
        regions=false,vectorfield=true,title="red: dx = 0; green: dv = 0")
    equilibria=beta>0 ? [-sqrt(beta),0,sqrt(beta)] : [0.0]
    scatter!(a,equilibria,equilibria.^4 ./4 .-beta.*equilibria.^2 ./2;color=:black)
    scatter!(b,equilibria,zero.(equilibria);color=:black)
    plot(a,b;layout=(1,2),size=(950,390),margin=5*Plots.mm)
end

# ╔═╡ a9d7edb0-3292-53ba-a33d-a16600dbb114
let
    gamma,beta=p
    equilibria=beta>0 ? [-sqrt(beta),0,sqrt(beta)] : [0.0]
    [(x=xstar, jacobian=[0.0 1.0;beta-3xstar^2 -gamma],
      eigenvalues=eigvals([0.0 1.0;beta-3xstar^2 -gamma])) for xstar in equilibria]
end

# ╔═╡ 57b93009-328d-507c-be88-11a5813edd7b
md"""
## A biological phase plane: Lotka–Volterra

Let $x>0$ be prey density and $y>0$ predator density. Encounters proportional
to $xy$ reduce prey and support predator reproduction:

```math
\dot x=(b-cy)x,\qquad \dot y=(ax-d)y,\qquad a,b,c,d>0.
```

The interior equilibrium is $(d/a,b/c)$. Its imaginary linear eigenvalues
alone do not establish a nonlinear center, but here the conserved quantity

```math
H(x,y)=ax-d\log x+cy-b\log y
```

does: different initial conditions remain on different closed level curves.
These persistent oscillations are **not** an attracting limit cycle.
"""

# ╔═╡ 8445ea32-88c9-5050-83cb-17f2d16254e5
function course_lv!(du,u,p,t)
    a,b,c,d=p
    x,y=u
    du[1]=(b-c*y)*x
    du[2]=(a*x-d)*y
    nothing
end

# ╔═╡ 0276de59-d82c-5da5-8b4d-e6ab1da20493
course_prey0 = 0.6  # initial prey density

# ╔═╡ 01cc4ff8-44f2-5cf3-923e-89c9a263898b
course_lv_u0 = [course_prey0,0.6]

# ╔═╡ 25529717-678f-560e-9544-827f39a1b036
course_lv_tspan = (0.0,40.0)

# ╔═╡ 7bf2512c-fce4-50bb-aba3-cb987938c2b2
course_lv_p = [1.0,1.0,1.0,1.0]

# ╔═╡ 82375086-388d-590c-a3d2-ce47d9f01902
course_lv_prob = ODEProblem(course_lv!, course_lv_u0, course_lv_tspan, course_lv_p)

# ╔═╡ 6d95ace9-4784-52b8-b049-46f8603c117d
course_lv_sol = solve(course_lv_prob, Tsit5(); abstol=1e-9, reltol=1e-8);

# ╔═╡ 8e14a3d3-22bd-5750-9bbf-43d04ca4efc8
let
    a,b,c,d=course_lv_p
    f=flow2d_nullclines(course_lv!,course_lv_p;xlims=[0.05,3.5],ylims=[0.05,3.5],
        regions=false,title="Lotka–Volterra: center orbits")
    plot!(f,course_lv_sol;idxs=(1,2),color=:black)
    for initial in ([0.8,0.8],[1.5,1.5])
        comparison=solve(remake(course_lv_prob;u0=initial),Tsit5();abstol=1e-9,reltol=1e-8)
        plot!(f,comparison;idxs=(1,2))
    end
    scatter!(f,[d/a],[b/c];color=:black)
end

# ╔═╡ 0512ed3c-4a8d-5f6a-9c57-a2e5130a2d70
let
    a,b,c,d=course_lv_p
    x=course_lv_sol[1,:]
    y=course_lv_sol[2,:]
    H=a.*x .-d.*log.(x) .+c.*y .-b.*log.(y)
    plot(course_lv_sol.t,H .-first(H);xlabel="t",ylabel="H(t) − H(0)",
        legend=false,title="Numerical drift of an exact invariant")
end

# ╔═╡ a753f229-51b3-5968-937a-6c90a21a252a
md"""
## A threshold without a limit cycle: reduced SIR

In a closed population, susceptible ($S$), infected ($I$), and removed ($R$)
fractions satisfy $S+I+R=1$. Eliminating $R$ gives the two-dimensional flow

```math
\dot S=-\beta SI,\qquad \dot I=\beta SI-\gamma I,\qquad
R=1-S-I,\quad \dot R=\gamma I.
```

For $I>0$, infections grow only while $S>\gamma/\beta$. The trajectory crosses
that nullcline at the peak of $I(t)$; it cannot form a closed orbit because
$S$ decreases. The whole line $I=0$ consists of equilibria, so a zero
eigenvalue there is expected. This is a simple mathematical epidemic model.
"""

# ╔═╡ 601d71eb-dc1b-5b7f-8620-b6547790299c
function course_sir!(du,u,p,t)
    β,γ=p
    S,I=u
    du[1]=-β*S*I
    du[2]=β*S*I-γ*I
    nothing
end

# ╔═╡ e88b1bb5-1dc2-58cd-a16f-7a5013fdd4a6
course_infection = 2.0  # β (γ = 1)

# ╔═╡ 052950bd-0f44-58de-b372-cbf8cad7360d
course_sir_u0 = [0.99,0.01]

# ╔═╡ abe95819-b95f-5b62-aad7-a61ff34ba3de
course_sir_tspan = (0.0,20.0)

# ╔═╡ 6507ea32-e353-5fdf-ba49-7ca80797f9fc
course_sir_p = [course_infection,1.0]

# ╔═╡ 4e6013c1-6102-50e2-a15a-23e028f8a176
course_sir_prob = ODEProblem(course_sir!, course_sir_u0, course_sir_tspan, course_sir_p)

# ╔═╡ 638ca57b-6c88-51a4-8b47-67ee195575e9
course_sir_sol = solve(course_sir_prob, Tsit5(); abstol=1e-9, reltol=1e-8);

# ╔═╡ 8b7b394a-7de8-53ae-9e49-79ac97affbfb
let S=course_sir_sol[1,:], I=course_sir_sol[2,:]
    a=plot(course_sir_sol.t,hcat(S,I,1 .-S .-I);label=["S" "I" "R"],xlabel="t",ylabel="fraction")
    b=plot(S,I;xlabel="S",ylabel="I",legend=false,xlims=(0,1),ylims=(0,0.6))
    threshold=course_sir_p[2]/course_sir_p[1]
    threshold<=1 && vline!(b,[threshold];linestyle=:dash,color=:green)
    plot!(b,[0,1],[1,0];color=:gray,label=false)
    plot(a,b;layout=(1,2),size=(900,370),margin=5*Plots.mm)
end

# ╔═╡ 67a3efe1-e8cf-5e32-8606-4f8706a72544
md"""
**Compare:** why does a perturbation change the Lotka–Volterra orbit, while
nearby trajectories of a self-oscillator recover a preferred amplitude?
Why does the SIR outbreak stop growing while susceptible individuals remain?
Use nullclines, invariants, and monotonic variables before relying on a time plot.
"""

# ╔═╡ Cell order:
# ╠═ff4ccfc6-c478-5eb8-ad55-6a263f056929
# ╟─cadba734-5b4b-510b-abd1-bddc2792cbec
# ╟─68460e98-3998-537c-8a5f-be8c7b6b1331
# ╠═4bc05396-18f8-59fb-be29-1f988fa376a9
# ╠═b9eb2d07-4791-520a-a3c8-6c2b1fbf55a4
# ╠═a1df769d-1131-5b35-86cc-2c969cb2d19c
# ╠═d1ba8bea-1d62-51f2-8044-4cf10095fae1
# ╠═c67fa65d-32d0-5cef-a854-33de4d66d906
# ╠═9a436f19-7c6a-5fe4-a884-4d0e268b0f8b
# ╠═3125de08-bb5b-5dee-b608-cda49b83d138
# ╠═a8bb9eb9-8562-5c7f-8656-ed6ab53c0dd3
# ╟─158b8eb6-a0a3-5827-a0e4-9222abe430a6
# ╠═b0e810f4-a0ea-5a31-9752-bd313f4d2a34
# ╠═a9d7edb0-3292-53ba-a33d-a16600dbb114
# ╟─57b93009-328d-507c-be88-11a5813edd7b
# ╠═8445ea32-88c9-5050-83cb-17f2d16254e5
# ╠═0276de59-d82c-5da5-8b4d-e6ab1da20493
# ╠═01cc4ff8-44f2-5cf3-923e-89c9a263898b
# ╠═25529717-678f-560e-9544-827f39a1b036
# ╠═7bf2512c-fce4-50bb-aba3-cb987938c2b2
# ╠═82375086-388d-590c-a3d2-ce47d9f01902
# ╠═6d95ace9-4784-52b8-b049-46f8603c117d
# ╠═8e14a3d3-22bd-5750-9bbf-43d04ca4efc8
# ╠═0512ed3c-4a8d-5f6a-9c57-a2e5130a2d70
# ╟─a753f229-51b3-5968-937a-6c90a21a252a
# ╠═601d71eb-dc1b-5b7f-8620-b6547790299c
# ╠═e88b1bb5-1dc2-58cd-a16f-7a5013fdd4a6
# ╠═052950bd-0f44-58de-b372-cbf8cad7360d
# ╠═abe95819-b95f-5b62-aad7-a61ff34ba3de
# ╠═6507ea32-e353-5fdf-ba49-7ca80797f9fc
# ╠═4e6013c1-6102-50e2-a15a-23e028f8a176
# ╠═638ca57b-6c88-51a4-8b47-67ee195575e9
# ╠═8b7b394a-7de8-53ae-9e49-79ac97affbfb
# ╟─67a3efe1-e8cf-5e32-8606-4f8706a72544
