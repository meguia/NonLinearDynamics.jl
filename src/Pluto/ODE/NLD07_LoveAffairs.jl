### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 21d57702-f5af-54c3-8a03-3a125964037c
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using NonLinearDynamics: classification_linear
    using DifferentialEquations, Plots
    using LinearAlgebra
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

# ╔═╡ 79a3579c-b488-5854-a2c6-d7c46efdabc9
md"""
## Linear flows: eigenvalues and geometry

For a constant matrix $A$, the equations and their characteristic polynomial are

```math
\dot{\mathbf u}=A\mathbf u,\qquad
\lambda^2-\tau\lambda+\Delta=0,\quad
\tau=\operatorname{tr}A,\quad\Delta=\det A.
```

For $\Delta<0$ the origin is a saddle. With $\Delta>0$, the sign of
$\tau$ distinguishes attraction and repulsion; $\tau^2-4\Delta<0$ gives a
focus. Real eigenvectors give invariant straight lines. On $\Delta=0$ or
$\tau=0,\Delta>0$, nonlinear terms may decide stability in a nonlinear system.

The last example has negative eigenvalues but can **grow transiently**:
nonorthogonal eigendirections allow components to reinforce one another
before eventual decay. Eigenvalues describe asymptotic behavior.
"""

# ╔═╡ 9215c291-887f-51ac-90df-08cb443dd576
function course_linear!(du,u,A,t)
    du[1]=A[1,1]*u[1]+A[1,2]*u[2]
    du[2]=A[2,1]*u[1]+A[2,2]*u[2]
    nothing
end

# ╔═╡ fd0b13aa-0cef-5a67-aae1-c4f9197f085c
course_matrices = Dict(
    "Stable node"=>[-0.5 0.0;0.0 -1.0],
    "Saddle"=>[0.5 0.0;0.0 -1.0],
    "Stable focus"=>[-0.1 -1.0;1.0 -0.1],
    "Center"=>[0.0 -1.0;1.0 0.0],
    "Transient growth"=>[-0.5 3.0;0.0 -1.0])

# ╔═╡ d8a2f0fd-e935-538d-bc07-3652a51e987a
course_linear_choice = "Stable focus"

# ╔═╡ 4dab07d0-14f0-59c3-860b-05730b8c6274
course_linear_u0 = [0.0,1.0]

# ╔═╡ 88bba503-9f28-59d4-9900-7e2a5f842a2f
course_linear_tspan = (0.0,12.0)

# ╔═╡ 6884da5e-0af1-5a10-84d8-a9bcd93a3645
course_linear_p = course_matrices[course_linear_choice]

# ╔═╡ 9408dd14-2e74-5be8-819d-9ac7617f5785
course_linear_prob = ODEProblem(course_linear!, course_linear_u0, course_linear_tspan, course_linear_p)

# ╔═╡ 54717269-122b-5f9d-a2fe-796511fd723c
course_linear_sol = solve(course_linear_prob, Tsit5(); abstol=1e-9, reltol=1e-8);

# ╔═╡ 00310fe9-7441-5c86-a815-bb228bac47ac
classification_linear(course_linear_p;Ngrid=3,tmax=2.0)

# ╔═╡ f59dec87-2b59-5f5a-bb6d-37b4c4a7ede5
(trace=tr(course_linear_p), determinant=det(course_linear_p),
 eigenvalues=eigvals(course_linear_p), eigenvectors=eigen(course_linear_p).vectors)

# ╔═╡ 83175c48-f08b-54aa-874e-cd6add8ac4b7
let a=plot(course_linear_sol;idxs=(1,2),xlabel="x",ylabel="y",legend=false),
    b=plot(course_linear_sol.t,norm.(course_linear_sol.u);xlabel="t",
        ylabel="‖u(t)‖",legend=false)
    plot(a,b;layout=(1,2),size=(900,370),margin=5*Plots.mm)
end

# ╔═╡ 25f8f2e4-9266-5033-b360-57ccfcdc8548
md"""
**Try:** classify each matrix from its trace and determinant before selecting
it. For a saddle, start exactly on its stable eigendirection and then perturb
that initial condition. For transient growth, explain why an increasing norm
over a finite interval does not contradict negative eigenvalues.
"""

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
# ╟─79a3579c-b488-5854-a2c6-d7c46efdabc9
# ╠═9215c291-887f-51ac-90df-08cb443dd576
# ╠═fd0b13aa-0cef-5a67-aae1-c4f9197f085c
# ╠═d8a2f0fd-e935-538d-bc07-3652a51e987a
# ╠═4dab07d0-14f0-59c3-860b-05730b8c6274
# ╠═88bba503-9f28-59d4-9900-7e2a5f842a2f
# ╠═6884da5e-0af1-5a10-84d8-a9bcd93a3645
# ╠═9408dd14-2e74-5be8-819d-9ac7617f5785
# ╠═54717269-122b-5f9d-a2fe-796511fd723c
# ╠═00310fe9-7441-5c86-a815-bb228bac47ac
# ╠═f59dec87-2b59-5f5a-bb6d-37b4c4a7ede5
# ╠═83175c48-f08b-54aa-874e-cd6add8ac4b7
# ╟─25f8f2e4-9266-5033-b360-57ccfcdc8548
