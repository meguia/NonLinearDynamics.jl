### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ 6cc0c64c-795c-4e97-aea7-9fb39ea0e2cd
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", "..", ".."))
    Pkg.instantiate()
    using PlutoUI, Plots, DifferentialEquations, ForwardDiff, BifurcationKit, Setfield
    using NonLinearDynamics
end

# ╔═╡ 0a2d892a-b321-4414-9528-a6bf38db0a60
gr();

# ╔═╡ 7a81862d-099a-4e92-9470-e54a338ceba0
md"""
Bifurcations correspond to qualitative changes in the dynamics, or topological structure, of dynamical systems. But when we refer to "changes" they are changes from the outside, by modifying the parameters of the system we modify the dynamical system.

We can adopt a more general point of view to understand how these changes are produced if we consider a dynamical system also as a function of a parameter (up to now all bifurcations we have seen were crossed by varying a single parameter). In that case we have to consider now the vector fields as a function of the variables and the bifurcation parameter ($\mu$) , i.e. we are going to write for 1D flows:

$\dot{x}=f(x,\mu)$

For example, for the case of the pitchfork normal form $f(x,mu)=mu x-x^3$. 

And let us also assume that the point $x=0$ is a fixed point and that the bifurcation occurs for a parameter value $\mu=0$ (this can always be achieved through a linear change of variables). In that case we have for the fixed point $f(0,0)=0$ the following conditions for bifurcations occurring on a one-dimensional manifold:

- the condition that an eigenvalue of the linearized flow becomes zero $\lambda=0$, which is equivalent to say that the derivative of the vector field with respect to the variable is zero , i.e. $\partial_x f(0,0)=0$. Note that we use the partial derivative because we now our vector field depends on the parameter explicitly. This is also known as the **non-hyperbolicity** condition. Recall that hyperbolic points were those with strictly positive or negative eigenvalues.
- There are two additional conditions. The first is called the **transversality** condition and requires that the function $f(x,mu)$ cuts the horizontal plane $f=0$ transversely, i.e. that its derivative with respect to the parameter does not also cancel out.
- the second is the **non-degeneracy** condition and requires that the second derivative of $f$ (with respect to x) be nonzero.
"""

# ╔═╡ 14dcf9c0-a184-467c-a48a-49ff86e492a4
md"""
If we take a fixed point of any dynamical system, with an arbitrary value of $\mu$, **the generic condition is that it is hyperbolic**. It is also said that in an environment of that fixed point the system is **robust** or **structurally stable**. If we perturb $\mu$ its dynamics does not change qualitatively. In that case the linearization determines the dynamics in its environment and there is nothing novel. 

More interesting things start to happen with as the genericity conditions are lost. The first loss of genericity is the condition of **non-hyperbolicity**. In that case a bifurcation arises. The condition is not generic because you have to tune to the parameter to hit the bifurcation. And at that point the system is not structurally stable anymore: if we disturb it the topological structure of its flow changes. 

If we just allow that loss of genericity we will always have a **saddle node** bifurcation. That is why it is made explicit that two other conditions of genericity must be met: transversality and non-degeneracy. 

But we can keep losing genericity and see what happens. If the position of the fixed point does not vary with the parameter the transversality condition is violated. And then we have the **transcritical** bifurcation (which is a bifurcation with a fixed point that does not move in parameter space).

If we now further impose a symmetry on our system, that violates the non-degeneracy condition because it cancels the second derivative of the vector field, and we have the **pitchfork** bifurcation.
                    
We can then summarize what we have so far for **local bifurcations** in terms of the conditions for them to occur and the corresponding normal form. 
"""

# ╔═╡ eb5d56d1-66ed-4dde-8d55-da2036fd2969
md"""
| Bifurcation 	| Condition | Non hyperbolicity   | Transversality 	| Non degeneracy | Normal Form      |
|---------------|:--------------:|:--------------:|:----------:|:----------:|:-------------:|
| Saddle Node 	| $\lambda=0$ |$\partial_x f(0,0)=0$ | $\partial_\mu f(0,0)\neq0$|$\partial_{xx} f(0,0)\neq0$ | $\dot{x}=\mu-x^2$|
| Transcritical 	| $\lambda=0$ |$\partial_x f(0,0)=0$ | $\partial_\mu f(0,0)=0$|$\partial_{xx} (0,0)\neq0$ | $\dot{x}=\mu x-x^2$|
| Pitchfork  	|  $\lambda=0$ | $\partial_x f(0,0)=0$ | $\partial_\mu f(0,0)=0$| $\partial_{xx} (0,0)=0$ | $\dot{x}=\mu x-x^3$|
| Hopf 	        | $Re(\lambda_{1,2})=0$ | $\lambda_{1,2}=\pm i\omega$| $\partial_\mu \beta(0)\neq0$|$l_1(0)\neq0$|$\dot{\rho}=\beta \rho-\rho^3$|
"""

# ╔═╡ Cell order:
# ╠═6cc0c64c-795c-4e97-aea7-9fb39ea0e2cd
# ╟─0a2d892a-b321-4414-9528-a6bf38db0a60
# ╟─7a81862d-099a-4e92-9470-e54a338ceba0
# ╟─14dcf9c0-a184-467c-a48a-49ff86e492a4
# ╟─eb5d56d1-66ed-4dde-8d55-da2036fd2969
