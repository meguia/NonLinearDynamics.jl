### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ fce7054d-cb51-45af-a3d1-e370fd3d9850
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", "..", ".."))
    Pkg.instantiate()
    using Plots, DifferentialEquations, Setfield, ForwardDiff, PlutoUI, IntervalRootFinding, StaticArrays
    import BifurcationKit as BK
    using NonLinearDynamics
end

# ╔═╡ ccd5d840-06c3-4138-999a-24c5ff7cb7aa
md"""
# Bogdanov-Takens bifurcation (codimension 2)

The bifurcations (of codimension 1) we have seen so far can be characterized in two groups:

- Those that happen when a real eigenvalue becomes zero. Generically we have a saddle-node bifurcation, but also a pitchfork or a transcritical one if other symmetry conditions are given.
- The one that happens when the real part of two conjugate complex eigenvalues becomes zero. In that case we have a Hopf bifurcation.

Clearly the first case can happen in a 1D system (or one of higher dimension along a particular direction), while the second we need at least a 2D system to have two eigenvalues, but it takes only one parameter to control it (to move the real part). 

But if we are in a 2D system, couldn't it happen that **both** eigenvalues become zero simultaneously? Clearly if we look at the eigenvalue expression we will need at least two parameters to adjust this point. With the usual nondegeneracy conditions and one eigenvector, this double zero gives a Bogdanov–Takens bifurcation. Its unfolding has nearby saddle-node, Hopf, and homoclinic bifurcation curves; it is distinct from a cusp.

This codimension-two bifurcation is known as Bogdanov–Takens (or Takens–Bogdanov). The normal form is characterized by having the following Jacobian:

$\begin{pmatrix}0 & 1\\0 & 0\end{pmatrix}$

which "induces" the following nonlinear terms to appear (in Bogdanov's version):

$\dot{x} = y$

$\dot{y} = x^2-xy$

Note that we did not introduce any parameters yet, this is the "pure" singularity. To extend this in parameter space (or dynamical systems to be more precise) it is necessary to do an **unfolding**, and here there are several possibilities, let's take the one done by Guckenheimer & Holmes:

$\dot{x} = y$

$\dot{y} = \mu_1+\mu_2x+ x^2 -xy$

Let's study the bifurcations directly without worrying about the solutions yet because this system (like the saddle node in the plane) has divergent trajectories.

Since we have terms up to quadratic order we will be able to have generically two fixed points or none. Always located on the horizontal axis $y=0$ and with the $x$ coordinate at:

$x_{\pm}=-\frac{\mu_2}{2}\pm \sqrt{\frac{\mu_2^2}{4}-\mu_1}$

The positive sign corresponds to the fixed point on the right and the negative sign to the fixed point on the left (when they exist).

the condition for existence of the fixed points is that the interior of the root is positive which gives us a condition to draw a saddle-node bifurcation curve in the plane $(\mu_1,\mu_2)$:

SN : $\mu_1=\frac{\mu_2^2}{4}$

On the other hand the Jacobian evaluated at the fixed points (preserving the order $pm$) is written as:

$\begin{pmatrix}0 & 1\\
\pm 2\sqrt{\frac{\mu_2^2}{4}-\mu_1} & \frac{\mu_2}{2}\mp \sqrt{\frac{\mu_2^2}{4}-\mu_1}
\end{pmatrix}$

Which gives us the determinant:

$\Delta = \mp 2\sqrt{\frac{\mu_2^2}{4}-\mu_1}$

which for the fixed point on the right is always negative (saddle point) and for the one on the left positive. For the latter we evaluate the trace (we only keep the sign below):

$\tau = \frac{\mu_2}{2}+\sqrt{\frac{\mu_2^2}{4}-\mu_1}$

since the root is positive when $\mu_2<0$, the trace will become zero when $\mu_1=0$. This gives us the condition to trace the Hopf bifurcation curve in the plane $(\mu_1,\mu_2)$:

Hopf: $\mu_2<0$  , $\mu_1=0$

This curve meets the SN parabola at the point $(0,0)$ so as we anticipated at this singular point we have an SN curve and a Hopf curve occurring simultaneously. The complete bifurcation diagram (taken from Scholarpedia) is:


Notice that there is an additional curve (in red) that corresponds to a global bifurcation (homoclinic connection). 

We have already seen this bifurcation! And with this same normal form, but arbitrarily fixing $\mu_2=-1$. This is a homoclinic bifurcation or saddle loop. The change from region (3) to (4) is the one we described when we saw this bifurcation. The limit cycle originating from the Hopf curve grows until it touches the saddle and saddle varieties and forms a loop. On the other side of the bifurcation there is no limit cycle and the unstable saddle manifold feeds an attractor focus. 

Explore how the moieties are modified in the graph below and try to locate when the homoclinic connection occurs. As a guide we show on the right the bifurcation diagram with the analytic SN and Hopf curves and in dotted line the homoclinic that occurs (we will not show the deduction of that) when:

HC: $\mu_1 \approx -\frac{6}{25}\mu_2^2$ ,  $\mu_2<0$
"""

# ╔═╡ 7dfefdc3-4d5c-4c5d-9855-66dbb2efc26c
md"""
## Bogdanov–Takens unfolding with cubic terms

As can be seen the BT bifurcation presents a very interesting and varied dynamics, with minimal alterations in the parameters we can go from oscillatory behaviors, creation of pairs of fixed points and infinite period orbits (HC connections). The problem with the above system is that it has diverging trajectories, so we need to add higher order terms (which will not alter the BT bifurcation but may change the bifurcation diagram outside that point) to ensure that the trajectories do not diverge. Back there are several alternatives, let's follow the one proposed by Mindlin:

$\dot{x} = y$

$\dot{y} = \mu_1+\mu_2x+ x^2 -xy - x^3 -x^2y$ 

In this case by having cubic terms we will have in general one or three fixed points, as in the case of the cusp the fixed points go from 1 to 3 through saddle-node bifurcations that occur in pairs of distinct points. In fact the cubic terms introduce a cusp in addition to the Bogdanov-Takens.

Let us see how the varieties are organized from these new terms.

"""

# ╔═╡ 335db94d-98ba-4dc5-ae09-0d4544c00ade
function takens3!(du,u,p,t)
    du[1]=u[2]
    du[2]=p[1]+u[1]*(p[2]-u[2]+u[1]*(1-u[1]-u[2]))
    du
end    

# ╔═╡ 8e470606-784f-4b02-9a43-2a351b7b8dbe
md"""
μ1 $(@bind μ1 Slider(-0.12:0.001:0.02,default=-0.01;show_value=true)) 
μ2 $(@bind μ2 Slider(-0.3:0.001:0.0,default=-0.15;show_value=true)) \
tmax $(@bind tmax Slider(100:10:400,default=100;show_value=true)) \
"""

# ╔═╡ 5852171c-8d21-454e-8121-c8fe9a9a9409
phase_portrait(takens3!,[μ1,μ2];tmax=tmax,xlims=[-1,1],ylims=[-0.5,0.5])

# ╔═╡ Cell order:
# ╠═fce7054d-cb51-45af-a3d1-e370fd3d9850
# ╟─ccd5d840-06c3-4138-999a-24c5ff7cb7aa
# ╟─7dfefdc3-4d5c-4c5d-9855-66dbb2efc26c
# ╠═335db94d-98ba-4dc5-ae09-0d4544c00ade
# ╟─8e470606-784f-4b02-9a43-2a351b7b8dbe
# ╠═5852171c-8d21-454e-8121-c8fe9a9a9409
