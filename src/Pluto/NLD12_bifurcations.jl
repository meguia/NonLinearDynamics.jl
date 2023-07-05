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

# ╔═╡ 87428079-6c7e-409a-988a-83a7bd59ae2c
import Pkg; Pkg.add(url="https://github.com/bifurcationkit/HclinicBifurcationKit.jl"); Pkg.add("Setfield")

# ╔═╡ 5e1d221c-18f7-11ee-20bc-b5624af0581e
using Plots, Interact, DifferentialEquations, Setfield, ForwardDiff, PlutoUI, IntervalRootFinding, StaticArrays

# ╔═╡ f1702b1b-4902-4e43-adb9-402adab2a4b5
include("../NLD_utils.jl")

# ╔═╡ 1e046183-5180-4e65-94ed-2ef6ee869c3c
import BifurcationKit as BK

# ╔═╡ 244b8a5e-ae4f-4965-9dcf-df860e33015e
import HclinicBifurcationKit as HBK

# ╔═╡ 77fd29b5-4cc6-4461-ab16-3f0720643dee
md"""
## Bogdanov Takens with cubic terms

This is the unforlding of the BT normal form with cubic terms added to prevent the divergence of the orbits in the form proposed by Mindlin:

$\dot{x} = y$

$\dot{y} = \mu_1+\mu_2x+ x^2 -xy - x^3 -x^2y$ 

In this case by having cubic terms we will have in general one or three fixed points, as in the case of the cusp the fixed points go from 1 to 3 through chair node bifurcations that occur in pairs of distinct points. In fact the cubic terms introduce a cusp in addition to the Bogdanov-Takens.

Let us see how the varieties are organized from these new terms.

"""

# ╔═╡ 4380968c-9421-4c34-82bd-850c3027b676
function takens3!(du,u,p,t)
    du[1]=u[2]
    du[2]=p[1]+u[1]*(p[2]-u[2]+u[1]*(1-u[1]-u[2]))
    du
end    

# ╔═╡ e65f74c0-edf0-4378-ab3d-5e848bda5d2f
md"""
# Bifurcation Analysis

## Single parameter Hopf
"""

# ╔═╡ b4574225-ff15-44f5-82e2-5eb779a16c49
begin
	p = (μ1=0.1,μ2=-0.1)
	μ1 = @lens _.μ1
	μ2 = @lens _.μ2
	takens3(u,p) = takens3!(similar(u),u,p,0) # out of place method
	u0 = [0.9,-0.1] #initial condition
end	

# ╔═╡ 6ed8af70-4918-4dc2-b419-d87c8f35cfb9
md"""
μ2 $(@bind mu_2 Slider(-0.4:0.005:0.05,default=-0.15;show_value=true)) 
"""

# ╔═╡ 90461200-bc32-4f01-ac06-a7e9b81966df
begin
	# define the bifurcation problem
	prob1 = BK.BifurcationProblem(takens3, u0, set(p,μ2,mu_2), μ1, recordFromSolution = (x, p) -> (x = x[1], y = x[2]))
	# continuation options
	opts_br = BK.ContinuationPar(pMin=-0.2,pMax=0.1, ds = -0.001, dsmax = 0.02, detectBifurcation=3, nInversion=8)
	# continuation of equilibria
	br1 = BK.continuation(prob1, BK.PALC(tangent=BK.Bordered()),opts_br)
	scene = plot(br1,xlabel="\\mu_1",title=string("BT Cubica \\mu_2 = ",mu_2));
end

# ╔═╡ 9eff4af0-8cef-4a9a-9e6c-092a1e046b96
md"""
1. Definition of the Bifurcation Problem, giving the out of place function, initial condition, parameters (wher we fixed μ2 value with set), and then the parameter to be explored (μ1). With `recordFromSolution = (x, p) -> (x = x[1], y = x[2])` we are recording only the variables x, y

2. Options for the continuation algorithm: starting from the higher value of μ1 (0.1) with a stable fixed point and going backwards (ds<0) down to μ1=-0.2. `Detectbifurcation =3` improves accuracy and `nInversions=8` the convergence.

3. Continuation algorithm using Pseudo Arc-Length Continuation with a Bordered predictor. It is based in computing the tangent to the curve of solutions in an extended plane (x,p) by solving a bordered linear system $F_x dx + F_p dp = 0$ ; $\theta dx_0 dx + (1-\theta) dp_0 dp = 1$
"""

# ╔═╡ fa557b12-9017-4166-ac5b-ead0c3e108b0
for n = 1:(length(br1.specialpoint)-1)
	print(BK.getNormalForm(br1, n))
	println()
end	

# ╔═╡ 1b94418e-60b7-4da4-8720-4ac37155c909
md"""
The branch obtained contain information about the bifurcation points `br.specialpoints` and `getNormalForm` can compute the Normal Form at that points. 
"""

# ╔═╡ 4b0d25fc-551b-48e6-b4db-154c12b8391f
md"""
## Continuation of the Fold (Saddle Node) > CUSP
"""

# ╔═╡ ffb2f68b-493e-4ae9-8626-3eb5a49b0e17
begin
	# we generate a new branch for a fixed value of μ2
	prob1b = BK.BifurcationProblem(takens3, u0, set(p,μ2,-0.1), μ1, recordFromSolution = (x, p) -> (x = x[1], y = x[2]))
	br1b = BK.continuation(prob1b, BK.PALC(tangent=BK.Bordered()),opts_br)
	opts_br2 = BK.ContinuationPar(pMin=-0.4,pMax=0.3, ds = -0.001, dsmax = 0.02,nInversion=6)
	# continuation of critical point 
	br2 = @time BK.continuation(br1b, 1, μ2, opts_br2,detectCodim2Bifurcation=2,bdlinsolver=BK.MatrixBLS(),updateMinAugEveryStep=1,startWithEigen=true)
	scene2 = plot(br2);
end

# ╔═╡ ccf38f3c-d73e-4174-84bd-1399d64982f4
md"""
We first generate a new branch (br1b) for a fixed value of μ2=-0.1 as a starting point (-0.0841231,-0.1) in the parameter space (μ1,μ2), and then give that branch as an initial guess for a Fold point to continuation in order to calculate a curve of Fold points in parameter space based on a Minimally Augmented formulation. We have to indicate index of the critical point (1) as a second argument, the second parameter to be varied (μ2) and the options. 
In this case the options also include the minimum and maximum value for the parameter to be varied and the initial step ds (note that in this case we are also going backwards at first).
The Matrix based Bordered Linear Solver is used as the bordered linear solver for the constraint equation (recommended for ODEs)
`detectCodim2Bifurcation=2` serves to detect Cusp/Bogdanov-Takens/Bautin codimension 2 critical points precisely.
And finally `updateMinAugEveryStep=1` update vectors a, b in Minimally Formulation at each step.
"""

# ╔═╡ 29b9b64a-19fa-470f-9596-45f67b30b983
BK.getNormalForm(br2, 1)

# ╔═╡ 41f4608f-fad2-4202-a9be-c0a3dc5a379b
md"""
## Continuation of the Hopf > BOGDANOV TAKENS
"""

# ╔═╡ b5e04981-da64-4e1d-b7d5-e9a34028b26b
begin
	# we generate a new branch for another fixed value of μ2
	prob1c = BK.BifurcationProblem(takens3, u0, set(p,μ2,-0.3), μ1, recordFromSolution = (x, p) -> (x = x[1], y = x[2]))
	br1c = BK.continuation(prob1c, BK.PALC(tangent=BK.Bordered()),opts_br)
	opts_br3 = BK.ContinuationPar(pMin=-0.3,pMax=0.01, ds = 0.001, dsmax = 0.02,nInversion=6)
	# continuation of critical point 3 (Hopf)
	br3 = @time BK.continuation(br1c, 3, μ2, opts_br3,detectCodim2Bifurcation=2,updateMinAugEveryStep = 1)
end

# ╔═╡ 7083a569-b724-4281-a7d9-0074a357f541
md"""
This case is similar to the previous one with the exception that now we start with a initial guess of the Hopf for μ2=-0.3 and go forward.
"""

# ╔═╡ b143e725-7d5f-45e1-ad32-523766b579bc
begin
	plot(br2,branchlabel="SN")
	plot!(br3,branchlabel="Hopf")
end	

# ╔═╡ 311c4ba1-be7f-465a-b38d-bbb635305ca1
br3.specialpoint[1]

# ╔═╡ be1e8c2b-cd43-493d-9e59-b7a591f37749
btpt = BK.getNormalForm(br3, 1; nev = 3, autodiff = false)

# ╔═╡ 79d047e0-0660-4b27-b1c7-0c24f46ed1d4
md"""
## Crossing the Hopf for a fixed μ2

We will fix a value of μ2 and start the continuation from a small positive value of μ1 before the Hopf (μ1=0.01) and the go down computing the period. As we will aproach an Homoclinic Orbit the period can diverge.

We first compute the branch for this trajectory in parameter space as before
"""

# ╔═╡ c6175b02-6711-45c8-b499-8ef8a0a2b82e
begin
	p2 = (μ1=0.01, μ2=-0.16)
	prob1d = BK.BifurcationProblem(takens3, u0, p2, μ1, recordFromSolution = (x, p) -> (x = x[1], y = x[2]))
	opts_br1d = BK.ContinuationPar(pMin=-0.2,pMax=p2[:μ1], ds = -0.001, dsmax = 0.02, detectBifurcation=3, nInversion=8)
	br1d = BK.continuation(prob1d, BK.PALC(tangent=BK.Bordered()),opts_br1d)
end	

# ╔═╡ 548518de-0251-4bee-be0c-06a96d952e43
md"""
## Continuation of Periodic Orbits using the Trapezoid Method

The options are the same than those used for the branch `br1d`, adding tolerance and maximum step options

"""

# ╔═╡ 77cdc1b9-d115-4bc5-aa00-5b28b3bdfae5
opts_po = setproperties(opts_br1d, maxSteps = 150, tolStability = 1e-8);

# ╔═╡ 21692850-ffa6-48ba-99a8-470460ae89bc
function recordPO(x, p)
	xtt = BK.getPeriodicOrbit(p.prob, x, p.p)
	period = BK.getPeriod(p.prob, x, p.p)
	return (max = maximum(xtt[1,:]), min = minimum(xtt[1,:]), period = period)
end

# ╔═╡ 6152a0c0-4900-473c-8528-54d04b08eb0f
md"""
Then we call continuation with the computed branch and the index corresponding to the Hopf critical point and a functional that is encoded in the composite type `PeriodicOrbitTrapProblem`
"""

# ╔═╡ e3039d9c-b06b-47ba-9262-014b7eba5b2a
br_pot = @time BK.continuation(br1d, 3, opts_po, BK.PeriodicOrbitTrapProblem(M = 150));

# ╔═╡ 8ce52366-7bea-4019-89d3-dec24648dde9
md"""
## Periodic orbits with Parallel Standard Shooting

For this method we need to provide the ODE problem because the PSS is based on finding a solution to the flow $\Phi^T(x)=x$ for some period T, where $\Phi^t(x)$ is the flow of the system at time t. 
"""

# ╔═╡ 5af26e47-3a0f-4c52-9a9b-803093652531
odeprob = ODEProblem(takens3!, copy(u0), (0., 1.), p2; abstol = 1e-10, reltol = 1e-9);

# ╔═╡ f9681d1d-ac9a-4186-b654-77c991bcb5fa
md"""
This is similar to the previous one with the exception of the functional and some additional parameters to ensure convergence.
"""

# ╔═╡ 660fba6e-f7e0-4f42-a572-3507afd24a7e
br_pos = @time BK.continuation(br1d, 3, opts_po, BK.ShootingProblem(35, odeprob, Rodas5(), parallel = true); δp = 0.0001);

# ╔═╡ f96a46fc-e2a0-413f-a04c-e083ae9b5d3c
begin
	plot(br_pot, label="Trapezoidal")
	plot!(br_pos, label="Standard Shooting")
end

# ╔═╡ 96e1eac9-5956-48b5-a76c-cfd5c7644f12
md"""
The period of the orbit increases as we approach the Homoclinic Bifurcation 
"""

# ╔═╡ 699bb082-68c0-446f-b11e-7281861b7f13
md"""
## Branch of homoclinic orbits with Orthogonal Collocation
"""

# ╔═╡ dce93919-7291-46b4-83cb-10ac297a502c
odeprob_hc = ODEProblem(takens3!, copy(u0), (0., 1.), set(p,μ2,0.0); abstol = 1e-10, reltol = 1e-9);

# ╔═╡ c053ce00-e4bb-4f41-8323-b9ce9bc61e3f
prob_hc = BK.BifurcationProblem(takens3, u0, set(p,μ2,0.0), μ1, recordFromSolution = (x, p) -> (x = x[1], y = x[2]))

# ╔═╡ 1c2a4ece-df85-4dbb-8cc3-4a293d724b61
opt_new = BK.NewtonPar(verbose = true, tol = 1e-9, maxIter = 9)

# ╔═╡ 48829aee-b1ad-4382-b0be-bde0a9259241
opt_hc = BK.ContinuationPar(newtonOptions = opt_new, maxSteps = 100, saveSolEveryStep = 1, dsmax = 0.01, pMin = -0.1, ds = 0.001, dsmin = 1e-5, 
	detectEvent = 2, detectBifurcation = 0);

# ╔═╡ 3599952e-d26b-4a66-856b-85ea05759a1f
br_hom_c = BK.continuation(prob_hc, btpt, BK.PeriodicOrbitOCollProblem(50, 3; meshadapt = true, K = 100), BK.PALC(tangent = BK.Bordered()), opt_hc, 
		 verbosity=2, ϵ0 = 1e-5, amplitude = 2e-3,  freeparams = ((@lens _.T), (@lens _.ϵ0)), updateEveryStep = 2)

# ╔═╡ 61088c17-f8fd-4404-80e2-c9e911d5a858
plot!(br_hom_c)

# ╔═╡ 4bb8c192-7dc0-4780-acfb-34f189a3df4a
# Esto es para ensancha la caja por defecto
html"""
<style>
	main {
		margin: 0 auto;
		max-width: 1800px;
    	padding-left: max(160px, 10%);
    	padding-right: max(160px, 10%);
	}
</style>
"""

# ╔═╡ Cell order:
# ╠═87428079-6c7e-409a-988a-83a7bd59ae2c
# ╠═5e1d221c-18f7-11ee-20bc-b5624af0581e
# ╠═1e046183-5180-4e65-94ed-2ef6ee869c3c
# ╠═244b8a5e-ae4f-4965-9dcf-df860e33015e
# ╠═f1702b1b-4902-4e43-adb9-402adab2a4b5
# ╟─77fd29b5-4cc6-4461-ab16-3f0720643dee
# ╠═4380968c-9421-4c34-82bd-850c3027b676
# ╟─e65f74c0-edf0-4378-ab3d-5e848bda5d2f
# ╠═b4574225-ff15-44f5-82e2-5eb779a16c49
# ╟─6ed8af70-4918-4dc2-b419-d87c8f35cfb9
# ╠═90461200-bc32-4f01-ac06-a7e9b81966df
# ╟─9eff4af0-8cef-4a9a-9e6c-092a1e046b96
# ╠═fa557b12-9017-4166-ac5b-ead0c3e108b0
# ╟─1b94418e-60b7-4da4-8720-4ac37155c909
# ╟─4b0d25fc-551b-48e6-b4db-154c12b8391f
# ╠═ffb2f68b-493e-4ae9-8626-3eb5a49b0e17
# ╟─ccf38f3c-d73e-4174-84bd-1399d64982f4
# ╠═29b9b64a-19fa-470f-9596-45f67b30b983
# ╟─41f4608f-fad2-4202-a9be-c0a3dc5a379b
# ╠═b5e04981-da64-4e1d-b7d5-e9a34028b26b
# ╟─7083a569-b724-4281-a7d9-0074a357f541
# ╠═b143e725-7d5f-45e1-ad32-523766b579bc
# ╠═311c4ba1-be7f-465a-b38d-bbb635305ca1
# ╠═be1e8c2b-cd43-493d-9e59-b7a591f37749
# ╟─79d047e0-0660-4b27-b1c7-0c24f46ed1d4
# ╠═c6175b02-6711-45c8-b499-8ef8a0a2b82e
# ╟─548518de-0251-4bee-be0c-06a96d952e43
# ╠═77cdc1b9-d115-4bc5-aa00-5b28b3bdfae5
# ╠═21692850-ffa6-48ba-99a8-470460ae89bc
# ╟─6152a0c0-4900-473c-8528-54d04b08eb0f
# ╠═e3039d9c-b06b-47ba-9262-014b7eba5b2a
# ╟─8ce52366-7bea-4019-89d3-dec24648dde9
# ╠═5af26e47-3a0f-4c52-9a9b-803093652531
# ╟─f9681d1d-ac9a-4186-b654-77c991bcb5fa
# ╠═660fba6e-f7e0-4f42-a572-3507afd24a7e
# ╠═f96a46fc-e2a0-413f-a04c-e083ae9b5d3c
# ╟─96e1eac9-5956-48b5-a76c-cfd5c7644f12
# ╟─699bb082-68c0-446f-b11e-7281861b7f13
# ╠═dce93919-7291-46b4-83cb-10ac297a502c
# ╠═c053ce00-e4bb-4f41-8323-b9ce9bc61e3f
# ╠═1c2a4ece-df85-4dbb-8cc3-4a293d724b61
# ╠═48829aee-b1ad-4382-b0be-bde0a9259241
# ╠═3599952e-d26b-4a66-856b-85ea05759a1f
# ╠═61088c17-f8fd-4404-80e2-c9e911d5a858
# ╟─4bb8c192-7dc0-4780-acfb-34f189a3df4a
