# ODE

Define an in-place function, initial conditions, time span, and parameters; build an ODEProblem, choose a solver, and plot the solution.

Open these files with Pluto; their first cell activates and installs the shared
repository environment automatically. See the [course index](../../../README.md)
for the matching topics in the other tracks.

`model!(du, u, p, t)` writes each derivative into `du`; `sol[i, :]` selects the
whole time series of state `i`, and `plot(sol; idxs=(1,2))` draws a phase portrait.

- [NLD03_Flows2D](NLD03_Flows2D.jl)
- [NLD03b_FlowsS1S2](NLD03b_FlowsS1S2.jl)
- [NLD04_BouncingBall](NLD04_BouncingBall.jl)
- [NLD05_Flows1D_Bifurcations](NLD05_Flows1D_Bifurcations.jl)
- [NLD06_Population_Models](NLD06_Population_Models.jl)
- [NLD07_LoveAffairs](NLD07_LoveAffairs.jl)
- [NLD08_Self_Oscillators](NLD08_Self_Oscillators.jl)
- [NLD09_Flows2D_Duffing](NLD09_Flows2D_Duffing.jl)
- [NLD10_Flows2D_Forced](NLD10_Flows2D_Forced.jl)
- [NLD11_Flows3D](NLD11_Flows3D.jl)
- [NLD12_bifurcations](NLD12_bifurcations.jl)
- [NLD13_Reed_Resonator](NLD13_Reed_Resonator.jl)
- [NLD14_Bowed_Resonator](NLD14_Bowed_Resonator.jl)
- [NLD15_Struck_Resonator](NLD15_Struck_Resonator.jl)
- [NLD16_Two_Resonators](NLD16_Two_Resonators.jl)

Original exploratory material is preserved in [Extras](Extras).
