# NonLinearDynamics

Pluto notebooks for a nonlinear dynamics course, following Steven Strogatz's
*Nonlinear Dynamics and Chaos*, with examples from musical instruments.

**Interactive** has sliders and visual experiments. **ODE** shows the equations,
function, initial conditions, parameters, problem, solver, and plots.
**RealTime** contains selected examples for live sound.

Install Julia 1.12 and download this repository. Keep the whole folder, then
start Pluto from Julia:

```julia
import Pkg
Pkg.add("Pluto")  # Once.
using Pluto
Pluto.run()
```

Open a notebook from the list below. Its first cell activates the shared
Project and Manifest and installs the required packages on the first run.
The musical examples in Interactive and ODE generate audio with WAV download
buttons. In RealTime, use Play to start and stop the sound.

| Notebook | Interactive | ODE | RealTime |
| --- | --- | --- | --- |
| 01 · Maps and iteration | [Open](src/Pluto/Interactive/NLD01_Maps.jl) | — | — |
| 02 · Flows in one dimension | [Open](src/Pluto/Interactive/NLD02_Flows1D.jl) | — | — |
| 03 · Flows in two dimensions | [Open](src/Pluto/Interactive/NLD03_Flows2D.jl) | [Open](src/Pluto/ODE/NLD03_Flows2D.jl) | — |
| 03b · Circle dynamics and coupled phases | [Open](src/Pluto/Interactive/NLD03b_FlowsS1S2.jl) | [Open](src/Pluto/ODE/NLD03b_FlowsS1S2.jl) | — |
| 04 · Bouncing ball | [Open](src/Pluto/Interactive/NLD04_BouncingBall.jl) | [Open](src/Pluto/ODE/NLD04_BouncingBall.jl) | — |
| 05 · One-dimensional bifurcations | [Open](src/Pluto/Interactive/NLD05_Flows1D_Bifurcations.jl) | [Open](src/Pluto/ODE/NLD05_Flows1D_Bifurcations.jl) | — |
| 06 · Population models | [Open](src/Pluto/Interactive/NLD06_Population_Models.jl) | [Open](src/Pluto/ODE/NLD06_Population_Models.jl) | — |
| 07 · Linear systems and love affairs | [Open](src/Pluto/Interactive/NLD07_LoveAffairs.jl) | [Open](src/Pluto/ODE/NLD07_LoveAffairs.jl) | — |
| 08 · Self-oscillators | [Open](src/Pluto/Interactive/NLD08_Self_Oscillators.jl) | [Open](src/Pluto/ODE/NLD08_Self_Oscillators.jl) | [Open](src/Pluto/RealTime/NLD08_Self_Oscillators.jl) |
| 09 · Duffing and other planar flows | [Open](src/Pluto/Interactive/NLD09_Flows2D_Duffing.jl) | [Open](src/Pluto/ODE/NLD09_Flows2D_Duffing.jl) | [Open](src/Pluto/RealTime/NLD09_Flows2D_Duffing.jl) |
| 10 · Forced oscillators | [Open](src/Pluto/Interactive/NLD10_Flows2D_Forced.jl) | [Open](src/Pluto/ODE/NLD10_Flows2D_Forced.jl) | [Open](src/Pluto/RealTime/NLD10_Flows2D_Forced.jl) |
| 11 · Lorenz and chaos | [Open](src/Pluto/Interactive/NLD11_Flows3D.jl) | [Open](src/Pluto/ODE/NLD11_Flows3D.jl) | [Open](src/Pluto/RealTime/NLD11_Flows3D.jl) |
| 12 · Local, global, and codimension-two bifurcations | [Open](src/Pluto/Interactive/NLD12_bifurcations.jl) | [Open](src/Pluto/ODE/NLD12_bifurcations.jl) | — |
| 13 · Reed, bowed, and struck resonators | [Open](src/Pluto/Interactive/NLD13_Source_Filter.jl) | [Open](src/Pluto/ODE/NLD13_Source_Filter.jl) | [Open](src/Pluto/RealTime/NLD13_Source_Filter.jl) |
| 14 · Two resonators and quasiperiodicity | [Open](src/Pluto/Interactive/NLD14_Two_Resonators.jl) | [Open](src/Pluto/ODE/NLD14_Two_Resonators.jl) | [Open](src/Pluto/RealTime/NLD14_Two_Resonators.jl) |
| 15 · Vocal folds and birdsong | [Open](src/Pluto/Interactive/NLD15_Examples.jl) | [Open](src/Pluto/ODE/NLD15_Examples.jl) | [Open](src/Pluto/RealTime/NLD15_Examples.jl) |

Additional notebooks:

- [Stochastic dynamics](src/Pluto/Interactive/Extras/NLD_Stochastic.jl)
- [Generalized Hopf](src/Pluto/ODE/Extras/NLD_Extras_GHopf.jl)
- [Takens unfolding](src/Pluto/ODE/Extras/NLD_Extras_Takens.jl)
- [Bell and clapper](src/Pluto/ODE/Extras/Others_Bell.jl)
- [Pumped oscillator](src/Pluto/ODE/Extras/Others_Pumped.jl)

Earlier notes remain in Interactive/.temp:
[codimension-two bifurcations](src/Pluto/Interactive/.temp/NLD10_Flows2D_Codimension2.jl)
and [Julia sets](src/Pluto/Interactive/.temp/NLD11_Fractals.jl).
