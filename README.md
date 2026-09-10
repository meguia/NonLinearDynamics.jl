# NonLinearDynamics

Pluto notebooks for learning nonlinear dynamics, following Steven Strogatz's
*Nonlinear Dynamics and Chaos*, with examples from musical instruments.

## Run the notebooks

1. Install [Julia 1.12](https://julialang.org/downloads/). The shared environment
   was generated with Julia 1.12.7 and uses Pluto 1.0.
2. Download this repository as a ZIP and extract it, or clone it. Keep the whole
   folder: every notebook uses the root `Project.toml`, `Manifest.toml`, and `src`.
3. Start Julia and run:

   ```julia
   import Pkg
   Pkg.add("Pluto")  # Only needed once.
   using Pluto
   Pluto.run()
   ```

4. Open a notebook from one of the four folders below. Students can start with
   `src/Pluto/01_Interactive/NLD01_Maps.jl`; the two coding tracks start at `NLD03`.

The first opening installs the shared dependencies automatically. An internet
connection is needed for this installation, and compilation can take several
minutes. Later openings reuse the environment. There is no need to install
individual numerical, plotting, or audio packages or change Julia's working
directory. Embedded remote images and videos still need internet access.

## Four tracks

| Folder | Approach |
| --- | --- |
| [01_Interactive](src/Pluto/01_Interactive) | Student lessons with explanations, sliders, and visual experiments. |
| [02_ModelingToolkit](src/Pluto/02_ModelingToolkit) | Visible variables → parameters → differential → equations → system → problem → solve → plot. |
| [03_ODE](src/Pluto/03_ODE) | Visible function → initial conditions → time span → parameters → problem → algorithm/solve → plot. |
| [04_RealTimeAudio](src/Pluto/04_RealTimeAudio) | The same ODEs connected to live audio, with parameter, pitch, gain, and Play controls. |

The MTK and ODE tracks contain 14 paired examples with matching equations,
parameters, initial conditions, and solver settings. Their numbers follow the
interactive sequence; maps and the introductory one-dimensional flows remain in
the interactive track. Later one-dimensional applications retain their lesson
numbers. The longer interactive lessons include additional models and analyses.

| Number | Topic | Interactive | MTK | ODE | Live audio |
| --- | --- | --- | --- | --- | --- |
| 01 | Maps | [Open](src/Pluto/01_Interactive/NLD01_Maps.jl) | — | — | — |
| 02 | One-dimensional flows | [Open](src/Pluto/01_Interactive/NLD02_Flows1D.jl) | — | — | — |
| 03 | Flows in 2D: damped oscillator | [Open](src/Pluto/01_Interactive/NLD03_Flows2D.jl) | [Open](src/Pluto/02_ModelingToolkit/NLD03_Flows2D.jl) | [Open](src/Pluto/03_ODE/NLD03_Flows2D.jl) | [Open](src/Pluto/04_RealTimeAudio/NLD03_Flows2D.jl) |
| 03b | Coupled phases on a torus | [Open](src/Pluto/01_Interactive/NLD03b_FlowsS1S2.jl) | [Open](src/Pluto/02_ModelingToolkit/NLD03b_FlowsS1S2.jl) | [Open](src/Pluto/03_ODE/NLD03b_FlowsS1S2.jl) | — |
| 04 | Bouncing ball: compliant contact | [Open](src/Pluto/01_Interactive/NLD04_BouncingBall.jl) | [Open](src/Pluto/02_ModelingToolkit/NLD04_BouncingBall.jl) | [Open](src/Pluto/03_ODE/NLD04_BouncingBall.jl) | — |
| 05 | One-dimensional bifurcations: pitchfork | [Open](src/Pluto/01_Interactive/NLD05_Flows1D_Bifurcations.jl) | [Open](src/Pluto/02_ModelingToolkit/NLD05_Flows1D_Bifurcations.jl) | [Open](src/Pluto/03_ODE/NLD05_Flows1D_Bifurcations.jl) | — |
| 06 | Population growth with harvesting | [Open](src/Pluto/01_Interactive/NLD06_Population_Models.jl) | [Open](src/Pluto/02_ModelingToolkit/NLD06_Population_Models.jl) | [Open](src/Pluto/03_ODE/NLD06_Population_Models.jl) | — |
| 07 | Love affairs: a linear system | [Open](src/Pluto/01_Interactive/NLD07_LoveAffairs.jl) | [Open](src/Pluto/02_ModelingToolkit/NLD07_LoveAffairs.jl) | [Open](src/Pluto/03_ODE/NLD07_LoveAffairs.jl) | — |
| 08 | Self-oscillation: Rayleigh reed | [Open](src/Pluto/01_Interactive/NLD08_Self_Oscillators.jl) | [Open](src/Pluto/02_ModelingToolkit/NLD08_Self_Oscillators.jl) | [Open](src/Pluto/03_ODE/NLD08_Self_Oscillators.jl) | — |
| 09 | The Duffing oscillator | [Open](src/Pluto/01_Interactive/NLD09_Flows2D_Duffing.jl) | [Open](src/Pluto/02_ModelingToolkit/NLD09_Flows2D_Duffing.jl) | [Open](src/Pluto/03_ODE/NLD09_Flows2D_Duffing.jl) | — |
| 10 | Periodically forced Duffing oscillator | [Open](src/Pluto/01_Interactive/NLD10_Flows2D_Forced.jl) | [Open](src/Pluto/02_ModelingToolkit/NLD10_Flows2D_Forced.jl) | [Open](src/Pluto/03_ODE/NLD10_Flows2D_Forced.jl) | [Open](src/Pluto/04_RealTimeAudio/NLD10_Flows2D_Forced.jl) |
| 11 | Flows in 3D: Lorenz system | [Open](src/Pluto/01_Interactive/NLD11_Flows3D.jl) | [Open](src/Pluto/02_ModelingToolkit/NLD11_Flows3D.jl) | [Open](src/Pluto/03_ODE/NLD11_Flows3D.jl) | — |
| 12 | Bifurcations: the Bogdanov–Takens model | [Open](src/Pluto/01_Interactive/NLD12_bifurcations.jl) | [Open](src/Pluto/02_ModelingToolkit/NLD12_bifurcations.jl) | [Open](src/Pluto/03_ODE/NLD12_bifurcations.jl) | — |
| 13 | Reed oscillator coupled to a bore resonance | [Open](src/Pluto/01_Interactive/NLD13_Reed_Resonator.jl) | [Open](src/Pluto/02_ModelingToolkit/NLD13_Reed_Resonator.jl) | [Open](src/Pluto/03_ODE/NLD13_Reed_Resonator.jl) | [Open](src/Pluto/04_RealTimeAudio/NLD13_Reed_Resonator.jl) |
| 14 | Bowed oscillator coupled to a body resonance | [Open](src/Pluto/01_Interactive/NLD14_Bowed_Resonator.jl) | [Open](src/Pluto/02_ModelingToolkit/NLD14_Bowed_Resonator.jl) | [Open](src/Pluto/03_ODE/NLD14_Bowed_Resonator.jl) | [Open](src/Pluto/04_RealTimeAudio/NLD14_Bowed_Resonator.jl) |
| 15 | Struck nonlinear oscillator and resonator | [Open](src/Pluto/01_Interactive/NLD15_Struck_Resonator.jl) | [Open](src/Pluto/02_ModelingToolkit/NLD15_Struck_Resonator.jl) | [Open](src/Pluto/03_ODE/NLD15_Struck_Resonator.jl) | [Open](src/Pluto/04_RealTimeAudio/NLD15_Struck_Resonator.jl) |

Interactive NLD12 also covers numerical continuation of equilibria, periodic
orbits, and homoclinic bifurcations. Some expensive Duffing cells remain disabled;
Pluto's cell menu can enable them. The stochastic extra has an optional slow scan.

The original explorations are preserved in `Extras`: stochastic dynamics in the
interactive track, double Hopf in MTK, and generalized Hopf, Takens, bell/clapper,
and pumped oscillator examples in ODE. Earlier drafts remain in
`01_Interactive/.temp`.

## Musical instruments and WAVs

NLD13–NLD15 appear in all four tracks. Each couples a nonlinear source to a damped
resonator through equal and opposite spring forces:

- **Reed and bore:** Rayleigh negative damping supplies energy at small velocity;
  nonlinear damping saturates the oscillation.
- **Bowed oscillator and body:** a smooth velocity-weakening friction law drives a
  string mode coupled to a body mode.
- **Struck oscillator and resonator:** an initial strike excites a hardening spring
  and a second, potentially inharmonic, resonance.

These dimensionless teaching models isolate excitation, coupling, and resonance.
They are not calibrated simulations of complete instruments. The resonator
coordinate `q` is the sound-output proxy, and `τ = 2π f₀ t` sets the frequency scale.

The offline notebooks show waveforms, source phase portraits, and resonator
amplitude envelopes. They include browser playback and a **WAV download button**.
Solutions are sampled at 4× the output sample rate, low-pass filtered and downsampled
to 48 kHz, stripped of DC, peak-normalized, and faded at the ends. Normalization
changes loudness; use the unnormalized plots to compare physical amplitudes.

Three generated, four-second mono PCM16 examples are included in [audio](audio):
[reed](audio/reed_resonator.wav), [bowed](audio/bowed_resonator.wav), and
[struck](audio/struck_resonator.wav). To regenerate them from the ODE notebooks:

```sh
julia --project=. scripts/render_audio.jl
```

Live audio starts only when **Play** is checked and requires an available local
audio output device. Uncheck it to stop; for a new strike, stop and play again.
Parameter sliders update the running source. Live output uses the selected gain;
the offline filtering and normalization apply to WAV exports. PlutoHooks stops
the source when its playback cell is rerun or removed, or its notebook is shut down.

## Shared environment

Every notebook in the four tracks follows Pluto's
[shared-environment pattern](https://plutojl.org/en/docs/packages-advanced/):
its first code cell activates the repository relative to the notebook file,
instantiates the manifest, and imports its packages. Download the whole repository
instead of opening a single raw notebook URL. The only local package path is the
repository itself (`.`); notebooks use no separate embedded environments.

The bounds retain the APIs used by the lessons: ModelingToolkit 9,
IntervalRootFinding 0.5, BifurcationKit 0.3.6 with HclinicBifurcationKit 0.1.1,
and compatible KrylovKit, DiffEqBase, SciMLBase, and DynamicQuantities versions.
[RealTimeAudioDiffEq](https://github.com/antonioortegabrook/RealTimeAudioDiffEq.jl)
0.1 shares this DifferentialEquations 7 stack; its newer 0.2 series requires a
newer SciMLBase. DSP 0.7.10 shares the existing Polynomials 3 dependency.
Third-party dependencies are registered and their resolved versions are fixed by the
root manifest.

To verify installation outside Pluto, run from the repository folder:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
```

Update `Project.toml` and `Manifest.toml` together. After changes, check enabled
cells in Pluto, compare the paired MTK/ODE trajectories, and regenerate the audio.

Validation on Julia 1.12.7: all 57 notebooks pass Pluto's structural checks; all
36 new notebooks and five preserved extras execute without cell errors. The
local numerical/audio regression suite passes 1,730 checks, including agreement
between the paired tracks. Realtime callbacks were tested in memory; speaker
playback and the intentionally disabled legacy calculations were not exercised.
