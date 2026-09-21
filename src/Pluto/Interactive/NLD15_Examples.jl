### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 0504907f-de44-547a-b346-fe28fbd333df
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using NonLinearDynamics: prepare_audio
    using DifferentialEquations, LinearAlgebra, Plots, PlutoUI, WAV, Base64
end

# ╔═╡ 5f8331a5-e154-5b74-a325-69b4f0355937
md"""
# Examples: vocal folds and birdsong

Two models connect nonlinear oscillation to vocal sound. Begin with phonation
onset and saturation in a reduced vocal-fold oscillator, then follow a
syringeal source through a delayed vocal tract and an acoustic filter.
"""

# ╔═╡ 460bdbd3-5ddf-5936-ac14-00b94f57c234
TableOfContents()

# ╔═╡ 063a0f42-b8ad-5bcf-804c-d8999e17445b
md"""
## Vocal folds: a self-sustained oscillator

This reduced vocal-fold model represents one effective vibrating tissue coordinate
$x$ and its velocity $v$. The net drive $\mu$ supplies energy at small
displacements; the nonlinear damping $x^2v$ limits the oscillation. The parameter
$k$ controls stiffness, while $F$ is a constant bias force that shifts the
equilibrium. This is a model of phonation onset and saturation, without a vocal
tract or a collision law for fold contact.

Adapted from `DNL_Vocal_Folds.jl` in `NonlinearDynamicsPluto`.
"""

# ╔═╡ 554dc2f7-3e8e-5774-9ea1-01ff45507877
md"""
With dimensionless time $\tau=2\pi f_0t$, define the effective stiffness
$K_{\mathrm{eff}}=k+0.1\mu(k+1)$. The equations are

```math
\begin{aligned}
x' &= v,\\
v' &= -K_{\mathrm{eff}}x+(\mu-x^2)v-F.
\end{aligned}
```

The state is `u = [x, v]` and the parameters are `p = [μ, k, F]`.
For $K_{\mathrm{eff}}>0$, the equilibrium is $x_*=-F/K_{\mathrm{eff}}$, $v_*=0$.
Its linearized damping changes sign at $\mu-x_*^2=0$; a bias can therefore
shift the onset of oscillation. The scale $f_0$ converts simulation time to
seconds; the resulting pitch also depends on stiffness and drive.
"""

# ╔═╡ 832ca063-2534-5fd7-be72-71cdb3be0d11
function vocal!(du, u, p, t)
    x, v = u
    μ, k, F = p
    effective_stiffness = k + 0.1μ*(k + 1)
    du[1] = v
    du[2] = -effective_stiffness*x + (μ - x^2)*v - F
    nothing
end

# ╔═╡ bb6538ff-7b33-57a0-b1db-299b1510a9bc
md"""
Drive μ $(@bind vocal_drive Slider(-0.2:0.01:0.8; default=0.3, show_value=true))

Stiffness k $(@bind vocal_stiffness Slider(0.2:0.05:2.0; default=1.0, show_value=true))

Bias F $(@bind vocal_bias Slider(-0.3:0.01:0.3; default=0.0, show_value=true))

x(0) $(@bind vocal_initial_x Slider(-1.0:0.05:1.0; default=0.1, show_value=true))

v(0) $(@bind vocal_initial_v Slider(-1.0:0.05:1.0; default=0.1, show_value=true))

Pitch scale (Hz) $(@bind vocal_f0 Slider(80.0:5.0:320.0; default=180.0, show_value=true))

Duration (s) $(@bind vocal_duration Slider(0.5:0.5:4.0; default=2.0, show_value=true))

Move μ through zero with F = 0 to compare decay and self-oscillation.
Then increase the bias and observe how the equilibrium and the onset change.
"""

# ╔═╡ 8807d519-e599-5ddd-a88f-0e05a2f793fb
begin
    vocal_sample_rate = 48_000
    vocal_oversample = 4
end

# ╔═╡ 20b9c5aa-c4a4-5630-9208-7d052d079d01
vocal_u0 = [vocal_initial_x, vocal_initial_v]

# ╔═╡ d18e831f-339d-53e2-b0b8-af74b39db6d4
vocal_tspan = (0.0, 2pi*vocal_f0*vocal_duration)  # Dimensionless time τ

# ╔═╡ 6c722944-2106-56e1-ac71-2416aeb51f05
vocal_p = [vocal_drive, vocal_stiffness, vocal_bias]  # μ, k, F

# ╔═╡ a681a729-be38-5962-9450-1e9cd4083111
vocal_prob = ODEProblem(vocal!, vocal_u0, vocal_tspan, vocal_p)

# ╔═╡ 848c93c4-453b-5ebc-9fee-3db58a01d87b
vocal_sample_times = (0:round(Int, vocal_duration*vocal_sample_rate*vocal_oversample)-1) .* (2pi*vocal_f0/(vocal_sample_rate*vocal_oversample));

# ╔═╡ 93db7593-3eb5-581c-a0ae-202554f0c98f
vocal_sol = solve(vocal_prob, Tsit5(); abstol=1e-9, reltol=1e-7,
    saveat=vocal_sample_times, save_end=false, dense=false);

# ╔═╡ ee67ab91-77e6-5934-87e0-af2b1d0d9008
plot(vocal_sol.t[1:40:end] ./ (2pi*vocal_f0), vocal_sol[1,1:40:end];
    xlabel="time (s)", ylabel="fold displacement x", legend=false, margin=5 * Plots.mm)

# ╔═╡ 8c23dccf-a9fc-5f73-93e8-9d884a71aaeb
let window = max(1, length(vocal_sol.t)-9600):length(vocal_sol.t)
    waveform = plot(vocal_sol.t[window] ./ (2pi*vocal_f0), vocal_sol[1,window];
        xlabel="time (s)", ylabel="fold displacement x", legend=false)
    orbit = plot(vocal_sol[1,window], vocal_sol[2,window]; xlabel="x", ylabel="v", legend=false)
    plot(waveform, orbit; layout=(1,2), size=(850,320), margin=5 * Plots.mm)
end

# ╔═╡ 14cf8dd6-b5fe-5bd6-9254-2f2d8926f439
md"""
A growing transient followed by a closed phase-plane orbit indicates
self-sustained vibration. The restoring force uses $K_{\mathrm{eff}}$, so even
the undamped small-motion period is $2\pi/\sqrt{K_{\mathrm{eff}}}$ in $\tau$,
not $2\pi/\sqrt{k}$. Listen to $x$ as a simple voice-source proxy below.
"""

# ╔═╡ 2dc5caee-1015-53ea-8534-2ad399dc3ec6
vocal_raw_audio = vocal_sol[1,:];

# ╔═╡ a1702c11-10b2-5d42-ab48-7c5913861535
vocal_audio = prepare_audio(vocal_raw_audio, vocal_sample_rate; oversample=vocal_oversample);

# ╔═╡ 0452ddc4-0371-55fb-9b4a-eeeba2a3060e
vocal_wav_data = let buffer = IOBuffer()
    wavwrite(vocal_audio, buffer; Fs=vocal_sample_rate, nbits=16, compression=WAVE_FORMAT_PCM)
    take!(buffer)
end;

# ╔═╡ 16b8ab67-6feb-5e96-97e6-5e32fdf1811f
md"""
The solution is sampled at 4× the output rate, low-pass filtered and reduced
to 48 kHz, then centered, peak-normalized, and faded at the ends. The player
uses a sound-output proxy; use the unnormalized plots to compare amplitudes.
"""

# ╔═╡ d0271cf2-a49a-579e-b9fc-8feaedd32276
Resource("data:audio/wav;base64," * base64encode(vocal_wav_data), MIME("audio/wav"), ())

# ╔═╡ c2d52ecf-214e-51cd-83f2-afbd799c18c2
DownloadButton(vocal_wav_data, "vocal_folds.wav")

# ╔═╡ 9c428c74-5ad1-51d6-aa8b-226350f44f1f
md"""
### Relate phonation onset to the bifurcation lessons

**Predict → listen → explain:** vary net drive through onset while holding
stiffness and bias fixed. Compare the origin's local stability with the
late-time amplitude. Does the oscillation disappear smoothly or after a
finite-amplitude jump? Measure amplitude in the original state plots:
the normalized WAV makes quiet and loud solutions easier to hear, so its
playback level is not a measure of the model's physical amplitude.
"""

# ╔═╡ aac839e8-8ba3-5f9d-b9a5-dbfb15286f4c
md"""
## Birdsong: syrinx and delayed vocal tract

The syringeal labia are represented by a nonlinear oscillator with displacement
$x$ and velocity $y$. The controls $\alpha$ and $\beta$ change the pressure-like
drive and stiffness, while $\gamma$ sets the source time scale. The source
feeds a travelling pressure wave into a tract of length $L$; reflections return
after the round-trip time $T=2L/c$, with assumed sound speed $c=343$ m/s.
A linear acoustic circuit then filters the transmitted pressure.

This preserves both the source and the delayed tract from `DNL_Birdsong.jl` in
`NonlinearDynamicsPluto`. In this formulation the source drives the tract in
one direction: acoustic pressure does not feed back into the labial equation.
Time $t$ is in seconds.
"""

# ╔═╡ 10829ce5-2776-5c54-a34e-c374b0f4fda6
md"""
The source and the two travelling pressures satisfy

```math
\begin{aligned}
\dot{x}&=y,\\
\dot{y}&=\gamma^2(-\alpha-\beta x-x^3+x^2)-\gamma x(x+1)y,\\
P_i(t)&=\nu y(t)-rP_i(t-T),\\
P_0(t)&=(1-r)P_i(t-T/2).
\end{aligned}
```

Here $\nu$ converts source velocity to pressure, $r$ is the reflection
coefficient, and $P_0$ drives the acoustic circuit. Its coordinates are
$i_1$, $\Omega=\dot{i}_1$, and $i_3$; $i_3$ is the output-flow proxy.
Using the original circuit parameters, define

```math
A=\frac{1}{C_hM_G},\quad
B=R_h\left(\frac{1}{M_B}+\frac{1}{M_G}\right),\quad
C=A+\frac{R_hR_B}{M_GM_B},\quad
D=\frac{R_hR_B}{M_GM_B}.
```

The circuit equations are

```math
\begin{aligned}
\dot{i}_1&=\Omega,\\
\dot{\Omega}&=-Ai_1-B\Omega+Ci_3+DP_0+\frac{\dot{P}_0}{M_G},\\
\dot{i}_3&=-\frac{M_G}{M_B}\Omega-\frac{R_h}{M_B}i_3+\frac{P_0}{M_B}.
\end{aligned}
```

$C_h$ is the compliance, $M_G,M_B$ are inertance parameters, and $R_h,R_B$
set the circuit's resistive terms.
"""

# ╔═╡ 5837a2b5-5a17-53d3-a3e7-01f493a4dae3
md"""
To avoid differentiating a delayed algebraic pressure, set
$w=\Omega-P_0/M_G$. The same circuit becomes

```math
\begin{aligned}
\dot{i}_1&=w+P_0/M_G,\\
\dot{w}&=-Ai_1-Bw+Ci_3+(D-B/M_G)P_0,\\
\dot{i}_3&=-(M_G/M_B)w-(R_h/M_B)i_3.
\end{aligned}
```

The state is `u = [x, y, i1, w, i3, Pi, P0]`. The last two equations are
algebraic constraints, so their mass-matrix entries are zero. The other five
entries are one: `M*u′ = f(u, history, p, t)`.
This keeps the original pressure relations and both propagation delays.
"""

# ╔═╡ c1378b49-80da-55b1-a3fd-2ad6220c595a
function bird!(du, u, history, p, t)
    x, y, i1, w, i3, Pi, P0 = u
    α, β, γ, Ch, MG, MB, RB, Rh, r, ν, T = p
    A = 1/(Ch*MG)
    B = Rh*(1/MB + 1/MG)
    D = Rh*RB/(MG*MB)
    C = A + D
    du[1] = y
    du[2] = γ^2*(-α - β*x - x^3 + x^2) - γ*x*(x+1)*y
    du[3] = w + P0/MG
    du[4] = -A*i1 - B*w + C*i3 + (D - B/MG)*P0
    du[5] = -(MG/MB)*w - (Rh/MB)*i3
    du[6] = ν*y - r*history(p, t-T; idxs=6) - Pi
    du[7] = (1-r)*history(p, t-T/2; idxs=6) - P0
    nothing
end

# ╔═╡ 473f4746-a728-58e8-acc5-60a5f5e07fea
md"""
Pressure control α $(@bind bird_alpha Slider(0.0:0.01:0.2; default=0.1, show_value=true))

Stiffness control β $(@bind bird_beta Slider(0.0:0.01:0.3; default=0.1, show_value=true))

Source rate γ (s⁻¹) $(@bind bird_rate Slider(10000.0:500.0:35000.0; default=23500.0, show_value=true))

Reflection r $(@bind bird_reflection Slider(0.0:0.05:0.75; default=0.65, show_value=true))

Tract length L (m) $(@bind bird_tract_length Slider(0.015:0.005:0.045; default=0.025, show_value=true))

Duration (s) $(@bind bird_duration Slider(0.1:0.05:0.5; default=0.25, show_value=true))

The default α = β = 0.1 gives sustained phonation. Try α = β = 0.05 to
compare a decaying transient, then vary tract length or reflection to hear
the filtering while retaining the same source controls.
"""

# ╔═╡ ecbcb82e-02c8-50e4-a35b-7615eaaa833b
bird_u0 = [0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

# ╔═╡ 5c664ed5-d8f5-5b1f-8e5d-e004fa7d0de9
bird_p = [bird_alpha, bird_beta, bird_rate, 1.43e-10, 0.001, 1e4, 5e6, 2.4e4, bird_reflection, 0.001, 2*bird_tract_length/343]
# α, β, γ, Ch, MG, MB, RB, Rh, r, ν, T

# ╔═╡ a43c5b8f-9012-546b-a538-808f0bb8e185
bird_mass_matrix = Diagonal([1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0])

# ╔═╡ 3e38dca4-871c-529d-b748-9966b92b13d0
bird_system = DDEFunction(bird!; mass_matrix=bird_mass_matrix)

# ╔═╡ cf668145-7b1f-52cf-95d9-2d5ea0f99bdd
bird_lags = [bird_p[11]/2, bird_p[11]]  # One-way and round-trip times

# ╔═╡ 586e10c1-b4ff-5d65-8a6f-071965206705
md"""
Before t = 0 there is no travelling pressure in the tract. The initial
labial velocity and both pressures are zero, so the two algebraic constraints
are consistent with this history.
"""

# ╔═╡ ba718b70-b32e-58f9-bc9d-60826a091a5a
bird_history(p, t; idxs=nothing) = isnothing(idxs) ? zeros(7) : zeros(7)[idxs]

# ╔═╡ 872921a4-58ce-5d99-9468-7b2030eed726
bird_tspan = (0.0, bird_duration)

# ╔═╡ 6c4132ae-7777-5cb1-bcc5-aa21cb40d6c2
bird_prob = DDEProblem(bird_system, bird_u0, bird_history, bird_tspan, bird_p; constant_lags=bird_lags)

# ╔═╡ b37fabfd-6d22-538e-8d91-e6b507506452
bird_sol = solve(bird_prob, MethodOfSteps(Rodas5P());
    abstol=1e-10, reltol=1e-6, dt=1e-6);

# ╔═╡ 4d2b9e39-3b70-5e57-b6bc-5f2da805c671
let
    times = range(max(0.0, last(bird_sol.t)-0.01), last(bird_sol.t); length=2000)
    values = bird_sol(times)
    source_plot = plot(times, values[1,:]; xlabel="time (s)", ylabel="labial x", legend=false)
    output_plot = plot(times, values[5,:]; xlabel="time (s)", ylabel="output flow i₃", legend=false)
    orbit = plot(values[1,:], values[2,:] ./ bird_p[3];
        xlabel="x", ylabel="y/γ", legend=false)
    pressure_plot = plot(times, [values[6,:] values[7,:]];
        xlabel="time (s)", ylabel="pressure", label=["Pᵢ" "P₀"])
    plot(source_plot, output_plot, orbit, pressure_plot;
        layout=(2,2), size=(850,550), margin=5 * Plots.mm)
end

# ╔═╡ ad4923aa-dfa2-5b17-83a1-74bd4188cb9b
md"""
The source phase portrait shows the labial vibration; the pressure traces
show propagation and reflection. Listen to the filtered flow i₃ below. These
constant controls produce one sustained syllable rather than a complete song.
"""

# ╔═╡ b6ec4d1b-ec53-50c9-84c2-dadfa43a9b44
begin
    bird_sample_rate = 48_000
    bird_oversample = 4
end

# ╔═╡ 6f967fb5-ca3b-5116-a389-55781d6d8fc1
bird_sample_times = (0:round(Int, bird_duration*bird_sample_rate*bird_oversample)-1) ./ (bird_sample_rate*bird_oversample);

# ╔═╡ 0180dd9f-7849-56de-9936-dae5e2c0af78
bird_raw_audio = bird_sol(bird_sample_times; idxs=5).u;

# ╔═╡ 280aae98-0f95-5691-830e-49e07279d9b6
bird_audio = prepare_audio(bird_raw_audio, bird_sample_rate; oversample=bird_oversample);

# ╔═╡ ae8834a3-7fa0-5b1f-b43a-71c36d473380
bird_wav_data = let buffer = IOBuffer()
    wavwrite(bird_audio, buffer; Fs=bird_sample_rate, nbits=16, compression=WAVE_FORMAT_PCM)
    take!(buffer)
end;

# ╔═╡ e7fb6197-9331-59b1-b93b-12779e63dca6
md"""
The solution is sampled at 4× the output rate, low-pass filtered and reduced
to 48 kHz, then centered, peak-normalized, and faded at the ends. The player
uses a sound-output proxy; use the unnormalized plots to compare amplitudes.
"""

# ╔═╡ 87647d40-a09f-545b-aa72-295212fc758a
Resource("data:audio/wav;base64," * base64encode(bird_wav_data), MIME("audio/wav"), ())

# ╔═╡ ad4d55ab-02ef-59cb-bc89-7047dc76df97
DownloadButton(bird_wav_data, "birdsong.wav")

# ╔═╡ 440acc3d-bcee-548f-9bd7-0a381ea1feec
md"""
### Separate source dynamics from tract memory

**Predict → listen → explain:** keep the source parameters fixed and change
tract length, then reverse the experiment. Compare source and tract signals
and their spectra. The propagation delays store past states: a projection
onto two plotted variables is not the complete phase space of this system.
Use the source's ODE phase plane to study its onset, and the full delayed
solution to study the filtered sound and echoes.
"""

# ╔═╡ Cell order:
# ╠═0504907f-de44-547a-b346-fe28fbd333df
# ╟─5f8331a5-e154-5b74-a325-69b4f0355937
# ╟─460bdbd3-5ddf-5936-ac14-00b94f57c234
# ╟─063a0f42-b8ad-5bcf-804c-d8999e17445b
# ╟─554dc2f7-3e8e-5774-9ea1-01ff45507877
# ╟─832ca063-2534-5fd7-be72-71cdb3be0d11
# ╟─bb6538ff-7b33-57a0-b1db-299b1510a9bc
# ╟─8807d519-e599-5ddd-a88f-0e05a2f793fb
# ╠═20b9c5aa-c4a4-5630-9208-7d052d079d01
# ╟─d18e831f-339d-53e2-b0b8-af74b39db6d4
# ╠═6c722944-2106-56e1-ac71-2416aeb51f05
# ╟─a681a729-be38-5962-9450-1e9cd4083111
# ╟─848c93c4-453b-5ebc-9fee-3db58a01d87b
# ╟─93db7593-3eb5-581c-a0ae-202554f0c98f
# ╟─ee67ab91-77e6-5934-87e0-af2b1d0d9008
# ╟─8c23dccf-a9fc-5f73-93e8-9d884a71aaeb
# ╟─14cf8dd6-b5fe-5bd6-9254-2f2d8926f439
# ╠═2dc5caee-1015-53ea-8534-2ad399dc3ec6
# ╠═a1702c11-10b2-5d42-ab48-7c5913861535
# ╠═0452ddc4-0371-55fb-9b4a-eeeba2a3060e
# ╟─16b8ab67-6feb-5e96-97e6-5e32fdf1811f
# ╟─d0271cf2-a49a-579e-b9fc-8feaedd32276
# ╟─c2d52ecf-214e-51cd-83f2-afbd799c18c2
# ╟─9c428c74-5ad1-51d6-aa8b-226350f44f1f
# ╟─aac839e8-8ba3-5f9d-b9a5-dbfb15286f4c
# ╟─10829ce5-2776-5c54-a34e-c374b0f4fda6
# ╟─5837a2b5-5a17-53d3-a3e7-01f493a4dae3
# ╟─c1378b49-80da-55b1-a3fd-2ad6220c595a
# ╟─473f4746-a728-58e8-acc5-60a5f5e07fea
# ╠═ecbcb82e-02c8-50e4-a35b-7615eaaa833b
# ╠═5c664ed5-d8f5-5b1f-8e5d-e004fa7d0de9
# ╟─a43c5b8f-9012-546b-a538-808f0bb8e185
# ╟─3e38dca4-871c-529d-b748-9966b92b13d0
# ╟─cf668145-7b1f-52cf-95d9-2d5ea0f99bdd
# ╟─586e10c1-b4ff-5d65-8a6f-071965206705
# ╟─ba718b70-b32e-58f9-bc9d-60826a091a5a
# ╟─872921a4-58ce-5d99-9468-7b2030eed726
# ╟─6c4132ae-7777-5cb1-bcc5-aa21cb40d6c2
# ╟─b37fabfd-6d22-538e-8d91-e6b507506452
# ╟─4d2b9e39-3b70-5e57-b6bc-5f2da805c671
# ╟─ad4923aa-dfa2-5b17-83a1-74bd4188cb9b
# ╟─b6ec4d1b-ec53-50c9-84c2-dadfa43a9b44
# ╟─6f967fb5-ca3b-5116-a389-55781d6d8fc1
# ╠═0180dd9f-7849-56de-9936-dae5e2c0af78
# ╠═280aae98-0f95-5691-830e-49e07279d9b6
# ╠═ae8834a3-7fa0-5b1f-b43a-71c36d473380
# ╟─e7fb6197-9331-59b1-b93b-12779e63dca6
# ╟─87647d40-a09f-545b-aa72-295212fc758a
# ╟─ad4d55ab-02ef-59cb-bc89-7047dc76df97
# ╟─440acc3d-bcee-548f-9bd7-0a381ea1feec
