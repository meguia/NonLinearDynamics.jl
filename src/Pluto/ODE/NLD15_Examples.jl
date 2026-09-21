### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 12e69ec2-2f96-5d53-9a85-3c18fe0618e0
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using NonLinearDynamics: prepare_audio
    using DifferentialEquations, LinearAlgebra, Plots, PlutoUI, WAV, Base64
end

# ╔═╡ 03c3f57a-1d63-5c3d-9fc0-b220ed20f6b9
md"""
# Examples: vocal folds and birdsong

Two models connect nonlinear oscillation to vocal sound. Begin with phonation
onset and saturation in a reduced vocal-fold oscillator, then follow a
syringeal source through a delayed vocal tract and an acoustic filter.
"""

# ╔═╡ 9a290496-dcb3-534b-b5c1-1e317c9ac5c2
TableOfContents()

# ╔═╡ b59b87e5-6d40-5f1b-a49b-0cd99569b9c9
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

# ╔═╡ 7f2ae8cd-3a0b-553e-9344-636e343b4c61
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

# ╔═╡ b35620c1-9444-5507-88da-2340c1ea1133
function vocal!(du, u, p, t)
    x, v = u
    μ, k, F = p
    effective_stiffness = k + 0.1μ*(k + 1)
    du[1] = v
    du[2] = -effective_stiffness*x + (μ - x^2)*v - F
    nothing
end

# ╔═╡ dca290a5-6866-5a3d-81d2-fcf057c762c7
md"""
Define the initial condition, time span, and parameters, build the problem,
then solve it with `Tsit5()`. Edit the values below to explore phonation onset.
"""

# ╔═╡ b476c770-ca32-5857-9c9c-978d81332775
vocal_f0 = 180.0  # Frequency scale in Hz

# ╔═╡ d9b6f231-c004-57fd-a090-5ea02b735fd2
vocal_duration = 2.0  # Seconds

# ╔═╡ 501c59bf-4155-566d-b0bb-d77457891544
begin
    vocal_sample_rate = 48_000
    vocal_oversample = 4
end

# ╔═╡ deed8ac1-6c08-544e-b87f-88d8cdeee7a6
vocal_u0 = [0.1, 0.1]  # x(0), v(0)

# ╔═╡ 4032c6b0-c19f-514f-88b3-b15c89efdb0b
vocal_tspan = (0.0, 2pi*vocal_f0*vocal_duration)  # Dimensionless time τ

# ╔═╡ 78a7b560-c87e-5555-ad00-78e1873af736
vocal_p = [1.3, 3.0, 0.0]  # μ, k, F

# ╔═╡ 353887b4-7cb3-5fb6-bb17-64a0489333e8
vocal_prob = ODEProblem(vocal!, vocal_u0, vocal_tspan, vocal_p)

# ╔═╡ dd1f4147-8810-513e-92ea-232031bd2910
vocal_sample_times = (0:round(Int, vocal_duration*vocal_sample_rate*vocal_oversample)-1) .* (2pi*vocal_f0/(vocal_sample_rate*vocal_oversample));

# ╔═╡ ba232a4f-55db-51ef-a64b-5b42f560f352
vocal_sol = solve(vocal_prob, Tsit5(); abstol=1e-9, reltol=1e-7,
    saveat=vocal_sample_times, save_end=false, dense=false);

# ╔═╡ 33b6239d-17e8-5b87-96ae-69088657413f
plot(vocal_sol.t[1:40:end] ./ (2pi*vocal_f0), vocal_sol[1,1:40:end];
    xlabel="time (s)", ylabel="fold displacement x", legend=false, margin=5 * Plots.mm)

# ╔═╡ c88c95f2-06ef-5e0f-b7e8-afbb47d9c6f8
let window = max(1, length(vocal_sol.t)-9600):length(vocal_sol.t)
    waveform = plot(vocal_sol.t[window] ./ (2pi*vocal_f0), vocal_sol[1,window];
        xlabel="time (s)", ylabel="fold displacement x", legend=false)
    orbit = plot(vocal_sol[1,window], vocal_sol[2,window]; xlabel="x", ylabel="v", legend=false)
    plot(waveform, orbit; layout=(1,2), size=(850,320), margin=5 * Plots.mm)
end

# ╔═╡ 9bf5facb-1f05-59b9-82ea-33dbc01bbaca
md"""
A growing transient followed by a closed phase-plane orbit indicates
self-sustained vibration. The restoring force uses $K_{\mathrm{eff}}$, so even
the undamped small-motion period is $2\pi/\sqrt{K_{\mathrm{eff}}}$ in $\tau$,
not $2\pi/\sqrt{k}$. Listen to $x$ as a simple voice-source proxy below.
"""

# ╔═╡ 8ea963e5-cadc-59f0-a5f1-c04ba7193626
vocal_raw_audio = vocal_sol[1,:];

# ╔═╡ df7b3c7e-43d4-573d-ba63-83cf07ce1696
vocal_audio = prepare_audio(vocal_raw_audio, vocal_sample_rate; oversample=vocal_oversample);

# ╔═╡ 123441d5-bcdd-544a-b0aa-15f5a319bc83
vocal_wav_data = let buffer = IOBuffer()
    wavwrite(vocal_audio, buffer; Fs=vocal_sample_rate, nbits=16, compression=WAVE_FORMAT_PCM)
    take!(buffer)
end;

# ╔═╡ d50b4e14-89f9-55ae-9c32-a174ba9961c2
md"""
The solution is sampled at 4× the output rate, low-pass filtered and reduced
to 48 kHz, then centered, peak-normalized, and faded at the ends. The player
uses a sound-output proxy; use the unnormalized plots to compare amplitudes.
"""

# ╔═╡ 5d691b2f-f016-53a2-a4a9-4c8602d49552
Resource("data:audio/wav;base64," * base64encode(vocal_wav_data), MIME("audio/wav"), ())

# ╔═╡ 2402302f-859e-56b5-b5fb-63c3ca4474a0
DownloadButton(vocal_wav_data, "vocal_folds.wav")

# ╔═╡ 18dcee18-1042-5ef5-811b-6f2fc0bfabc1
md"""
### Relate phonation onset to the bifurcation lessons

**Predict → listen → explain:** vary net drive through onset while holding
stiffness and bias fixed. Compare the origin's local stability with the
late-time amplitude. Does the oscillation disappear smoothly or after a
finite-amplitude jump? Measure amplitude in the original state plots:
the normalized WAV makes quiet and loud solutions easier to hear, so its
playback level is not a measure of the model's physical amplitude.
"""

# ╔═╡ 8eddbb24-cf07-5ffc-848a-3c2f67bc869b
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

# ╔═╡ 6d9fa23c-3f29-547f-930b-48aa3bbe2c1e
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

# ╔═╡ 749509b6-cde1-50f5-893d-351993684606
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

# ╔═╡ d511bf6c-30c6-5f18-9246-d3528e1d02cd
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

# ╔═╡ aa388e93-d4f8-5767-ad9f-c26cd0b08018
md"""
Propagation requires a `DDEProblem`: add a history function and the two
constant delays to the usual function → initial condition → parameters → problem → solve workflow.
Use `MethodOfSteps(Rodas5P())` because the acoustic circuit is stiff and has algebraic constraints.
"""

# ╔═╡ 42854355-c30a-5097-9f10-503a86439f6e
bird_duration = 0.25  # Seconds

# ╔═╡ 53a487b7-0bf0-5c6a-8777-080c8cfdbfa8
bird_u0 = [0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

# ╔═╡ 27a6d723-2ed9-5fcd-9495-c4ea45d163bc
bird_p = [0.1, 0.1, 23500.0, 1.43e-10, 0.001, 1e4, 5e6, 2.4e4, 0.65, 0.001, 2*0.025/343]
# α, β, γ, Ch, MG, MB, RB, Rh, r, ν, T

# ╔═╡ e5ea78fa-03e3-55fb-935d-9c5cf1d36e32
bird_mass_matrix = Diagonal([1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0])

# ╔═╡ 5041f362-3deb-50d6-a412-97e64c658b99
bird_system = DDEFunction(bird!; mass_matrix=bird_mass_matrix)

# ╔═╡ 187d31be-3a11-510c-b9d4-0a8822acd380
bird_lags = [bird_p[11]/2, bird_p[11]]  # One-way and round-trip times

# ╔═╡ 8e7e5933-bf53-5b26-8fac-c6195000fc85
md"""
Before t = 0 there is no travelling pressure in the tract. The initial
labial velocity and both pressures are zero, so the two algebraic constraints
are consistent with this history.
"""

# ╔═╡ 824f4713-86b6-5d8b-b569-b596a0a14246
bird_history(p, t; idxs=nothing) = isnothing(idxs) ? zeros(7) : zeros(7)[idxs]

# ╔═╡ 6d9c4f51-a8e0-5159-b300-ca8e649ee6ba
bird_tspan = (0.0, bird_duration)

# ╔═╡ 4a8909d4-c3be-5d5d-bb29-713465544614
bird_prob = DDEProblem(bird_system, bird_u0, bird_history, bird_tspan, bird_p; constant_lags=bird_lags)

# ╔═╡ 7540c981-1712-5ed8-8f4a-b48ed4fcfb6f
bird_sol = solve(bird_prob, MethodOfSteps(Rodas5P());
    abstol=1e-10, reltol=1e-6, dt=1e-6);

# ╔═╡ bed9d115-a1e6-5b0c-a24b-0ca4204eca95
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

# ╔═╡ cb5a6eae-968c-54f5-b7bd-138efa51530b
md"""
The source phase portrait shows the labial vibration; the pressure traces
show propagation and reflection. Listen to the filtered flow i₃ below. These
constant controls produce one sustained syllable rather than a complete song.
"""

# ╔═╡ abaaf933-2bf3-597a-8057-700fd95faf43
begin
    bird_sample_rate = 48_000
    bird_oversample = 4
end

# ╔═╡ eb729210-3a51-5f44-ab31-409a43fb5df2
bird_sample_times = (0:round(Int, bird_duration*bird_sample_rate*bird_oversample)-1) ./ (bird_sample_rate*bird_oversample);

# ╔═╡ e0362ad0-8bff-54e9-957e-4fd9f79b810d
bird_raw_audio = bird_sol(bird_sample_times; idxs=5).u;

# ╔═╡ 49b8df1b-7f17-5d8f-9825-439142fad86f
bird_audio = prepare_audio(bird_raw_audio, bird_sample_rate; oversample=bird_oversample);

# ╔═╡ 37ce00ab-375d-50fb-b273-724cacf3d31f
bird_wav_data = let buffer = IOBuffer()
    wavwrite(bird_audio, buffer; Fs=bird_sample_rate, nbits=16, compression=WAVE_FORMAT_PCM)
    take!(buffer)
end;

# ╔═╡ 1046e011-47ea-5e4e-ba06-d1b8c74d5963
md"""
The solution is sampled at 4× the output rate, low-pass filtered and reduced
to 48 kHz, then centered, peak-normalized, and faded at the ends. The player
uses a sound-output proxy; use the unnormalized plots to compare amplitudes.
"""

# ╔═╡ 634a7abb-2bf6-5221-90f1-60f5e1b7bdd1
Resource("data:audio/wav;base64," * base64encode(bird_wav_data), MIME("audio/wav"), ())

# ╔═╡ 6e8708d5-50fb-5724-88a9-266d62daf5fd
DownloadButton(bird_wav_data, "birdsong.wav")

# ╔═╡ d61f0370-89c6-59b4-8b3f-eaf92dd434a6
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
# ╠═12e69ec2-2f96-5d53-9a85-3c18fe0618e0
# ╟─03c3f57a-1d63-5c3d-9fc0-b220ed20f6b9
# ╟─9a290496-dcb3-534b-b5c1-1e317c9ac5c2
# ╟─b59b87e5-6d40-5f1b-a49b-0cd99569b9c9
# ╟─7f2ae8cd-3a0b-553e-9344-636e343b4c61
# ╠═b35620c1-9444-5507-88da-2340c1ea1133
# ╟─dca290a5-6866-5a3d-81d2-fcf057c762c7
# ╠═b476c770-ca32-5857-9c9c-978d81332775
# ╠═d9b6f231-c004-57fd-a090-5ea02b735fd2
# ╠═501c59bf-4155-566d-b0bb-d77457891544
# ╠═deed8ac1-6c08-544e-b87f-88d8cdeee7a6
# ╠═4032c6b0-c19f-514f-88b3-b15c89efdb0b
# ╠═78a7b560-c87e-5555-ad00-78e1873af736
# ╠═353887b4-7cb3-5fb6-bb17-64a0489333e8
# ╠═dd1f4147-8810-513e-92ea-232031bd2910
# ╠═ba232a4f-55db-51ef-a64b-5b42f560f352
# ╠═33b6239d-17e8-5b87-96ae-69088657413f
# ╠═c88c95f2-06ef-5e0f-b7e8-afbb47d9c6f8
# ╟─9bf5facb-1f05-59b9-82ea-33dbc01bbaca
# ╠═8ea963e5-cadc-59f0-a5f1-c04ba7193626
# ╠═df7b3c7e-43d4-573d-ba63-83cf07ce1696
# ╠═123441d5-bcdd-544a-b0aa-15f5a319bc83
# ╟─d50b4e14-89f9-55ae-9c32-a174ba9961c2
# ╟─5d691b2f-f016-53a2-a4a9-4c8602d49552
# ╟─2402302f-859e-56b5-b5fb-63c3ca4474a0
# ╟─18dcee18-1042-5ef5-811b-6f2fc0bfabc1
# ╟─8eddbb24-cf07-5ffc-848a-3c2f67bc869b
# ╟─6d9fa23c-3f29-547f-930b-48aa3bbe2c1e
# ╟─749509b6-cde1-50f5-893d-351993684606
# ╠═d511bf6c-30c6-5f18-9246-d3528e1d02cd
# ╟─aa388e93-d4f8-5767-ad9f-c26cd0b08018
# ╠═42854355-c30a-5097-9f10-503a86439f6e
# ╠═53a487b7-0bf0-5c6a-8777-080c8cfdbfa8
# ╠═27a6d723-2ed9-5fcd-9495-c4ea45d163bc
# ╠═e5ea78fa-03e3-55fb-935d-9c5cf1d36e32
# ╠═5041f362-3deb-50d6-a412-97e64c658b99
# ╠═187d31be-3a11-510c-b9d4-0a8822acd380
# ╟─8e7e5933-bf53-5b26-8fac-c6195000fc85
# ╠═824f4713-86b6-5d8b-b569-b596a0a14246
# ╠═6d9c4f51-a8e0-5159-b300-ca8e649ee6ba
# ╠═4a8909d4-c3be-5d5d-bb29-713465544614
# ╠═7540c981-1712-5ed8-8f4a-b48ed4fcfb6f
# ╠═bed9d115-a1e6-5b0c-a24b-0ca4204eca95
# ╟─cb5a6eae-968c-54f5-b7bd-138efa51530b
# ╠═abaaf933-2bf3-597a-8057-700fd95faf43
# ╠═eb729210-3a51-5f44-ab31-409a43fb5df2
# ╠═e0362ad0-8bff-54e9-957e-4fd9f79b810d
# ╠═49b8df1b-7f17-5d8f-9825-439142fad86f
# ╠═37ce00ab-375d-50fb-b273-724cacf3d31f
# ╟─1046e011-47ea-5e4e-ba06-d1b8c74d5963
# ╟─634a7abb-2bf6-5221-90f1-60f5e1b7bdd1
# ╟─6e8708d5-50fb-5724-88a9-266d62daf5fd
# ╟─d61f0370-89c6-59b4-8b3f-eaf92dd434a6
