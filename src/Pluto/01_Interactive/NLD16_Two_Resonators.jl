### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ c54d8cf4-1652-52bf-a0dc-8c3b46699d62
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots, PlutoUI, DSP, WAV, Base64
    using NonLinearDynamics: prepare_audio
end

# ╔═╡ 7bd05538-4ac4-5288-841f-c6a9be240c1b
md"""
# Nonlinear reed oscillator coupled to two resonators

Two acoustic modes interact through the reed and its nonlinear airflow, allowing a periodic tone or a two-frequency oscillation.
"""

# ╔═╡ 3f8538a3-e648-58f7-9c98-f76884432be8
md"""
```math
\begin{aligned}
P &= p_1+p_2, & \Delta &= \gamma-P, & U &= \zeta\max(h,0)\frac{\Delta}{(\Delta^2+\epsilon^2)^{1/4}},\\
h' &= v, & v' &= -d_r\Omega_r v-\Omega_r^2(h-1+\Delta),\\
q_1' &= p_1, & p_1' &= -q_1-\frac{p_1}{Q_1}+F_1U,\\
q_2' &= p_2, & p_2' &= -\Omega_2^2q_2-\frac{\Omega_2}{Q_2}p_2+F_2U.
\end{aligned}
```

h is the reed opening, p₁ and p₂ are acoustic pressures, and q₁ and q₂ are their time integrals.
Constant blowing pressure γ supplies energy through the nonlinear airflow U; both resonators feed pressure back to the reed.
Primes use ``τ=2πf_0t`` and the first resonance has ``Ω_1=1``.
The small ε smooths the square-root slope at zero pressure drop. Closure stops flow; reed-contact forces are omitted.
"""

# ╔═╡ 16c4fa9b-1bb1-56e4-82bc-277ac2474474
function model!(du, u, p, t)
    h, v, q1, p1, q2, p2 = u
    γ, ζ, Ωr, dr, Ω2, Q1, Q2, F1, F2, ε = p
    Δ = γ - (p1+p2)
    flow = ζ*max(h, 0)*Δ/sqrt(sqrt(Δ^2+ε^2))
    du[1] = v
    du[2] = -dr*Ωr*v - Ωr^2*(h-1+Δ)
    du[3] = p1
    du[4] = -q1 - p1/Q1 + F1*flow
    du[5] = p2
    du[6] = -Ω2^2*q2 - Ω2*p2/Q2 + F2*flow
    nothing
end

# ╔═╡ 1da1d949-c4ab-57f1-9e10-3eedf6d8c983
md"""
Blowing preset $(@bind blowing Select([0.45 => "Periodic (γ = 0.45)", 0.47 => "Quasiperiodic (γ = 0.47)"]; default=0.47))

Second resonance Ω₂ $(@bind ratio Slider(2.0:0.01:2.6; default=2.5, show_value=true))

Pitch scale (Hz) $(@bind f0 Slider(110.0:1.0:330.0; default=185.0, show_value=true))

Keep Ω₂ = 2.5 to compare the two presets, then explore the resonance spacing.
"""

# ╔═╡ 04408a27-c9f5-59c4-9a30-a9a83dc792c0
# Static equilibrium at the default parameters, with p₁ perturbed by 0.01.
u0 = [0.53, 0.0, 0.11582631989741293, 0.01, 0.03344769734042085, 0.0]

# ╔═╡ ab104a33-e5bf-5de2-999e-6b75cb5d75e5
begin
    duration = 4.0
    sample_rate = 48_000
    oversample = 4
    transient = 5000.0  # Dimensionless time discarded before listening and analysis
end

# ╔═╡ b21ee3a4-43cc-5ec8-a667-d5d0fcfadcc8
tspan = (0.0, transient + 2pi*f0*duration)

# ╔═╡ 2d69a2fa-f094-5862-aa02-b7da6d6f8fec
p = [blowing, 0.28, 2pi*2000/1161.2, 1.0, ratio, 36.6, 41.2, 1322/1161.2, 2386/1161.2, 1e-4]
# γ, ζ, Ωr, dr, Ω2, Q1, Q2, F1, F2, ε

# ╔═╡ 149c7b46-be88-55ba-bada-6679b9deb453
prob = ODEProblem(model!, u0, tspan, p)

# ╔═╡ ce94914d-8037-59e4-aed9-dbc8583bb0a1
sample_times = transient .+ (0:round(Int, duration*sample_rate*oversample)-1) .* (2pi*f0/(sample_rate*oversample));

# ╔═╡ 75d0629c-0f08-5cbe-82af-2b90994dbcb6
sol = solve(prob, Tsit5(); abstol=1e-10, reltol=1e-8, maxiters=10^7,
    saveat=sample_times, save_start=false, save_end=false, dense=false);

# ╔═╡ b0bb733e-5e51-5808-a376-8012e569b788
raw_audio = sol[4,:] .+ sol[6,:];

# ╔═╡ 67d20c1f-205a-5146-9c35-77900edc2ed9
md"""
After discarding τ < 5000, a periodic orbit gives isolated Poincaré points; a two-frequency torus gives a closed curve.
Compare the section with the spectrum: two peaks alone do not establish quasiperiodicity, and the oscillation frequencies need not equal the bare resonances.
"""

# ╔═╡ 85d3d97d-1774-5316-82c2-c3a26af6f3b6
section = let
    keep = findall(t -> t >= 5000.0, sol.t)
    a, b, c = sol[3,:][keep], sol[5,:][keep], sol[6,:][keep]
    times = sol.t[keep]
    level = sum(a)/length(a)
    crossings = [i for i in 1:length(a)-1 if a[i] <= level < a[i+1]]
    fraction = [(level-a[i])/(a[i+1]-a[i]) for i in crossings]
    (q2 = [b[i]+s*(b[i+1]-b[i]) for (i,s) in zip(crossings,fraction)],
     p2 = [c[i]+s*(c[i+1]-c[i]) for (i,s) in zip(crossings,fraction)],
     times = [times[i]+s*(times[i+1]-times[i]) for (i,s) in zip(crossings,fraction)],
     level = level)
end;

# ╔═╡ d6ebb728-77bb-5eb7-8898-4ed6c14652eb
let
    window = findall(t -> t >= last(sol.t)-80.0, sol.t)
    waveform = plot((sol.t[window] .- first(sol.t[window])) ./ (2pi*f0), raw_audio[window];
        xlabel="time (s)", ylabel="pressure p₁ + p₂", legend=false)
    section_plot = scatter(section.q2, section.p2; markersize=2, markerstrokewidth=0,
        xlabel="q₂", ylabel="p₂", title="Poincaré: q₁ = mean(q₁), p₁ > 0", legend=false)
    plot(waveform, section_plot; layout=(1,2), size=(850,320))
end

# ╔═╡ 3068f59a-2fa0-58bf-8f53-329c3fb90fa4
let
    spectrum = periodogram(raw_audio; fs=sample_rate*oversample, window=hanning)
    relative_power = spectrum.power ./ max(maximum(spectrum.power), eps(Float64))
    plot(spectrum.freq, 10 .* log10.(max.(relative_power, 1e-10));
        xlims=(0, 6f0), ylims=(-80, 5), xlabel="frequency (Hz)", ylabel="relative power (dB)", legend=false)
end

# ╔═╡ 4a889945-7778-53db-b56a-a3d0c1f9d32d
audio = prepare_audio(raw_audio, sample_rate; oversample);

# ╔═╡ 1f98fd30-53c2-511d-97d1-d8206e3e29f7
wav_data = let buffer = IOBuffer()
    wavwrite(audio, buffer; Fs=sample_rate, nbits=16, compression=WAVE_FORMAT_PCM)
    take!(buffer)
end;

# ╔═╡ f3e14eb9-82d8-5f27-8161-4999b285d0db
md"""
Listen to the sum of the two modal pressures after the transient.
The 4× sampled solution is low-pass filtered to 48 kHz, peak-normalized, and faded at both ends.
"""

# ╔═╡ d9e8385b-fcb0-55d2-a6e0-e3d3a3baaa1a
Resource("data:audio/wav;base64," * base64encode(wav_data), MIME("audio/wav"), ())

# ╔═╡ 14a5ef05-7ce0-5f2a-96a1-bc32bd31fad3
DownloadButton(wav_data, "two_resonators.wav")

# ╔═╡ 750638f4-6eb3-5be5-ba2e-74308e40bb4c
md"""
The bore parameters are based on [Doc, Vergez and Missoum (2014)](https://doi.org/10.3813/AAA.918734), with detuned resonances and stronger reed damping for this teaching example.
"""

# ╔═╡ Cell order:
# ╠═c54d8cf4-1652-52bf-a0dc-8c3b46699d62
# ╟─7bd05538-4ac4-5288-841f-c6a9be240c1b
# ╟─3f8538a3-e648-58f7-9c98-f76884432be8
# ╠═16c4fa9b-1bb1-56e4-82bc-277ac2474474
# ╟─1da1d949-c4ab-57f1-9e10-3eedf6d8c983
# ╠═04408a27-c9f5-59c4-9a30-a9a83dc792c0
# ╟─ab104a33-e5bf-5de2-999e-6b75cb5d75e5
# ╟─b21ee3a4-43cc-5ec8-a667-d5d0fcfadcc8
# ╠═2d69a2fa-f094-5862-aa02-b7da6d6f8fec
# ╟─149c7b46-be88-55ba-bada-6679b9deb453
# ╟─ce94914d-8037-59e4-aed9-dbc8583bb0a1
# ╟─75d0629c-0f08-5cbe-82af-2b90994dbcb6
# ╟─b0bb733e-5e51-5808-a376-8012e569b788
# ╟─67d20c1f-205a-5146-9c35-77900edc2ed9
# ╟─85d3d97d-1774-5316-82c2-c3a26af6f3b6
# ╟─d6ebb728-77bb-5eb7-8898-4ed6c14652eb
# ╟─3068f59a-2fa0-58bf-8f53-329c3fb90fa4
# ╟─4a889945-7778-53db-b56a-a3d0c1f9d32d
# ╟─1f98fd30-53c2-511d-97d1-d8206e3e29f7
# ╟─f3e14eb9-82d8-5f27-8161-4999b285d0db
# ╟─d9e8385b-fcb0-55d2-a6e0-e3d3a3baaa1a
# ╟─14a5ef05-7ce0-5f2a-96a1-bc32bd31fad3
# ╟─750638f4-6eb3-5be5-ba2e-74308e40bb4c
