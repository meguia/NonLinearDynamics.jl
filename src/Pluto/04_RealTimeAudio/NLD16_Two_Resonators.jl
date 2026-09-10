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

# ╔═╡ 8225fd70-94c3-50ae-8c68-c787f4d16b32
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots, PlutoUI, DSP
    using RealTimeAudioDiffEq, PlutoHooks
end

# ╔═╡ 52cd1bec-3895-5b08-8827-77f3bf7e1dab
md"""
# Nonlinear reed oscillator coupled to two resonators: live audio

Two acoustic modes interact through the reed and its nonlinear airflow, allowing a periodic tone or a two-frequency oscillation.

Choose Play to start; after the attack, change the blowing preset to hear the transition.
"""

# ╔═╡ bb6fdb95-5a64-59a2-bbf4-f60ce8df9edc
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

# ╔═╡ 625dcba9-0680-52b3-a12a-7d0d51de9927
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

# ╔═╡ f3e270d1-127f-50b3-a72c-f5398fd60e8d
# Static equilibrium at the default parameters, with p₁ perturbed by 0.01.
u0 = [0.53, 0.0, 0.11582631989741293, 0.01, 0.03344769734042085, 0.0]

# ╔═╡ d1190997-d723-5251-920c-c9b894b1a22e
p = [0.47, 0.28, 2pi*2000/1161.2, 1.0, 2.5, 36.6, 41.2, 1322/1161.2, 2386/1161.2, 1e-4]
# γ, ζ, Ωr, dr, Ω2, Q1, Q2, F1, F2, ε

# ╔═╡ d57b16a7-8efe-5236-ae1d-12f21d848f58
source = let
    audio_source = DESource(model!, copy(u0), copy(p); alg=Tsit5(), channel_map=[[4,6], [4,6]])
    # Version 0.1 keeps solver tolerances in the source problem; set them before playback.
    audio_source.data.problem = remake(audio_source.data.problem; abstol=1e-9, reltol=1e-7)
    audio_source
end;

# ╔═╡ 9f0a9b29-5a74-5659-b209-131cce4633fe
md"""
Blowing preset $(@bind live_parameter Select([0.45 => "Periodic (γ = 0.45)", 0.47 => "Quasiperiodic (γ = 0.47)"]; default=0.47))

Second resonance Ω₂ $(@bind ratio Slider(2.0:0.01:2.6; default=2.5, show_value=true))

Pitch scale (Hz) $(@bind f0 Slider(110.0:1.0:330.0; default=185.0, show_value=true))

Keep Ω₂ = 2.5 to compare the two presets, then explore the resonance spacing.

Gain $(@bind gain Slider(0.0:0.01:0.5; default=0.05, show_value=true))
"""

# ╔═╡ befb536e-4d2d-57ad-99d6-eba1de009abf
controls = begin
    set_param!(source, 1, Float64(live_parameter))
    set_param!(source, 5, Float64(ratio))
    set_ts!(source, 2pi*f0)
    set_gain!(source, Float64(gain))
    nothing
end

# ╔═╡ 990825b3-74c4-5581-bfb4-b40cbc5118f5
# Solve a preview through the transient before starting live playback.
preview = let
    controls
    prob = ODEProblem(model!, copy(u0), (0.0, 10000.0), copy(get_params(source)))
    solve(prob, Tsit5(); saveat=0.05, abstol=1e-10, reltol=1e-8, maxiters=10^7)
end;

# ╔═╡ 05c861a2-aadc-5ab6-9eec-10c1d9f59514
pressure = preview[4,:] .+ preview[6,:];

# ╔═╡ 07d2b877-6739-562d-9efb-584c67db5ff7
md"""
After discarding τ < 5000, a periodic orbit gives isolated Poincaré points; a two-frequency torus gives a closed curve.
Compare the section with the spectrum: two peaks alone do not establish quasiperiodicity, and the oscillation frequencies need not equal the bare resonances.
"""

# ╔═╡ 5407b598-deb0-5bed-85d2-255f43c0921f
section = let
    keep = findall(t -> t >= 5000.0, preview.t)
    a, b, c = preview[3,:][keep], preview[5,:][keep], preview[6,:][keep]
    times = preview.t[keep]
    level = sum(a)/length(a)
    crossings = [i for i in 1:length(a)-1 if a[i] <= level < a[i+1]]
    fraction = [(level-a[i])/(a[i+1]-a[i]) for i in crossings]
    (q2 = [b[i]+s*(b[i+1]-b[i]) for (i,s) in zip(crossings,fraction)],
     p2 = [c[i]+s*(c[i+1]-c[i]) for (i,s) in zip(crossings,fraction)],
     times = [times[i]+s*(times[i+1]-times[i]) for (i,s) in zip(crossings,fraction)],
     level = level)
end;

# ╔═╡ 2d9d5722-06b6-59b1-91ea-335e6628327d
let
    window = findall(t -> t >= last(preview.t)-80.0, preview.t)
    waveform = plot((preview.t[window] .- first(preview.t[window])) ./ (2pi*f0), pressure[window];
        xlabel="time (s)", ylabel="pressure p₁ + p₂", legend=false)
    section_plot = scatter(section.q2, section.p2; markersize=2, markerstrokewidth=0,
        xlabel="q₂", ylabel="p₂", title="Poincaré: q₁ = mean(q₁), p₁ > 0", legend=false)
    plot(waveform, section_plot; layout=(1,2), size=(850,320))
end

# ╔═╡ f1e8eb10-c560-5516-bc10-25c1bca7315e
let
    spectrum = periodogram(pressure[findall(t -> t >= 5000.0, preview.t)]; fs=2pi*f0/0.05, window=hanning)
    relative_power = spectrum.power ./ max(maximum(spectrum.power), eps(Float64))
    plot(spectrum.freq, 10 .* log10.(max.(relative_power, 1e-10));
        xlims=(0, 6f0), ylims=(-80, 5), xlabel="frequency (Hz)", ylabel="relative power (dB)", legend=false)
end

# ╔═╡ d70e1a3c-e834-5248-9452-678f052cc089
md"""
Play $(@bind playing CheckBox(default=false))
"""

# ╔═╡ b2565914-7b3e-59f5-b495-a919cc7f326d
@use_effect([source, playing]) do
    if playing
        device = get_default_output_device()
        device < 0 && error("No output device is available. Use the WAV example or select a device with list_devices().")
        reset_state!(source)
        start_DESource(source, device; buffer_size=UInt32(2048))
    end
    return () -> begin
        isactive(source) && stop_DESource(source)
    end
end

# ╔═╡ 972cdaab-3fd0-599a-bfec-f1bcb37e5a6d
md"""
Uncheck Play to stop. Both channels contain p₁ + p₂. Live gain controls the output level; the offline notebook also exports a filtered, normalized WAV.
"""

# ╔═╡ 4404e3d1-bf7b-56bb-a4cb-6e8b375f24e3
md"""
The bore parameters are based on [Doc, Vergez and Missoum (2014)](https://doi.org/10.3813/AAA.918734), with detuned resonances and stronger reed damping for this teaching example.
"""

# ╔═╡ Cell order:
# ╠═8225fd70-94c3-50ae-8c68-c787f4d16b32
# ╟─52cd1bec-3895-5b08-8827-77f3bf7e1dab
# ╟─bb6fdb95-5a64-59a2-bbf4-f60ce8df9edc
# ╠═625dcba9-0680-52b3-a12a-7d0d51de9927
# ╠═f3e270d1-127f-50b3-a72c-f5398fd60e8d
# ╠═d1190997-d723-5251-920c-c9b894b1a22e
# ╠═d57b16a7-8efe-5236-ae1d-12f21d848f58
# ╟─9f0a9b29-5a74-5659-b209-131cce4633fe
# ╠═befb536e-4d2d-57ad-99d6-eba1de009abf
# ╠═990825b3-74c4-5581-bfb4-b40cbc5118f5
# ╠═05c861a2-aadc-5ab6-9eec-10c1d9f59514
# ╟─07d2b877-6739-562d-9efb-584c67db5ff7
# ╟─5407b598-deb0-5bed-85d2-255f43c0921f
# ╠═2d9d5722-06b6-59b1-91ea-335e6628327d
# ╠═f1e8eb10-c560-5516-bc10-25c1bca7315e
# ╟─d70e1a3c-e834-5248-9452-678f052cc089
# ╠═b2565914-7b3e-59f5-b495-a919cc7f326d
# ╟─972cdaab-3fd0-599a-bfec-f1bcb37e5a6d
# ╟─4404e3d1-bf7b-56bb-a4cb-6e8b375f24e3
