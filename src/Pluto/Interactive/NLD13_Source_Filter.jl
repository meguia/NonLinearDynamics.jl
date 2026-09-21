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

# ╔═╡ 7b801516-e388-5740-8b96-6ac2afe673ec
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using NonLinearDynamics: prepare_audio
    using DifferentialEquations, Plots, PlutoUI, WAV, Base64
end

# ╔═╡ 644bc060-a288-5c73-a7c9-7118eeaf3058
md"""
# Source–filter models

A nonlinear source exchanges energy with a damped resonator. Compare reed,
bowed, and struck excitation: the first two can sustain oscillation, while
the strike supplies its energy through the initial condition. In these
coupled models, the resonator also acts back on the source.
"""

# ╔═╡ 0df328c3-dce4-506a-9a10-0d81e9431d5a
TableOfContents()

# ╔═╡ cc6008d5-43fb-5573-9f71-0357022f5f45
md"""
## Reed oscillator coupled to a bore resonance

The nonlinear damping supplies energy at small velocity and removes it at large velocity.
The coupling forces are equal and opposite, so the resonator also acts back on the source.

```math
\begin{aligned}
x' &= v, & v' &= -x+\mu(1-v^2)v+\kappa(q-x),\\
q' &= w, & w' &= -\Omega^2q-2\zeta\Omega w+\kappa(x-q).
\end{aligned}
```

With dimensionless time τ, the total energy satisfies
``E'=\mu(v^2-v^4)-2\zeta\Omega w^2``. Blowing supplies energy and bore losses remove it.
This is a teaching model of excitation and resonance, with one mode and no reed-contact or air-flow calculation.
"""

# ╔═╡ 63eeec30-d79c-5393-a833-0b60ba752742
function reed!(du, u, p, t)
    x, v, q, w = u
    μ, κ, Ω, ζ = p
    du[1] = v
    du[2] = -x + μ*(1-v^2)*v + κ*(q-x)
    du[3] = w
    du[4] = -Ω^2*q - 2ζ*Ω*w + κ*(x-q)
    nothing
end

# ╔═╡ 07e9803f-d742-5f20-8517-584ed8216fd1
md"""
μ $(@bind reed_drive Slider(0.0:0.05:0.8; default=0.3, show_value=true))

Coupling κ $(@bind reed_coupling Slider(0.0:0.02:0.3; default=0.12, show_value=true))

Resonator ratio Ω $(@bind reed_ratio Slider(0.5:0.02:3.0; default=1.02, show_value=true))

Pitch scale (Hz) $(@bind reed_f0 Slider(110.0:1.0:330.0; default=220.0, show_value=true))

Lower μ until oscillations die out, then vary Ω to hear the effect of the resonator.
"""

# ╔═╡ 03083e00-4677-5bbd-a5c8-a9737e680413
begin
    reed_duration = 4.0
    reed_sample_rate = 48_000
    reed_oversample = 4
end

# ╔═╡ 2d3f7ae1-4a97-5f21-b4d1-be9e54b2e728
reed_u0 = [0.05, 0.0, 0.0, 0.0]

# ╔═╡ 904063aa-5ca4-5afa-9783-1cc4f9c5123b
reed_tspan = (0.0, 2pi*reed_f0*reed_duration)  # Dimensionless time τ

# ╔═╡ 0bc5c78a-eb57-5380-85e9-fe6880dfac3f
reed_p = [reed_drive, reed_coupling, reed_ratio, 0.025]  # μ, κ, Ω, ζ

# ╔═╡ ad10e763-0b83-5f29-b7ab-f05a6fdc9b0d
reed_prob = ODEProblem(reed!, reed_u0, reed_tspan, reed_p)

# ╔═╡ bbc9bde3-a3c5-58ba-9a1a-2326345ae228
reed_sample_times = (0:round(Int, reed_duration*reed_sample_rate*reed_oversample)-1) .* (2pi*reed_f0/(reed_sample_rate*reed_oversample));

# ╔═╡ abab7691-cabb-52b8-8203-60e256467c52
reed_sol = solve(reed_prob, Tsit5(); abstol=1e-9, reltol=1e-7,
    saveat=reed_sample_times, save_end=false, dense=false);

# ╔═╡ 18540938-de0a-5162-9538-701c475fb13b
let window = max(1, length(reed_sol.t)-4800):length(reed_sol.t)
    waveform = plot(reed_sol.t[window] ./ (2pi*reed_f0), reed_sol[3,:][window];
        xlabel="time (s)", ylabel="resonator q", legend=false)
    orbit = plot(reed_sol[1,:][window], reed_sol[2,:][window];
        xlabel="source x", ylabel="source v", legend=false)
    plot(waveform, orbit; layout=(1,2), size=(800,300))
end

# ╔═╡ 6734c4c3-4e8d-53c4-b72c-442f6e9ea64d
reed_raw_audio = reed_sol[3,:];

# ╔═╡ 23557e4c-8447-5782-98e8-4614cf89248d
let block = round(Int, 0.02*reed_sample_rate*reed_oversample)
    n = length(reed_raw_audio) ÷ block
    envelope = [sqrt(sum(abs2, view(reed_raw_audio, (k-1)*block+1:k*block))/block) for k in 1:n]
    times = ((0:n-1) .+ 0.5) .* (block/(reed_sample_rate*reed_oversample))
    plot(times, envelope; xlabel="time (s)", ylabel="resonator RMS", legend=false)
end

# ╔═╡ 670122bc-c363-551e-b2e3-f1939c9e61d6
reed_audio = prepare_audio(reed_raw_audio, reed_sample_rate; oversample=reed_oversample);

# ╔═╡ b14922b2-7c94-545c-916f-5be9562d4295
reed_wav_data = let buffer = IOBuffer()
    wavwrite(reed_audio, buffer; Fs=reed_sample_rate, nbits=16, compression=WAVE_FORMAT_PCM)
    take!(buffer)
end;

# ╔═╡ b324eeb9-ef55-5fc1-b578-e1b88f065a8b
md"""
Listen to the resonator coordinate q, with ``τ=2πf_0 t`` setting the pitch scale.
The 4× sampled solution is low-pass filtered and reduced to 48 kHz, then
peak-normalized with short fades; compare the plotted amplitudes for loudness.
"""

# ╔═╡ a513a223-adf0-5bb8-a30a-897806f43905
Resource("data:audio/wav;base64," * base64encode(reed_wav_data), MIME("audio/wav"), ())

# ╔═╡ df80b8d3-6420-5c65-9822-d55f4e7492bf
DownloadButton(reed_wav_data, "reed_resonator.wav")

# ╔═╡ 9f117512-16b4-5978-bc1d-ce3363fc7240
md"""
### Connect the sound to the phase plane

**Predict → listen → explain:** first reduce the coupling and identify the
self-oscillator's limit cycle. Restore the coupling and move the resonator
frequency through the source frequency. Compare displacement, radiated
output, and spectral peaks. A resonator can emphasize a frequency without
supplying energy; identify the nonlinear term that sustains oscillation.
"""

# ╔═╡ 5399e2f3-e291-58d3-8cb8-eecab2804042
md"""
## Bowed oscillator coupled to a body resonance

The bow moves at speed V, its force changes sign with the relative speed and weakens during fast sliding.
The smoothing scale ε keeps the equations continuous, so this model approximates sliding friction rather than enforcing exact sticking.

```math
\begin{aligned}
x' &= v, & v' &= -x-\delta v+F_b(V-v)+\kappa(q-x),\\
q' &= w, & w' &= -\Omega^2q-2\zeta\Omega w+\kappa(x-q),\\
F_b(s) &= F\,\frac{\tanh(s/\epsilon)}{1+(s/v_s)^2}.
\end{aligned}
```

The bow supplies work through ``F_b(V-v)v``; the body acts back on the string through the coupling spring.
One string mode and one body mode illustrate friction-driven oscillation without representing a complete violin.
"""

# ╔═╡ 8c3e256d-d516-5c81-adc5-61b32112ed55
function bowed!(du, u, p, t)
    x, v, q, w = u
    F, V, ε, vs, δ, κ, Ω, ζ = p
    du[1] = v
    du[2] = -x - δ*v + F*tanh((V-v)/ε)/(1+((V-v)/vs)^2) + κ*(q-x)
    du[3] = w
    du[4] = -Ω^2*q - 2ζ*Ω*w + κ*(x-q)
    nothing
end

# ╔═╡ 9c1103ad-9477-5da9-8618-a265dab90aa7
md"""
V $(@bind bowed_drive Slider(0.1:0.05:0.6; default=0.3, show_value=true))

Coupling κ $(@bind bowed_coupling Slider(0.0:0.02:0.3; default=0.08, show_value=true))

Resonator ratio Ω $(@bind bowed_ratio Slider(0.5:0.02:3.0; default=1.4, show_value=true))

Pitch scale (Hz) $(@bind bowed_f0 Slider(110.0:1.0:330.0; default=196.0, show_value=true))

Change the bow speed V and body frequency Ω, then compare the waveform and its sound.
"""

# ╔═╡ 804fdef5-91f3-544c-8a97-8aa30ea12738
begin
    bowed_duration = 4.0
    bowed_sample_rate = 48_000
    bowed_oversample = 4
end

# ╔═╡ 670572a6-ebb9-5d65-a2f1-750987a631c5
bowed_u0 = [0.0, 0.0, 0.0, 0.0]

# ╔═╡ 26f5a210-2596-5046-81f1-f12a0b1285e1
bowed_tspan = (0.0, 2pi*bowed_f0*bowed_duration)  # Dimensionless time τ

# ╔═╡ f1b3328e-c3a5-5ba0-b23b-98b956872f8d
bowed_p = [0.7, bowed_drive, 0.02, 0.25, 0.02, bowed_coupling, bowed_ratio, 0.02]  # F, V, ε, vs, δ, κ, Ω, ζ

# ╔═╡ 50e48c6f-5c8a-5e89-b66d-f1f369277c90
bowed_prob = ODEProblem(bowed!, bowed_u0, bowed_tspan, bowed_p)

# ╔═╡ 095309b8-9b13-51e6-a020-3f1db7e9a026
bowed_sample_times = (0:round(Int, bowed_duration*bowed_sample_rate*bowed_oversample)-1) .* (2pi*bowed_f0/(bowed_sample_rate*bowed_oversample));

# ╔═╡ 13f02eca-5bef-5d4b-a095-d821edd555bd
bowed_sol = solve(bowed_prob, Tsit5(); abstol=1e-9, reltol=1e-7,
    saveat=bowed_sample_times, save_end=false, dense=false);

# ╔═╡ a7f9fd5a-0ec9-547f-b1d0-8eb441d867d4
let window = max(1, length(bowed_sol.t)-9600):length(bowed_sol.t)
    waveform = plot(bowed_sol.t[window] ./ (2pi*bowed_f0), bowed_sol[3,:][window];
        xlabel="time (s)", ylabel="resonator q", legend=false)
    orbit = plot(bowed_sol[1,:][window], bowed_sol[2,:][window];
        xlabel="source x", ylabel="source v", legend=false)
    plot(waveform, orbit; layout=(1,2), size=(800,300))
end

# ╔═╡ 2264a931-173d-5c72-a9ea-15e01f436f05
bowed_raw_audio = bowed_sol[3,:];

# ╔═╡ d8a00dfe-3923-5fb5-a3ae-e50551ac108a
let block = round(Int, 0.02*bowed_sample_rate*bowed_oversample)
    n = length(bowed_raw_audio) ÷ block
    envelope = [sqrt(sum(abs2, view(bowed_raw_audio, (k-1)*block+1:k*block))/block) for k in 1:n]
    times = ((0:n-1) .+ 0.5) .* (block/(bowed_sample_rate*bowed_oversample))
    plot(times, envelope; xlabel="time (s)", ylabel="resonator RMS", legend=false)
end

# ╔═╡ 4979235f-c765-5336-bcba-066e8cc08d16
bowed_audio = prepare_audio(bowed_raw_audio, bowed_sample_rate; oversample=bowed_oversample);

# ╔═╡ 0e8d1c6f-81d6-5b39-b36f-98996b645d18
bowed_wav_data = let buffer = IOBuffer()
    wavwrite(bowed_audio, buffer; Fs=bowed_sample_rate, nbits=16, compression=WAVE_FORMAT_PCM)
    take!(buffer)
end;

# ╔═╡ cda8c010-e2a4-5926-aa38-031a87ef9acf
md"""
Listen to the resonator coordinate q, with ``τ=2πf_0 t`` setting the pitch scale.
The 4× sampled solution is low-pass filtered and reduced to 48 kHz, then
peak-normalized with short fades; compare the plotted amplitudes for loudness.
"""

# ╔═╡ 65ed3282-b664-5cb1-9763-76ac03117357
Resource("data:audio/wav;base64," * base64encode(bowed_wav_data), MIME("audio/wav"), ())

# ╔═╡ 49ab782f-f01f-5ad2-adcf-e91a8889ff67
DownloadButton(bowed_wav_data, "bowed_resonator.wav")

# ╔═╡ 01db4f31-0407-5409-b3a6-58bec7cc7455
md"""
### Connect friction to self-oscillation

**Predict → listen → explain:** locate the part of the friction curve that
supplies energy, then compare it with the bow and oscillator velocities.
Vary bow speed and coupling separately. Does the strongest spectral peak
follow the source, the body resonance, or a harmonic? The smooth friction
law approximates slipping; exact sticking would require an additional rule.
"""

# ╔═╡ fe002ece-6262-5b8a-8b6a-373073f06d46
md"""
## Struck nonlinear oscillator and resonator

An initial velocity gives the source a single strike, with no further external forcing.
The cubic spring stiffens at large displacement, and the second mode introduces another resonant frequency.

```math
\begin{aligned}
x' &= v, & v' &= -x-\alpha x^3-2\delta v+\kappa(q-x),\\
q' &= w, & w' &= -\Omega^2q-2\zeta\Omega w+\kappa(x-q).
\end{aligned}
```

For positive α and κ, the energy
``E=(v^2+x^2+w^2+\Omega^2q^2+\kappa(x-q)^2)/2+\alpha x^4/4``
satisfies ``E'=-2\delta v^2-2\zeta\Omega w^2\leq0``.
This is a two-mode teaching model for a struck bar or plate, not a complete bell-impact model.
"""

# ╔═╡ 392a981c-f44e-53a0-91b5-fcda561b3643
function struck!(du, u, p, t)
    x, v, q, w = u
    α, δ, κ, Ω, ζ = p
    du[1] = v
    du[2] = -x - α*x^3 - 2δ*v + κ*(q-x)
    du[3] = w
    du[4] = -Ω^2*q - 2ζ*Ω*w + κ*(x-q)
    nothing
end

# ╔═╡ 09003168-2ece-5b03-994f-df323eab6d3c
md"""
α $(@bind struck_drive Slider(0.0:0.1:1.5; default=0.8, show_value=true))

Coupling κ $(@bind struck_coupling Slider(0.0:0.02:0.3; default=0.08, show_value=true))

Resonator ratio Ω $(@bind struck_ratio Slider(0.5:0.02:3.0; default=2.4, show_value=true))

Pitch scale (Hz) $(@bind struck_f0 Slider(110.0:1.0:330.0; default=220.0, show_value=true))

Increase α to hear the pitch change during the decay, and vary Ω to compare harmonic and inharmonic resonances.
"""

# ╔═╡ 789aea5e-9156-5142-98b9-1a45f8867ab8
begin
    struck_duration = 4.0
    struck_sample_rate = 48_000
    struck_oversample = 4
end

# ╔═╡ b6d52cd7-e312-5807-ad98-8f9eed2e28ed
struck_u0 = [0.0, 1.0, 0.0, 0.0]

# ╔═╡ 1a49cdc0-7e14-5702-8292-20d60f4d688c
struck_tspan = (0.0, 2pi*struck_f0*struck_duration)  # Dimensionless time τ

# ╔═╡ 6aafba6c-9100-54ee-8f5e-38f7bc1e6b0d
struck_p = [struck_drive, 0.001, struck_coupling, struck_ratio, 0.0015]  # α, δ, κ, Ω, ζ

# ╔═╡ ddc9b0d5-eee8-541e-93c0-a47bb23863b4
struck_prob = ODEProblem(struck!, struck_u0, struck_tspan, struck_p)

# ╔═╡ e075e7b8-4007-5e83-ab16-9ad0b7c90492
struck_sample_times = (0:round(Int, struck_duration*struck_sample_rate*struck_oversample)-1) .* (2pi*struck_f0/(struck_sample_rate*struck_oversample));

# ╔═╡ ebc4d026-0b43-5c3b-8c99-46414649d968
struck_sol = solve(struck_prob, Tsit5(); abstol=1e-9, reltol=1e-7,
    saveat=struck_sample_times, save_end=false, dense=false);

# ╔═╡ 72c35fb1-46bb-566b-9494-713c1be745d6
let window = 1:min(9601, length(struck_sol.t))
    waveform = plot(struck_sol.t[window] ./ (2pi*struck_f0), struck_sol[3,:][window];
        xlabel="time (s)", ylabel="resonator q", legend=false)
    orbit = plot(struck_sol[1,:][window], struck_sol[2,:][window];
        xlabel="source x", ylabel="source v", legend=false)
    plot(waveform, orbit; layout=(1,2), size=(800,300))
end

# ╔═╡ 76424dc2-4374-586e-8ef8-e286bcb24571
struck_raw_audio = struck_sol[3,:];

# ╔═╡ 5d0c1755-385d-5bc4-8bdb-f61976ad9635
let block = round(Int, 0.02*struck_sample_rate*struck_oversample)
    n = length(struck_raw_audio) ÷ block
    envelope = [sqrt(sum(abs2, view(struck_raw_audio, (k-1)*block+1:k*block))/block) for k in 1:n]
    times = ((0:n-1) .+ 0.5) .* (block/(struck_sample_rate*struck_oversample))
    plot(times, envelope; xlabel="time (s)", ylabel="resonator RMS", legend=false)
end

# ╔═╡ 62783ade-795f-5cb8-9af1-245c6e7683b7
struck_audio = prepare_audio(struck_raw_audio, struck_sample_rate; oversample=struck_oversample);

# ╔═╡ 4234c63c-b5df-5db7-95d3-609800eea854
struck_wav_data = let buffer = IOBuffer()
    wavwrite(struck_audio, buffer; Fs=struck_sample_rate, nbits=16, compression=WAVE_FORMAT_PCM)
    take!(buffer)
end;

# ╔═╡ d5b1bb31-34dd-592a-a86a-413cbe172d45
md"""
Listen to the resonator coordinate q, with ``τ=2πf_0 t`` setting the pitch scale.
The 4× sampled solution is low-pass filtered and reduced to 48 kHz, then
peak-normalized with short fades; compare the plotted amplitudes for loudness.
"""

# ╔═╡ 87cb660a-a189-580f-b180-a484e0ba9376
Resource("data:audio/wav;base64," * base64encode(struck_wav_data), MIME("audio/wav"), ())

# ╔═╡ 91718610-02a8-55a5-af4f-d70671a92272
DownloadButton(struck_wav_data, "struck_resonator.wav")

# ╔═╡ 01f8b3ee-9818-5525-a5e5-535f30dc6210
md"""
### A transient is not an attracting cycle

**Predict → listen → explain:** increase the strike while keeping the
parameters fixed. Follow the changing instantaneous oscillation period and
the decaying envelope. Compare this with the sustained reed example.
Two peaks in a transient spectrum can represent two decaying modes; they do
not establish a quasiperiodic attractor.
"""

# ╔═╡ Cell order:
# ╠═7b801516-e388-5740-8b96-6ac2afe673ec
# ╟─644bc060-a288-5c73-a7c9-7118eeaf3058
# ╟─0df328c3-dce4-506a-9a10-0d81e9431d5a
# ╟─cc6008d5-43fb-5573-9f71-0357022f5f45
# ╠═63eeec30-d79c-5393-a833-0b60ba752742
# ╟─07e9803f-d742-5f20-8517-584ed8216fd1
# ╟─18540938-de0a-5162-9538-701c475fb13b
# ╟─03083e00-4677-5bbd-a5c8-a9737e680413
# ╠═2d3f7ae1-4a97-5f21-b4d1-be9e54b2e728
# ╟─904063aa-5ca4-5afa-9783-1cc4f9c5123b
# ╠═0bc5c78a-eb57-5380-85e9-fe6880dfac3f
# ╟─ad10e763-0b83-5f29-b7ab-f05a6fdc9b0d
# ╟─bbc9bde3-a3c5-58ba-9a1a-2326345ae228
# ╟─abab7691-cabb-52b8-8203-60e256467c52
# ╟─6734c4c3-4e8d-53c4-b72c-442f6e9ea64d
# ╟─23557e4c-8447-5782-98e8-4614cf89248d
# ╟─670122bc-c363-551e-b2e3-f1939c9e61d6
# ╟─b14922b2-7c94-545c-916f-5be9562d4295
# ╟─b324eeb9-ef55-5fc1-b578-e1b88f065a8b
# ╟─a513a223-adf0-5bb8-a30a-897806f43905
# ╟─df80b8d3-6420-5c65-9822-d55f4e7492bf
# ╟─9f117512-16b4-5978-bc1d-ce3363fc7240
# ╟─5399e2f3-e291-58d3-8cb8-eecab2804042
# ╠═8c3e256d-d516-5c81-adc5-61b32112ed55
# ╟─9c1103ad-9477-5da9-8618-a265dab90aa7
# ╠═a7f9fd5a-0ec9-547f-b1d0-8eb441d867d4
# ╟─804fdef5-91f3-544c-8a97-8aa30ea12738
# ╠═670572a6-ebb9-5d65-a2f1-750987a631c5
# ╟─26f5a210-2596-5046-81f1-f12a0b1285e1
# ╠═f1b3328e-c3a5-5ba0-b23b-98b956872f8d
# ╟─50e48c6f-5c8a-5e89-b66d-f1f369277c90
# ╟─095309b8-9b13-51e6-a020-3f1db7e9a026
# ╟─13f02eca-5bef-5d4b-a095-d821edd555bd
# ╟─2264a931-173d-5c72-a9ea-15e01f436f05
# ╟─d8a00dfe-3923-5fb5-a3ae-e50551ac108a
# ╟─4979235f-c765-5336-bcba-066e8cc08d16
# ╟─0e8d1c6f-81d6-5b39-b36f-98996b645d18
# ╟─cda8c010-e2a4-5926-aa38-031a87ef9acf
# ╟─65ed3282-b664-5cb1-9763-76ac03117357
# ╟─49ab782f-f01f-5ad2-adcf-e91a8889ff67
# ╟─01db4f31-0407-5409-b3a6-58bec7cc7455
# ╟─fe002ece-6262-5b8a-8b6a-373073f06d46
# ╠═392a981c-f44e-53a0-91b5-fcda561b3643
# ╟─09003168-2ece-5b03-994f-df323eab6d3c
# ╟─789aea5e-9156-5142-98b9-1a45f8867ab8
# ╠═b6d52cd7-e312-5807-ad98-8f9eed2e28ed
# ╟─1a49cdc0-7e14-5702-8292-20d60f4d688c
# ╠═6aafba6c-9100-54ee-8f5e-38f7bc1e6b0d
# ╟─ddc9b0d5-eee8-541e-93c0-a47bb23863b4
# ╟─e075e7b8-4007-5e83-ab16-9ad0b7c90492
# ╟─ebc4d026-0b43-5c3b-8c99-46414649d968
# ╠═72c35fb1-46bb-566b-9494-713c1be745d6
# ╟─76424dc2-4374-586e-8ef8-e286bcb24571
# ╟─5d0c1755-385d-5bc4-8bdb-f61976ad9635
# ╟─62783ade-795f-5cb8-9af1-245c6e7683b7
# ╟─4234c63c-b5df-5db7-95d3-609800eea854
# ╟─d5b1bb31-34dd-592a-a86a-413cbe172d45
# ╟─87cb660a-a189-580f-b180-a484e0ba9376
# ╟─91718610-02a8-55a5-af4f-d70671a92272
# ╟─01f8b3ee-9818-5525-a5e5-535f30dc6210
