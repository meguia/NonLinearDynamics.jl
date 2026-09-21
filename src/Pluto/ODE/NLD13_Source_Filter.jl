### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 23bfa7b2-9b8f-572f-bb20-7901fdd84161
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using NonLinearDynamics: prepare_audio
    using DifferentialEquations, Plots, PlutoUI, WAV, Base64
end

# ╔═╡ dd87b0f0-fdf3-5858-8058-b730a11360ff
md"""
# Source–filter models

A nonlinear source exchanges energy with a damped resonator. Compare reed,
bowed, and struck excitation: the first two can sustain oscillation, while
the strike supplies its energy through the initial condition. In these
coupled models, the resonator also acts back on the source.
"""

# ╔═╡ 349c8bfa-8317-5164-9118-f762ee1abcbd
TableOfContents()

# ╔═╡ 5813756a-2aa6-5c5a-99bb-bffd1996ed6e
md"""
## Reed oscillator coupled to a bore resonance

A Rayleigh source exchanges energy with one damped bore mode; the resonator coordinate q is a simple sound-output proxy.
"""

# ╔═╡ b925cdd8-32c6-5320-ad37-d304bec8de47
md"""
```math
\begin{aligned}
x' &= v, & v' &= -x+\mu(1-v^2)v+\kappa(q-x),\\
q' &= w, & w' &= -\Omega^2q-2\zeta\Omega w+\kappa(x-q).
\end{aligned}
```

Primes denote derivatives with respect to dimensionless time ``τ=2πf_0 t``.
"""

# ╔═╡ e8db5332-0778-57eb-a419-89cbb06dbd75
function reed!(du, u, p, t)
    x, v, q, w = u
    μ, κ, Ω, ζ = p
    du[1] = v
    du[2] = -x + μ*(1-v^2)*v + κ*(q-x)
    du[3] = w
    du[4] = -Ω^2*q - 2ζ*Ω*w + κ*(x-q)
    nothing
end

# ╔═╡ efa1b4d5-aeef-5eb6-b5fe-aa43374f2f31
reed_f0 = 440.0  # Frequency scale in Hz

# ╔═╡ e1b94c94-e7cb-5c7d-8512-9ec0410d3084
reed_duration = 4.0  # Seconds

# ╔═╡ 1e981501-2642-562a-b065-cf820f51d0e6
begin
    reed_sample_rate = 48_000
    reed_oversample = 4
end

# ╔═╡ 11e9129e-2d52-5140-a3ba-924b85bd3dd9
reed_u0 = [0.05, 0.0, 0.0, 0.0]

# ╔═╡ bdff993b-83a6-56a9-9d8c-24494f798851
reed_tspan = (0.0, 2pi*reed_f0*reed_duration)  # Dimensionless time τ

# ╔═╡ 9d49ecc1-fb2d-5e2d-8c9f-1c6184e1a0cc
reed_p = [7, 0.1, 0.91, 0.15]  # μ, κ, Ω, ζ

# ╔═╡ b206bd76-4107-5393-8144-2e649c090168
reed_prob = ODEProblem(reed!, reed_u0, reed_tspan, reed_p)

# ╔═╡ 2f05e51c-6681-5fb6-a938-b3c132bdc48d
reed_sample_times = (0:round(Int, reed_duration*reed_sample_rate*reed_oversample)-1) .* (2pi*reed_f0/(reed_sample_rate*reed_oversample));

# ╔═╡ b88d1334-d885-5400-8496-e39e69572428
reed_sol = solve(reed_prob, Tsit5(); abstol=1e-9, reltol=1e-7,
    saveat=reed_sample_times, save_end=false, dense=false);

# ╔═╡ 32f6fba9-e83e-5ebc-9aab-d63a9ffe6c44
let window = max(1, length(reed_sol.t)-9600):length(reed_sol.t)
    waveform = plot(reed_sol.t[window] ./ (2pi*reed_f0), reed_sol[3,:][window];
        xlabel="time (s)", ylabel="resonator q", legend=false)
    orbit = plot(reed_sol[1,:][window], reed_sol[2,:][window];
        xlabel="source x", ylabel="source v", legend=false)
    plot(waveform, orbit; layout=(1,2), size=(800,300))
end

# ╔═╡ 3ba443fb-d410-5de3-96ac-3e5a34834947
reed_raw_audio = reed_sol[3,:];

# ╔═╡ aab39aee-3460-5140-bd34-0a339eeaa92e
let block = round(Int, 0.02*reed_sample_rate*reed_oversample)
    n = length(reed_raw_audio) ÷ block
    envelope = [sqrt(sum(abs2, view(reed_raw_audio, (k-1)*block+1:k*block))/block) for k in 1:n]
    times = ((0:n-1) .+ 0.5) .* (block/(reed_sample_rate*reed_oversample))
    plot(times, envelope; xlabel="time (s)", ylabel="resonator RMS", legend=false)
end

# ╔═╡ 71505b27-7e34-5f55-933e-73c67189ba9a
reed_audio = prepare_audio(reed_raw_audio, reed_sample_rate; oversample=reed_oversample);

# ╔═╡ 4acfca4f-b51b-5d85-9e6f-c033dce73e71
reed_wav_data = let buffer = IOBuffer()
    wavwrite(reed_audio, buffer; Fs=reed_sample_rate, nbits=16, compression=WAVE_FORMAT_PCM)
    take!(buffer)
end;

# ╔═╡ 4dc10528-93ec-5041-b9a4-f52b88aa76af
md"""
Listen to the resonator coordinate q, with ``τ=2πf_0 t`` setting the pitch scale.
The 4× sampled solution is low-pass filtered and reduced to 48 kHz, then
peak-normalized with short fades; compare the plotted amplitudes for loudness.
"""

# ╔═╡ f417718f-736f-56db-bcee-f95ebf0e2b2b
Resource("data:audio/wav;base64," * base64encode(reed_wav_data), MIME("audio/wav"), ())

# ╔═╡ 3c7ac3bf-e652-5eab-95e0-37c9a9ebd2d6
DownloadButton(reed_wav_data, "reed_resonator.wav")

# ╔═╡ 6ebb53f3-e3c5-5a18-8803-74a21a1c5e89
md"""
### Connect the sound to the phase plane

**Predict → listen → explain:** first reduce the coupling and identify the
self-oscillator's limit cycle. Restore the coupling and move the resonator
frequency through the source frequency. Compare displacement, radiated
output, and spectral peaks. A resonator can emphasize a frequency without
supplying energy; identify the nonlinear term that sustains oscillation.
"""

# ╔═╡ 6b883fff-6f7e-598f-9b45-c3e69fd54f5f
md"""
## Bowed oscillator coupled to a body resonance

A smooth sliding-friction law drives a string mode coupled to a damped body mode; q is the sound-output proxy.
"""

# ╔═╡ ab523375-2d5c-535f-863c-e5a13cf8d419
md"""
```math
\begin{aligned}
x' &= v, & v' &= -x-\delta v+F_b(V-v)+\kappa(q-x),\\
q' &= w, & w' &= -\Omega^2q-2\zeta\Omega w+\kappa(x-q),\\
F_b(s) &= F\frac{\tanh(s/\epsilon)}{1+(s/v_s)^2}.
\end{aligned}
```

Primes denote derivatives with respect to dimensionless time ``τ=2πf_0 t``.
"""

# ╔═╡ 6aee948d-dde3-5758-8902-2c7cd7e54150
function bowed!(du, u, p, t)
    x, v, q, w = u
    F, V, ε, vs, δ, κ, Ω, ζ = p
    du[1] = v
    du[2] = -x - δ*v + F*tanh((V-v)/ε)/(1+((V-v)/vs)^2) + κ*(q-x)
    du[3] = w
    du[4] = -Ω^2*q - 2ζ*Ω*w + κ*(x-q)
    nothing
end

# ╔═╡ 16b18521-fabb-5281-9297-0131a250779d
bowed_f0 = 196.0  # Frequency scale in Hz

# ╔═╡ b6f85afb-123a-5180-97b5-30d29bb18319
bowed_duration = 6.0  # Seconds

# ╔═╡ b80c51cf-fd17-5770-a8d2-1998fdf8dfb4
begin
    bowed_sample_rate = 48_000
    bowed_oversample = 4
end

# ╔═╡ 0a58b2ce-b9e5-5be0-ba18-72f82774de03
bowed_u0 = [0.0, 0.0, 0.0, 0.0]

# ╔═╡ 2ebe1af9-99ab-5021-95f7-9190a6d01914
bowed_tspan = (0.0, 2pi*bowed_f0*bowed_duration)  # Dimensionless time τ

# ╔═╡ d89d73ab-d021-597a-bfb6-aaf696e82be4
bowed_p = [0.7, 0.3, 0.02, 0.25, 0.02, 0.01, 1.4, 0.02]  # F, V, ε, vs, δ, κ, Ω, ζ

# ╔═╡ ccd6bd9c-77af-5636-9c0c-907fa640b216
bowed_prob = ODEProblem(bowed!, bowed_u0, bowed_tspan, bowed_p)

# ╔═╡ 3e6e6228-9607-575c-9f7b-881cdef11c8f
bowed_sample_times = (0:round(Int, bowed_duration*bowed_sample_rate*bowed_oversample)-1) .* (2pi*bowed_f0/(bowed_sample_rate*bowed_oversample));

# ╔═╡ 842144ba-df5f-5024-8d57-35824957a2f2
bowed_sol = solve(bowed_prob, Tsit5(); abstol=1e-9, reltol=1e-7,
    saveat=bowed_sample_times, save_end=false, dense=false);

# ╔═╡ 1029fda6-0708-571d-83e3-93f6e0f3e2ec
let window = max(1, length(bowed_sol.t)-9600):length(bowed_sol.t)
    waveform = plot(bowed_sol.t[window] ./ (2pi*bowed_f0), bowed_sol[3,:][window];
        xlabel="time (s)", ylabel="resonator q", legend=false)
    orbit = plot(bowed_sol[1,:][window], bowed_sol[2,:][window];
        xlabel="source x", ylabel="source v", legend=false)
    plot(waveform, orbit; layout=(1,2), size=(800,300))
end

# ╔═╡ e29097e9-5ac9-5661-9dd4-c6e5959f8e00
bowed_raw_audio = bowed_sol[3,:];

# ╔═╡ af18bf0a-2688-52c5-80f5-44a72189d9e1
let block = round(Int, 0.02*bowed_sample_rate*bowed_oversample)
    n = length(bowed_raw_audio) ÷ block
    envelope = [sqrt(sum(abs2, view(bowed_raw_audio, (k-1)*block+1:k*block))/block) for k in 1:n]
    times = ((0:n-1) .+ 0.5) .* (block/(bowed_sample_rate*bowed_oversample))
    plot(times, envelope; xlabel="time (s)", ylabel="resonator RMS", legend=false)
end

# ╔═╡ 4edcd655-f76d-5227-a3cc-3c61400e84e9
bowed_audio = prepare_audio(bowed_raw_audio, bowed_sample_rate; oversample=bowed_oversample);

# ╔═╡ 01198787-5459-5d22-b628-6dc25d45a031
bowed_wav_data = let buffer = IOBuffer()
    wavwrite(bowed_audio, buffer; Fs=bowed_sample_rate, nbits=16, compression=WAVE_FORMAT_PCM)
    take!(buffer)
end;

# ╔═╡ fe02f631-31e4-5729-ad98-1af32a276ab9
md"""
Listen to the resonator coordinate q, with ``τ=2πf_0 t`` setting the pitch scale.
The 4× sampled solution is low-pass filtered and reduced to 48 kHz, then
peak-normalized with short fades; compare the plotted amplitudes for loudness.
"""

# ╔═╡ 262aa357-60f2-5d61-a02e-c0ad15801f23
Resource("data:audio/wav;base64," * base64encode(bowed_wav_data), MIME("audio/wav"), ())

# ╔═╡ 80b71a3f-3176-5583-932e-0e8746cd03c4
DownloadButton(bowed_wav_data, "bowed_resonator.wav")

# ╔═╡ 63591ce3-f86e-51da-a7e7-d8bc53cd4eb2
md"""
### Connect friction to self-oscillation

**Predict → listen → explain:** locate the part of the friction curve that
supplies energy, then compare it with the bow and oscillator velocities.
Vary bow speed and coupling separately. Does the strongest spectral peak
follow the source, the body resonance, or a harmonic? The smooth friction
law approximates slipping; exact sticking would require an additional rule.
"""

# ╔═╡ 21b7f344-c2de-56d3-846b-d71e79bab6d2
md"""
## Struck nonlinear oscillator and resonator

An initial velocity represents a strike; cubic stiffness shifts the source frequency as energy decays into a second resonant mode.
"""

# ╔═╡ ebc1fc29-d6a0-5fc5-b650-70da96b4e2a4
md"""
```math
\begin{aligned}
x' &= v, & v' &= -x-\alpha x^3-2\delta v+\kappa(q-x),\\
q' &= w, & w' &= -\Omega^2q-2\zeta\Omega w+\kappa(x-q).
\end{aligned}
```

Primes denote derivatives with respect to dimensionless time ``τ=2πf_0 t``.
"""

# ╔═╡ 999074c4-8cd7-50b4-9350-26fdf3debb07
function struck!(du, u, p, t)
    x, v, q, w = u
    α, δ, κ, Ω, ζ = p
    du[1] = v
    du[2] = -x - α*x^3 - 2δ*v + κ*(q-x)
    du[3] = w
    du[4] = -Ω^2*q - 2ζ*Ω*w + κ*(x-q)
    nothing
end

# ╔═╡ 7689abf0-c36a-5b10-9555-2357b45cc71b
struck_f0 = 220.0  # Frequency scale in Hz

# ╔═╡ b9aa321b-0f6f-50c6-bfc9-7c7da506d6fc
struck_duration = 4.0  # Seconds

# ╔═╡ 4c6ab4d3-3b77-5139-90a5-8194501c4e52
begin
    struck_sample_rate = 48_000
    struck_oversample = 4
end

# ╔═╡ 90173f0a-4ed1-5325-8577-4d3f2e01df85
struck_u0 = [0.0, 1.0, 0.0, 0.0]

# ╔═╡ 26e20636-59f9-5a8d-b2c0-bf07ff9d15f3
struck_tspan = (0.0, 2pi*struck_f0*struck_duration)  # Dimensionless time τ

# ╔═╡ 0846c33f-f3d4-56ad-87e6-c09543f61eef
struck_p = [0.8, 0.001, 0.08, 2.4, 0.0015]  # α, δ, κ, Ω, ζ

# ╔═╡ 42609673-51ab-5edb-84b6-af737a7631e9
struck_prob = ODEProblem(struck!, struck_u0, struck_tspan, struck_p)

# ╔═╡ dbe9f04d-3f31-5107-9a62-550a6443b676
struck_sample_times = (0:round(Int, struck_duration*struck_sample_rate*struck_oversample)-1) .* (2pi*struck_f0/(struck_sample_rate*struck_oversample));

# ╔═╡ 884b4126-ffb3-572f-b608-6cb17b8663f7
struck_sol = solve(struck_prob, Tsit5(); abstol=1e-9, reltol=1e-7,
    saveat=struck_sample_times, save_end=false, dense=false);

# ╔═╡ add74f4c-8105-58cd-9863-3cb32eb89c41
let window = 1:min(4801, length(struck_sol.t))
    waveform = plot(struck_sol.t[window] ./ (2pi*struck_f0), struck_sol[3,:][window];
        xlabel="time (s)", ylabel="resonator q", legend=false)
    orbit = plot(struck_sol[1,:][window], struck_sol[2,:][window];
        xlabel="source x", ylabel="source v", legend=false)
    plot(waveform, orbit; layout=(1,2), size=(800,300))
end

# ╔═╡ 5a7e4d61-2cb9-5bce-9633-7ddaf910341d
struck_raw_audio = struck_sol[3,:];

# ╔═╡ 89a42758-5c94-507a-947b-a042bf080c02
let block = round(Int, 0.02*struck_sample_rate*struck_oversample)
    n = length(struck_raw_audio) ÷ block
    envelope = [sqrt(sum(abs2, view(struck_raw_audio, (k-1)*block+1:k*block))/block) for k in 1:n]
    times = ((0:n-1) .+ 0.5) .* (block/(struck_sample_rate*struck_oversample))
    plot(times, envelope; xlabel="time (s)", ylabel="resonator RMS", legend=false)
end

# ╔═╡ 4b4d5767-3616-5768-baae-d5fe93fe5a70
struck_audio = prepare_audio(struck_raw_audio, struck_sample_rate; oversample=struck_oversample);

# ╔═╡ 98ce9178-74c6-594c-a1a0-0aaaf447261f
struck_wav_data = let buffer = IOBuffer()
    wavwrite(struck_audio, buffer; Fs=struck_sample_rate, nbits=16, compression=WAVE_FORMAT_PCM)
    take!(buffer)
end;

# ╔═╡ ab791b74-16a5-554e-9a25-efa7663630ec
md"""
Listen to the resonator coordinate q, with ``τ=2πf_0 t`` setting the pitch scale.
The 4× sampled solution is low-pass filtered and reduced to 48 kHz, then
peak-normalized with short fades; compare the plotted amplitudes for loudness.
"""

# ╔═╡ aafc1191-d96f-5e50-a0ff-e98b7b61c07c
Resource("data:audio/wav;base64," * base64encode(struck_wav_data), MIME("audio/wav"), ())

# ╔═╡ 039fbf21-3e2a-518b-938d-2423ef4eeff0
DownloadButton(struck_wav_data, "struck_resonator.wav")

# ╔═╡ ccff33d4-b252-52ef-9f38-2a3fb6961948
md"""
### A transient is not an attracting cycle

**Predict → listen → explain:** increase the strike while keeping the
parameters fixed. Follow the changing instantaneous oscillation period and
the decaying envelope. Compare this with the sustained reed example.
Two peaks in a transient spectrum can represent two decaying modes; they do
not establish a quasiperiodic attractor.
"""

# ╔═╡ Cell order:
# ╠═23bfa7b2-9b8f-572f-bb20-7901fdd84161
# ╟─dd87b0f0-fdf3-5858-8058-b730a11360ff
# ╟─349c8bfa-8317-5164-9118-f762ee1abcbd
# ╟─5813756a-2aa6-5c5a-99bb-bffd1996ed6e
# ╟─b925cdd8-32c6-5320-ad37-d304bec8de47
# ╠═e8db5332-0778-57eb-a419-89cbb06dbd75
# ╠═efa1b4d5-aeef-5eb6-b5fe-aa43374f2f31
# ╠═e1b94c94-e7cb-5c7d-8512-9ec0410d3084
# ╠═1e981501-2642-562a-b065-cf820f51d0e6
# ╠═11e9129e-2d52-5140-a3ba-924b85bd3dd9
# ╠═bdff993b-83a6-56a9-9d8c-24494f798851
# ╠═9d49ecc1-fb2d-5e2d-8c9f-1c6184e1a0cc
# ╠═b206bd76-4107-5393-8144-2e649c090168
# ╠═2f05e51c-6681-5fb6-a938-b3c132bdc48d
# ╠═b88d1334-d885-5400-8496-e39e69572428
# ╠═32f6fba9-e83e-5ebc-9aab-d63a9ffe6c44
# ╠═3ba443fb-d410-5de3-96ac-3e5a34834947
# ╠═aab39aee-3460-5140-bd34-0a339eeaa92e
# ╠═71505b27-7e34-5f55-933e-73c67189ba9a
# ╠═4acfca4f-b51b-5d85-9e6f-c033dce73e71
# ╟─4dc10528-93ec-5041-b9a4-f52b88aa76af
# ╟─f417718f-736f-56db-bcee-f95ebf0e2b2b
# ╟─3c7ac3bf-e652-5eab-95e0-37c9a9ebd2d6
# ╟─6ebb53f3-e3c5-5a18-8803-74a21a1c5e89
# ╟─6b883fff-6f7e-598f-9b45-c3e69fd54f5f
# ╟─ab523375-2d5c-535f-863c-e5a13cf8d419
# ╠═6aee948d-dde3-5758-8902-2c7cd7e54150
# ╠═16b18521-fabb-5281-9297-0131a250779d
# ╠═b6f85afb-123a-5180-97b5-30d29bb18319
# ╠═b80c51cf-fd17-5770-a8d2-1998fdf8dfb4
# ╠═0a58b2ce-b9e5-5be0-ba18-72f82774de03
# ╠═2ebe1af9-99ab-5021-95f7-9190a6d01914
# ╠═d89d73ab-d021-597a-bfb6-aaf696e82be4
# ╠═ccd6bd9c-77af-5636-9c0c-907fa640b216
# ╠═3e6e6228-9607-575c-9f7b-881cdef11c8f
# ╠═842144ba-df5f-5024-8d57-35824957a2f2
# ╠═1029fda6-0708-571d-83e3-93f6e0f3e2ec
# ╠═e29097e9-5ac9-5661-9dd4-c6e5959f8e00
# ╠═af18bf0a-2688-52c5-80f5-44a72189d9e1
# ╠═4edcd655-f76d-5227-a3cc-3c61400e84e9
# ╠═01198787-5459-5d22-b628-6dc25d45a031
# ╟─fe02f631-31e4-5729-ad98-1af32a276ab9
# ╟─262aa357-60f2-5d61-a02e-c0ad15801f23
# ╟─80b71a3f-3176-5583-932e-0e8746cd03c4
# ╟─63591ce3-f86e-51da-a7e7-d8bc53cd4eb2
# ╟─21b7f344-c2de-56d3-846b-d71e79bab6d2
# ╟─ebc1fc29-d6a0-5fc5-b650-70da96b4e2a4
# ╠═999074c4-8cd7-50b4-9350-26fdf3debb07
# ╠═7689abf0-c36a-5b10-9555-2357b45cc71b
# ╠═b9aa321b-0f6f-50c6-bfc9-7c7da506d6fc
# ╠═4c6ab4d3-3b77-5139-90a5-8194501c4e52
# ╠═90173f0a-4ed1-5325-8577-4d3f2e01df85
# ╠═26e20636-59f9-5a8d-b2c0-bf07ff9d15f3
# ╠═0846c33f-f3d4-56ad-87e6-c09543f61eef
# ╠═42609673-51ab-5edb-84b6-af737a7631e9
# ╠═dbe9f04d-3f31-5107-9a62-550a6443b676
# ╠═884b4126-ffb3-572f-b608-6cb17b8663f7
# ╠═add74f4c-8105-58cd-9863-3cb32eb89c41
# ╠═5a7e4d61-2cb9-5bce-9633-7ddaf910341d
# ╠═89a42758-5c94-507a-947b-a042bf080c02
# ╠═4b4d5767-3616-5768-baae-d5fe93fe5a70
# ╠═98ce9178-74c6-594c-a1a0-0aaaf447261f
# ╟─ab791b74-16a5-554e-9a25-efa7663630ec
# ╟─aafc1191-d96f-5e50-a0ff-e98b7b61c07c
# ╟─039fbf21-3e2a-518b-938d-2423ef4eeff0
# ╟─ccff33d4-b252-52ef-9f38-2a3fb6961948
