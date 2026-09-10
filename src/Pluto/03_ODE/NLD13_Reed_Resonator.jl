### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 23bfa7b2-9b8f-572f-bb20-7901fdd84161
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots, PlutoUI, WAV, Base64
    using NonLinearDynamics: prepare_audio
end

# ╔═╡ 5813756a-2aa6-5c5a-99bb-bffd1996ed6e
md"""
# Reed oscillator coupled to a bore resonance

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
function model!(du, u, p, t)
    x, v, q, w = u
    μ, κ, Ω, ζ = p
    du[1] = v
    du[2] = -x + μ*(1-v^2)*v + κ*(q-x)
    du[3] = w
    du[4] = -Ω^2*q - 2ζ*Ω*w + κ*(x-q)
    nothing
end

# ╔═╡ efa1b4d5-aeef-5eb6-b5fe-aa43374f2f31
f0 = 440.0  # Frequency scale in Hz

# ╔═╡ e1b94c94-e7cb-5c7d-8512-9ec0410d3084
duration = 4.0  # Seconds

# ╔═╡ 1e981501-2642-562a-b065-cf820f51d0e6
begin
    sample_rate = 48_000
    oversample = 4
end

# ╔═╡ 11e9129e-2d52-5140-a3ba-924b85bd3dd9
u0 = [0.05, 0.0, 0.0, 0.0]

# ╔═╡ bdff993b-83a6-56a9-9d8c-24494f798851
tspan = (0.0, 2pi*f0*duration)  # Dimensionless time τ

# ╔═╡ 9d49ecc1-fb2d-5e2d-8c9f-1c6184e1a0cc
p = [0.3, 0.12, 1.02, 0.025]  # μ, κ, Ω, ζ

# ╔═╡ b206bd76-4107-5393-8144-2e649c090168
prob = ODEProblem(model!, u0, tspan, p)

# ╔═╡ 2f05e51c-6681-5fb6-a938-b3c132bdc48d
sample_times = (0:round(Int, duration*sample_rate*oversample)-1) .* (2pi*f0/(sample_rate*oversample));

# ╔═╡ b88d1334-d885-5400-8496-e39e69572428
sol = solve(prob, Tsit5(); abstol=1e-9, reltol=1e-7,
    saveat=sample_times, save_end=false, dense=false);

# ╔═╡ 32f6fba9-e83e-5ebc-9aab-d63a9ffe6c44
let window = max(1, length(sol.t)-4800):length(sol.t)
    waveform = plot(sol.t[window] ./ (2pi*f0), sol[3,:][window];
        xlabel="time (s)", ylabel="resonator q", legend=false)
    orbit = plot(sol[1,:][window], sol[2,:][window];
        xlabel="source x", ylabel="source v", legend=false)
    plot(waveform, orbit; layout=(1,2), size=(800,300))
end

# ╔═╡ 3ba443fb-d410-5de3-96ac-3e5a34834947
raw_audio = sol[3,:];

# ╔═╡ aab39aee-3460-5140-bd34-0a339eeaa92e
let block = round(Int, 0.02*sample_rate*oversample)
    n = length(raw_audio) ÷ block
    envelope = [sqrt(sum(abs2, view(raw_audio, (k-1)*block+1:k*block))/block) for k in 1:n]
    times = ((0:n-1) .+ 0.5) .* (block/(sample_rate*oversample))
    plot(times, envelope; xlabel="time (s)", ylabel="resonator RMS", legend=false)
end

# ╔═╡ 71505b27-7e34-5f55-933e-73c67189ba9a
audio = prepare_audio(raw_audio, sample_rate; oversample);

# ╔═╡ 4acfca4f-b51b-5d85-9e6f-c033dce73e71
wav_data = let buffer = IOBuffer()
    wavwrite(audio, buffer; Fs=sample_rate, nbits=16, compression=WAVE_FORMAT_PCM)
    take!(buffer)
end;

# ╔═╡ 4dc10528-93ec-5041-b9a4-f52b88aa76af
md"""
Listen to the resonator coordinate q, with ``τ=2πf_0 t`` setting the pitch scale.
The 4× sampled solution is low-pass filtered and reduced to 48 kHz, then
peak-normalized with short fades; compare the plotted amplitudes for loudness.
"""

# ╔═╡ f417718f-736f-56db-bcee-f95ebf0e2b2b
Resource("data:audio/wav;base64," * base64encode(wav_data), MIME("audio/wav"), ())

# ╔═╡ 3c7ac3bf-e652-5eab-95e0-37c9a9ebd2d6
DownloadButton(wav_data, "reed_resonator.wav")

# ╔═╡ Cell order:
# ╠═23bfa7b2-9b8f-572f-bb20-7901fdd84161
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
