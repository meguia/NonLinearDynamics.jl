### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 3f8c2b19-ce81-59f6-82c2-63663e0088b6
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots, PlutoUI, WAV, Base64
    using NonLinearDynamics: prepare_audio
end

# ╔═╡ 21b7f344-c2de-56d3-846b-d71e79bab6d2
md"""
# Struck nonlinear oscillator and resonator

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
function model!(du, u, p, t)
    x, v, q, w = u
    α, δ, κ, Ω, ζ = p
    du[1] = v
    du[2] = -x - α*x^3 - 2δ*v + κ*(q-x)
    du[3] = w
    du[4] = -Ω^2*q - 2ζ*Ω*w + κ*(x-q)
    nothing
end

# ╔═╡ 7689abf0-c36a-5b10-9555-2357b45cc71b
f0 = 220.0  # Frequency scale in Hz

# ╔═╡ b9aa321b-0f6f-50c6-bfc9-7c7da506d6fc
duration = 4.0  # Seconds

# ╔═╡ 4c6ab4d3-3b77-5139-90a5-8194501c4e52
begin
    sample_rate = 48_000
    oversample = 4
end

# ╔═╡ 90173f0a-4ed1-5325-8577-4d3f2e01df85
u0 = [0.0, 1.0, 0.0, 0.0]

# ╔═╡ 26e20636-59f9-5a8d-b2c0-bf07ff9d15f3
tspan = (0.0, 2pi*f0*duration)  # Dimensionless time τ

# ╔═╡ 0846c33f-f3d4-56ad-87e6-c09543f61eef
p = [0.8, 0.001, 0.08, 2.4, 0.0015]  # α, δ, κ, Ω, ζ

# ╔═╡ 42609673-51ab-5edb-84b6-af737a7631e9
prob = ODEProblem(model!, u0, tspan, p)

# ╔═╡ dbe9f04d-3f31-5107-9a62-550a6443b676
sample_times = (0:round(Int, duration*sample_rate*oversample)-1) .* (2pi*f0/(sample_rate*oversample));

# ╔═╡ 884b4126-ffb3-572f-b608-6cb17b8663f7
sol = solve(prob, Tsit5(); abstol=1e-9, reltol=1e-7,
    saveat=sample_times, save_end=false, dense=false);

# ╔═╡ add74f4c-8105-58cd-9863-3cb32eb89c41
let window = 1:min(4801, length(sol.t))
    waveform = plot(sol.t[window] ./ (2pi*f0), sol[3,:][window];
        xlabel="time (s)", ylabel="resonator q", legend=false)
    orbit = plot(sol[1,:][window], sol[2,:][window];
        xlabel="source x", ylabel="source v", legend=false)
    plot(waveform, orbit; layout=(1,2), size=(800,300))
end

# ╔═╡ 5a7e4d61-2cb9-5bce-9633-7ddaf910341d
raw_audio = sol[3,:];

# ╔═╡ 89a42758-5c94-507a-947b-a042bf080c02
let block = round(Int, 0.02*sample_rate*oversample)
    n = length(raw_audio) ÷ block
    envelope = [sqrt(sum(abs2, view(raw_audio, (k-1)*block+1:k*block))/block) for k in 1:n]
    times = ((0:n-1) .+ 0.5) .* (block/(sample_rate*oversample))
    plot(times, envelope; xlabel="time (s)", ylabel="resonator RMS", legend=false)
end

# ╔═╡ 4b4d5767-3616-5768-baae-d5fe93fe5a70
audio = prepare_audio(raw_audio, sample_rate; oversample);

# ╔═╡ 98ce9178-74c6-594c-a1a0-0aaaf447261f
wav_data = let buffer = IOBuffer()
    wavwrite(audio, buffer; Fs=sample_rate, nbits=16, compression=WAVE_FORMAT_PCM)
    take!(buffer)
end;

# ╔═╡ ab791b74-16a5-554e-9a25-efa7663630ec
md"""
Listen to the resonator coordinate q, with ``τ=2πf_0 t`` setting the pitch scale.
The 4× sampled solution is low-pass filtered and reduced to 48 kHz, then
peak-normalized with short fades; compare the plotted amplitudes for loudness.
"""

# ╔═╡ aafc1191-d96f-5e50-a0ff-e98b7b61c07c
Resource("data:audio/wav;base64," * base64encode(wav_data), MIME("audio/wav"), ())

# ╔═╡ 039fbf21-3e2a-518b-938d-2423ef4eeff0
DownloadButton(wav_data, "struck_resonator.wav")

# ╔═╡ Cell order:
# ╠═3f8c2b19-ce81-59f6-82c2-63663e0088b6
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
