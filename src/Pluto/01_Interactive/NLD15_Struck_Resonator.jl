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

# ╔═╡ ccd639c9-8881-52fa-9703-f40c2dd0a2a5
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots, PlutoUI, WAV, Base64
    using NonLinearDynamics: prepare_audio
end

# ╔═╡ fe002ece-6262-5b8a-8b6a-373073f06d46
md"""
# Struck nonlinear oscillator and resonator

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
function model!(du, u, p, t)
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
α $(@bind drive Slider(0.0:0.1:1.5; default=0.8, show_value=true))

Coupling κ $(@bind coupling Slider(0.0:0.02:0.3; default=0.08, show_value=true))

Resonator ratio Ω $(@bind ratio Slider(0.5:0.02:3.0; default=2.4, show_value=true))

Pitch scale (Hz) $(@bind f0 Slider(110.0:1.0:330.0; default=220.0, show_value=true))

Increase α to hear the pitch change during the decay, and vary Ω to compare harmonic and inharmonic resonances.
"""

# ╔═╡ 789aea5e-9156-5142-98b9-1a45f8867ab8
begin
    duration = 4.0
    sample_rate = 48_000
    oversample = 4
end

# ╔═╡ b6d52cd7-e312-5807-ad98-8f9eed2e28ed
u0 = [0.0, 1.0, 0.0, 0.0]

# ╔═╡ 1a49cdc0-7e14-5702-8292-20d60f4d688c
tspan = (0.0, 2pi*f0*duration)  # Dimensionless time τ

# ╔═╡ 6aafba6c-9100-54ee-8f5e-38f7bc1e6b0d
p = [drive, 0.001, coupling, ratio, 0.0015]  # α, δ, κ, Ω, ζ

# ╔═╡ ddc9b0d5-eee8-541e-93c0-a47bb23863b4
prob = ODEProblem(model!, u0, tspan, p)

# ╔═╡ e075e7b8-4007-5e83-ab16-9ad0b7c90492
sample_times = (0:round(Int, duration*sample_rate*oversample)-1) .* (2pi*f0/(sample_rate*oversample));

# ╔═╡ ebc4d026-0b43-5c3b-8c99-46414649d968
sol = solve(prob, Tsit5(); abstol=1e-9, reltol=1e-7,
    saveat=sample_times, save_end=false, dense=false);

# ╔═╡ 72c35fb1-46bb-566b-9494-713c1be745d6
let window = 1:min(4801, length(sol.t))
    waveform = plot(sol.t[window] ./ (2pi*f0), sol[3,:][window];
        xlabel="time (s)", ylabel="resonator q", legend=false)
    orbit = plot(sol[1,:][window], sol[2,:][window];
        xlabel="source x", ylabel="source v", legend=false)
    plot(waveform, orbit; layout=(1,2), size=(800,300))
end

# ╔═╡ 76424dc2-4374-586e-8ef8-e286bcb24571
raw_audio = sol[3,:];

# ╔═╡ 5d0c1755-385d-5bc4-8bdb-f61976ad9635
let block = round(Int, 0.02*sample_rate*oversample)
    n = length(raw_audio) ÷ block
    envelope = [sqrt(sum(abs2, view(raw_audio, (k-1)*block+1:k*block))/block) for k in 1:n]
    times = ((0:n-1) .+ 0.5) .* (block/(sample_rate*oversample))
    plot(times, envelope; xlabel="time (s)", ylabel="resonator RMS", legend=false)
end

# ╔═╡ 62783ade-795f-5cb8-9af1-245c6e7683b7
audio = prepare_audio(raw_audio, sample_rate; oversample);

# ╔═╡ 4234c63c-b5df-5db7-95d3-609800eea854
wav_data = let buffer = IOBuffer()
    wavwrite(audio, buffer; Fs=sample_rate, nbits=16, compression=WAVE_FORMAT_PCM)
    take!(buffer)
end;

# ╔═╡ d5b1bb31-34dd-592a-a86a-413cbe172d45
md"""
Listen to the resonator coordinate q, with ``τ=2πf_0 t`` setting the pitch scale.
The 4× sampled solution is low-pass filtered and reduced to 48 kHz, then
peak-normalized with short fades; compare the plotted amplitudes for loudness.
"""

# ╔═╡ 87cb660a-a189-580f-b180-a484e0ba9376
Resource("data:audio/wav;base64," * base64encode(wav_data), MIME("audio/wav"), ())

# ╔═╡ 91718610-02a8-55a5-af4f-d70671a92272
DownloadButton(wav_data, "struck_resonator.wav")

# ╔═╡ Cell order:
# ╠═ccd639c9-8881-52fa-9703-f40c2dd0a2a5
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
# ╟─72c35fb1-46bb-566b-9494-713c1be745d6
# ╟─76424dc2-4374-586e-8ef8-e286bcb24571
# ╟─5d0c1755-385d-5bc4-8bdb-f61976ad9635
# ╟─62783ade-795f-5cb8-9af1-245c6e7683b7
# ╟─4234c63c-b5df-5db7-95d3-609800eea854
# ╟─d5b1bb31-34dd-592a-a86a-413cbe172d45
# ╟─87cb660a-a189-580f-b180-a484e0ba9376
# ╟─91718610-02a8-55a5-af4f-d70671a92272
