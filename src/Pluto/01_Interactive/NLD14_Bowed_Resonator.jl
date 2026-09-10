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

# ╔═╡ 9926f49d-6949-5427-bf59-f8734ea410f5
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots, PlutoUI, WAV, Base64
    using NonLinearDynamics: prepare_audio
end

# ╔═╡ 5399e2f3-e291-58d3-8cb8-eecab2804042
md"""
# Bowed oscillator coupled to a body resonance

The bow moves at speed V; its force changes sign with the relative speed and weakens during fast sliding.
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
function model!(du, u, p, t)
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
V $(@bind drive Slider(0.1:0.05:0.6; default=0.3, show_value=true))

Coupling κ $(@bind coupling Slider(0.0:0.02:0.3; default=0.08, show_value=true))

Resonator ratio Ω $(@bind ratio Slider(0.5:0.02:3.0; default=1.4, show_value=true))

Pitch scale (Hz) $(@bind f0 Slider(110.0:1.0:330.0; default=196.0, show_value=true))

Change the bow speed V and body frequency Ω, then compare the waveform and its sound.
"""

# ╔═╡ 804fdef5-91f3-544c-8a97-8aa30ea12738
begin
    duration = 4.0
    sample_rate = 48_000
    oversample = 4
end

# ╔═╡ 670572a6-ebb9-5d65-a2f1-750987a631c5
u0 = [0.0, 0.0, 0.0, 0.0]

# ╔═╡ 26f5a210-2596-5046-81f1-f12a0b1285e1
tspan = (0.0, 2pi*f0*duration)  # Dimensionless time τ

# ╔═╡ f1b3328e-c3a5-5ba0-b23b-98b956872f8d
p = [0.7, drive, 0.02, 0.25, 0.02, coupling, ratio, 0.02]  # F, V, ε, vs, δ, κ, Ω, ζ

# ╔═╡ 50e48c6f-5c8a-5e89-b66d-f1f369277c90
prob = ODEProblem(model!, u0, tspan, p)

# ╔═╡ 095309b8-9b13-51e6-a020-3f1db7e9a026
sample_times = (0:round(Int, duration*sample_rate*oversample)-1) .* (2pi*f0/(sample_rate*oversample));

# ╔═╡ 13f02eca-5bef-5d4b-a095-d821edd555bd
sol = solve(prob, Tsit5(); abstol=1e-9, reltol=1e-7,
    saveat=sample_times, save_end=false, dense=false);

# ╔═╡ a7f9fd5a-0ec9-547f-b1d0-8eb441d867d4
let window = max(1, length(sol.t)-4800):length(sol.t)
    waveform = plot(sol.t[window] ./ (2pi*f0), sol[3,:][window];
        xlabel="time (s)", ylabel="resonator q", legend=false)
    orbit = plot(sol[1,:][window], sol[2,:][window];
        xlabel="source x", ylabel="source v", legend=false)
    plot(waveform, orbit; layout=(1,2), size=(800,300))
end

# ╔═╡ 2264a931-173d-5c72-a9ea-15e01f436f05
raw_audio = sol[3,:];

# ╔═╡ d8a00dfe-3923-5fb5-a3ae-e50551ac108a
let block = round(Int, 0.02*sample_rate*oversample)
    n = length(raw_audio) ÷ block
    envelope = [sqrt(sum(abs2, view(raw_audio, (k-1)*block+1:k*block))/block) for k in 1:n]
    times = ((0:n-1) .+ 0.5) .* (block/(sample_rate*oversample))
    plot(times, envelope; xlabel="time (s)", ylabel="resonator RMS", legend=false)
end

# ╔═╡ 4979235f-c765-5336-bcba-066e8cc08d16
audio = prepare_audio(raw_audio, sample_rate; oversample);

# ╔═╡ 0e8d1c6f-81d6-5b39-b36f-98996b645d18
wav_data = let buffer = IOBuffer()
    wavwrite(audio, buffer; Fs=sample_rate, nbits=16, compression=WAVE_FORMAT_PCM)
    take!(buffer)
end;

# ╔═╡ cda8c010-e2a4-5926-aa38-031a87ef9acf
md"""
Listen to the resonator coordinate q, with ``τ=2πf_0 t`` setting the pitch scale.
The 4× sampled solution is low-pass filtered and reduced to 48 kHz, then
peak-normalized with short fades; compare the plotted amplitudes for loudness.
"""

# ╔═╡ 65ed3282-b664-5cb1-9763-76ac03117357
Resource("data:audio/wav;base64," * base64encode(wav_data), MIME("audio/wav"), ())

# ╔═╡ 49ab782f-f01f-5ad2-adcf-e91a8889ff67
DownloadButton(wav_data, "bowed_resonator.wav")

# ╔═╡ Cell order:
# ╠═9926f49d-6949-5427-bf59-f8734ea410f5
# ╟─5399e2f3-e291-58d3-8cb8-eecab2804042
# ╠═8c3e256d-d516-5c81-adc5-61b32112ed55
# ╟─9c1103ad-9477-5da9-8618-a265dab90aa7
# ╟─804fdef5-91f3-544c-8a97-8aa30ea12738
# ╠═670572a6-ebb9-5d65-a2f1-750987a631c5
# ╟─26f5a210-2596-5046-81f1-f12a0b1285e1
# ╠═f1b3328e-c3a5-5ba0-b23b-98b956872f8d
# ╟─50e48c6f-5c8a-5e89-b66d-f1f369277c90
# ╟─095309b8-9b13-51e6-a020-3f1db7e9a026
# ╟─13f02eca-5bef-5d4b-a095-d821edd555bd
# ╟─a7f9fd5a-0ec9-547f-b1d0-8eb441d867d4
# ╟─2264a931-173d-5c72-a9ea-15e01f436f05
# ╟─d8a00dfe-3923-5fb5-a3ae-e50551ac108a
# ╟─4979235f-c765-5336-bcba-066e8cc08d16
# ╟─0e8d1c6f-81d6-5b39-b36f-98996b645d18
# ╟─cda8c010-e2a4-5926-aa38-031a87ef9acf
# ╟─65ed3282-b664-5cb1-9763-76ac03117357
# ╟─49ab782f-f01f-5ad2-adcf-e91a8889ff67
