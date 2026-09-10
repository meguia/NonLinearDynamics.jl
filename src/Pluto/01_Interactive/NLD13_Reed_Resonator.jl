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

# ╔═╡ 7b801516-e388-5740-8b96-6ac2afe673ec
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots, PlutoUI, WAV, Base64
    using NonLinearDynamics: prepare_audio
end

# ╔═╡ cc6008d5-43fb-5573-9f71-0357022f5f45
md"""
# Reed oscillator coupled to a bore resonance

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
function model!(du, u, p, t)
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
μ $(@bind drive Slider(0.0:0.05:0.8; default=0.3, show_value=true))

Coupling κ $(@bind coupling Slider(0.0:0.02:0.3; default=0.12, show_value=true))

Resonator ratio Ω $(@bind ratio Slider(0.5:0.02:3.0; default=1.02, show_value=true))

Pitch scale (Hz) $(@bind f0 Slider(110.0:1.0:330.0; default=220.0, show_value=true))

Lower μ until oscillations die out, then vary Ω to hear the effect of the resonator.
"""

# ╔═╡ 03083e00-4677-5bbd-a5c8-a9737e680413
begin
    duration = 4.0
    sample_rate = 48_000
    oversample = 4
end

# ╔═╡ 2d3f7ae1-4a97-5f21-b4d1-be9e54b2e728
u0 = [0.05, 0.0, 0.0, 0.0]

# ╔═╡ 904063aa-5ca4-5afa-9783-1cc4f9c5123b
tspan = (0.0, 2pi*f0*duration)  # Dimensionless time τ

# ╔═╡ 0bc5c78a-eb57-5380-85e9-fe6880dfac3f
p = [drive, coupling, ratio, 0.025]  # μ, κ, Ω, ζ

# ╔═╡ ad10e763-0b83-5f29-b7ab-f05a6fdc9b0d
prob = ODEProblem(model!, u0, tspan, p)

# ╔═╡ bbc9bde3-a3c5-58ba-9a1a-2326345ae228
sample_times = (0:round(Int, duration*sample_rate*oversample)-1) .* (2pi*f0/(sample_rate*oversample));

# ╔═╡ abab7691-cabb-52b8-8203-60e256467c52
sol = solve(prob, Tsit5(); abstol=1e-9, reltol=1e-7,
    saveat=sample_times, save_end=false, dense=false);

# ╔═╡ 18540938-de0a-5162-9538-701c475fb13b
let window = max(1, length(sol.t)-4800):length(sol.t)
    waveform = plot(sol.t[window] ./ (2pi*f0), sol[3,:][window];
        xlabel="time (s)", ylabel="resonator q", legend=false)
    orbit = plot(sol[1,:][window], sol[2,:][window];
        xlabel="source x", ylabel="source v", legend=false)
    plot(waveform, orbit; layout=(1,2), size=(800,300))
end

# ╔═╡ 6734c4c3-4e8d-53c4-b72c-442f6e9ea64d
raw_audio = sol[3,:];

# ╔═╡ 23557e4c-8447-5782-98e8-4614cf89248d
let block = round(Int, 0.02*sample_rate*oversample)
    n = length(raw_audio) ÷ block
    envelope = [sqrt(sum(abs2, view(raw_audio, (k-1)*block+1:k*block))/block) for k in 1:n]
    times = ((0:n-1) .+ 0.5) .* (block/(sample_rate*oversample))
    plot(times, envelope; xlabel="time (s)", ylabel="resonator RMS", legend=false)
end

# ╔═╡ 670122bc-c363-551e-b2e3-f1939c9e61d6
audio = prepare_audio(raw_audio, sample_rate; oversample);

# ╔═╡ b14922b2-7c94-545c-916f-5be9562d4295
wav_data = let buffer = IOBuffer()
    wavwrite(audio, buffer; Fs=sample_rate, nbits=16, compression=WAVE_FORMAT_PCM)
    take!(buffer)
end;

# ╔═╡ b324eeb9-ef55-5fc1-b578-e1b88f065a8b
md"""
Listen to the resonator coordinate q, with ``τ=2πf_0 t`` setting the pitch scale.
The 4× sampled solution is low-pass filtered and reduced to 48 kHz, then
peak-normalized with short fades; compare the plotted amplitudes for loudness.
"""

# ╔═╡ a513a223-adf0-5bb8-a30a-897806f43905
Resource("data:audio/wav;base64," * base64encode(wav_data), MIME("audio/wav"), ())

# ╔═╡ df80b8d3-6420-5c65-9822-d55f4e7492bf
DownloadButton(wav_data, "reed_resonator.wav")

# ╔═╡ Cell order:
# ╠═7b801516-e388-5740-8b96-6ac2afe673ec
# ╟─cc6008d5-43fb-5573-9f71-0357022f5f45
# ╠═63eeec30-d79c-5393-a833-0b60ba752742
# ╟─07e9803f-d742-5f20-8517-584ed8216fd1
# ╟─03083e00-4677-5bbd-a5c8-a9737e680413
# ╠═2d3f7ae1-4a97-5f21-b4d1-be9e54b2e728
# ╟─904063aa-5ca4-5afa-9783-1cc4f9c5123b
# ╠═0bc5c78a-eb57-5380-85e9-fe6880dfac3f
# ╟─ad10e763-0b83-5f29-b7ab-f05a6fdc9b0d
# ╟─bbc9bde3-a3c5-58ba-9a1a-2326345ae228
# ╟─abab7691-cabb-52b8-8203-60e256467c52
# ╟─18540938-de0a-5162-9538-701c475fb13b
# ╟─6734c4c3-4e8d-53c4-b72c-442f6e9ea64d
# ╟─23557e4c-8447-5782-98e8-4614cf89248d
# ╟─670122bc-c363-551e-b2e3-f1939c9e61d6
# ╟─b14922b2-7c94-545c-916f-5be9562d4295
# ╟─b324eeb9-ef55-5fc1-b578-e1b88f065a8b
# ╟─a513a223-adf0-5bb8-a30a-897806f43905
# ╟─df80b8d3-6420-5c65-9822-d55f4e7492bf
