### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ a7d22e91-e01c-5de1-869c-ffeeb3a7d421
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots, PlutoUI, WAV, Base64
    using NonLinearDynamics: prepare_audio
end

# ╔═╡ 6b883fff-6f7e-598f-9b45-c3e69fd54f5f
md"""
# Bowed oscillator coupled to a body resonance

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
function model!(du, u, p, t)
    x, v, q, w = u
    F, V, ε, vs, δ, κ, Ω, ζ = p
    du[1] = v
    du[2] = -x - δ*v + F*tanh((V-v)/ε)/(1+((V-v)/vs)^2) + κ*(q-x)
    du[3] = w
    du[4] = -Ω^2*q - 2ζ*Ω*w + κ*(x-q)
    nothing
end

# ╔═╡ 16b18521-fabb-5281-9297-0131a250779d
f0 = 196.0  # Frequency scale in Hz

# ╔═╡ b6f85afb-123a-5180-97b5-30d29bb18319
duration = 4.0  # Seconds

# ╔═╡ b80c51cf-fd17-5770-a8d2-1998fdf8dfb4
begin
    sample_rate = 48_000
    oversample = 4
end

# ╔═╡ 0a58b2ce-b9e5-5be0-ba18-72f82774de03
u0 = [0.0, 0.0, 0.0, 0.0]

# ╔═╡ 2ebe1af9-99ab-5021-95f7-9190a6d01914
tspan = (0.0, 2pi*f0*duration)  # Dimensionless time τ

# ╔═╡ d89d73ab-d021-597a-bfb6-aaf696e82be4
p = [0.7, 0.3, 0.02, 0.25, 0.02, 0.08, 1.4, 0.02]  # F, V, ε, vs, δ, κ, Ω, ζ

# ╔═╡ ccd6bd9c-77af-5636-9c0c-907fa640b216
prob = ODEProblem(model!, u0, tspan, p)

# ╔═╡ 3e6e6228-9607-575c-9f7b-881cdef11c8f
sample_times = (0:round(Int, duration*sample_rate*oversample)-1) .* (2pi*f0/(sample_rate*oversample));

# ╔═╡ 842144ba-df5f-5024-8d57-35824957a2f2
sol = solve(prob, Tsit5(); abstol=1e-9, reltol=1e-7,
    saveat=sample_times, save_end=false, dense=false);

# ╔═╡ 1029fda6-0708-571d-83e3-93f6e0f3e2ec
let window = max(1, length(sol.t)-4800):length(sol.t)
    waveform = plot(sol.t[window] ./ (2pi*f0), sol[3,:][window];
        xlabel="time (s)", ylabel="resonator q", legend=false)
    orbit = plot(sol[1,:][window], sol[2,:][window];
        xlabel="source x", ylabel="source v", legend=false)
    plot(waveform, orbit; layout=(1,2), size=(800,300))
end

# ╔═╡ e29097e9-5ac9-5661-9dd4-c6e5959f8e00
raw_audio = sol[3,:];

# ╔═╡ af18bf0a-2688-52c5-80f5-44a72189d9e1
let block = round(Int, 0.02*sample_rate*oversample)
    n = length(raw_audio) ÷ block
    envelope = [sqrt(sum(abs2, view(raw_audio, (k-1)*block+1:k*block))/block) for k in 1:n]
    times = ((0:n-1) .+ 0.5) .* (block/(sample_rate*oversample))
    plot(times, envelope; xlabel="time (s)", ylabel="resonator RMS", legend=false)
end

# ╔═╡ 4edcd655-f76d-5227-a3cc-3c61400e84e9
audio = prepare_audio(raw_audio, sample_rate; oversample);

# ╔═╡ 01198787-5459-5d22-b628-6dc25d45a031
wav_data = let buffer = IOBuffer()
    wavwrite(audio, buffer; Fs=sample_rate, nbits=16, compression=WAVE_FORMAT_PCM)
    take!(buffer)
end;

# ╔═╡ fe02f631-31e4-5729-ad98-1af32a276ab9
md"""
Listen to the resonator coordinate q, with ``τ=2πf_0 t`` setting the pitch scale.
The 4× sampled solution is low-pass filtered and reduced to 48 kHz, then
peak-normalized with short fades; compare the plotted amplitudes for loudness.
"""

# ╔═╡ 262aa357-60f2-5d61-a02e-c0ad15801f23
Resource("data:audio/wav;base64," * base64encode(wav_data), MIME("audio/wav"), ())

# ╔═╡ 80b71a3f-3176-5583-932e-0e8746cd03c4
DownloadButton(wav_data, "bowed_resonator.wav")

# ╔═╡ Cell order:
# ╠═a7d22e91-e01c-5de1-869c-ffeeb3a7d421
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
