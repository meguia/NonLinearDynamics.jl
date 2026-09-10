### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 1025b1f4-39b6-5dea-876a-1f62995e7382
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots, ModelingToolkit, PlutoUI, WAV, Base64
    using NonLinearDynamics: prepare_audio
end

# ╔═╡ af86684c-8b17-5257-a29a-718182cf801d
md"""
# Reed oscillator coupled to a bore resonance

A Rayleigh source exchanges energy with one damped bore mode; the resonator coordinate q is a simple sound-output proxy.
"""

# ╔═╡ df88175f-43d1-5029-b203-ec81cf539c55
md"""
```math
\begin{aligned}
x' &= v, & v' &= -x+\mu(1-v^2)v+\kappa(q-x),\\
q' &= w, & w' &= -\Omega^2q-2\zeta\Omega w+\kappa(x-q).
\end{aligned}
```

Primes denote derivatives with respect to dimensionless time ``τ=2πf_0 t``.
"""

# ╔═╡ 80af0cd1-c678-5128-b2f3-20acd5df4e85
# Variables
begin
    @independent_variables t
    @variables x(t) v(t) q(t) w(t)
end

# ╔═╡ b5a8dbc3-84a5-509c-929c-3a4b34a56605
# Parameters
@parameters μ κ Ω ζ

# ╔═╡ ccae731f-0dd3-5332-bc0f-870eb29215aa
# Time derivative
D = Differential(t)

# ╔═╡ 6ab4ae3a-00a2-5de7-8152-9648727c43a6
equations = [
    D(x) ~ v,
    D(v) ~ -x + μ*(1-v^2)*v + κ*(q-x),
    D(q) ~ w,
    D(w) ~ -Ω^2*q - 2ζ*Ω*w + κ*(x-q)
]

# ╔═╡ b8cbfca6-405b-58c9-902d-a1b16373bbc2
@named system = ODESystem(equations, t)

# ╔═╡ d131a4f6-d636-56fa-be7b-3fe28bba8324
simplified = structural_simplify(system)

# ╔═╡ 7ba8aaf4-e687-506b-ada5-bd164ceb79ab
f0 = 440.0  # Frequency scale in Hz

# ╔═╡ c8d0d87d-82db-5124-ba0f-d44488938635
duration = 4.0  # Seconds

# ╔═╡ f31a6e13-b16b-5544-8b62-9ab0707f4959
begin
    sample_rate = 48_000
    oversample = 4
end

# ╔═╡ 371fc3b8-d9c2-58bf-bd7e-ab0956954e37
u0 = [x => 0.05, v => 0.0, q => 0.0, w => 0.0]

# ╔═╡ 5ce51030-2cfe-513b-81ce-1bb52ba541ea
tspan = (0.0, 2pi*f0*duration)  # Dimensionless time τ

# ╔═╡ d50bc6e9-4dc3-5af3-a06f-75a9bcdefa4f
p = [μ => 0.3, κ => 0.12, Ω => 1.02, ζ => 0.025]

# ╔═╡ 343001f0-5ebd-50d6-b96e-d059bb322497
prob = ODEProblem(simplified, u0, tspan, p)

# ╔═╡ caefc896-92d7-5493-a18d-5bc63cc9e33b
sample_times = (0:round(Int, duration*sample_rate*oversample)-1) .* (2pi*f0/(sample_rate*oversample));

# ╔═╡ 1f7db0f4-d3f1-5b15-b6ed-b4e89eecd5a3
sol = solve(prob, Tsit5(); abstol=1e-9, reltol=1e-7,
    saveat=sample_times, save_end=false, dense=false);

# ╔═╡ d035ffdb-e622-556e-a3f5-340bc71c0e63
let window = max(1, length(sol.t)-4800):length(sol.t)
    waveform = plot(sol.t[window] ./ (2pi*f0), sol[q][window];
        xlabel="time (s)", ylabel="resonator q", legend=false)
    orbit = plot(sol[x][window], sol[v][window];
        xlabel="source x", ylabel="source v", legend=false)
    plot(waveform, orbit; layout=(1,2), size=(800,300))
end

# ╔═╡ 3b924980-1f3c-5abe-9e15-d0c4761afc14
raw_audio = sol[q];

# ╔═╡ 2c1f1e9f-4965-5baa-9e98-43848916eae1
let block = round(Int, 0.02*sample_rate*oversample)
    n = length(raw_audio) ÷ block
    envelope = [sqrt(sum(abs2, view(raw_audio, (k-1)*block+1:k*block))/block) for k in 1:n]
    times = ((0:n-1) .+ 0.5) .* (block/(sample_rate*oversample))
    plot(times, envelope; xlabel="time (s)", ylabel="resonator RMS", legend=false)
end

# ╔═╡ 7d87281d-b489-59cd-8576-61cea5b837c1
audio = prepare_audio(raw_audio, sample_rate; oversample);

# ╔═╡ 2d51f942-58a9-5e6d-8307-8a84bcd4702e
wav_data = let buffer = IOBuffer()
    wavwrite(audio, buffer; Fs=sample_rate, nbits=16, compression=WAVE_FORMAT_PCM)
    take!(buffer)
end;

# ╔═╡ b97d5d76-7cdd-519c-8731-e8a0c419384a
md"""
Listen to the resonator coordinate q, with ``τ=2πf_0 t`` setting the pitch scale.
The 4× sampled solution is low-pass filtered and reduced to 48 kHz, then
peak-normalized with short fades; compare the plotted amplitudes for loudness.
"""

# ╔═╡ 6048a4ae-45e4-50e5-b927-6d239398409e
Resource("data:audio/wav;base64," * base64encode(wav_data), MIME("audio/wav"), ())

# ╔═╡ b1e867f5-9af6-539a-96c4-00f9013a928e
DownloadButton(wav_data, "reed_resonator.wav")

# ╔═╡ Cell order:
# ╠═1025b1f4-39b6-5dea-876a-1f62995e7382
# ╟─af86684c-8b17-5257-a29a-718182cf801d
# ╟─df88175f-43d1-5029-b203-ec81cf539c55
# ╠═80af0cd1-c678-5128-b2f3-20acd5df4e85
# ╠═b5a8dbc3-84a5-509c-929c-3a4b34a56605
# ╠═ccae731f-0dd3-5332-bc0f-870eb29215aa
# ╠═6ab4ae3a-00a2-5de7-8152-9648727c43a6
# ╠═b8cbfca6-405b-58c9-902d-a1b16373bbc2
# ╠═d131a4f6-d636-56fa-be7b-3fe28bba8324
# ╠═7ba8aaf4-e687-506b-ada5-bd164ceb79ab
# ╠═c8d0d87d-82db-5124-ba0f-d44488938635
# ╠═f31a6e13-b16b-5544-8b62-9ab0707f4959
# ╠═371fc3b8-d9c2-58bf-bd7e-ab0956954e37
# ╠═5ce51030-2cfe-513b-81ce-1bb52ba541ea
# ╠═d50bc6e9-4dc3-5af3-a06f-75a9bcdefa4f
# ╠═343001f0-5ebd-50d6-b96e-d059bb322497
# ╠═caefc896-92d7-5493-a18d-5bc63cc9e33b
# ╠═1f7db0f4-d3f1-5b15-b6ed-b4e89eecd5a3
# ╠═d035ffdb-e622-556e-a3f5-340bc71c0e63
# ╠═3b924980-1f3c-5abe-9e15-d0c4761afc14
# ╠═2c1f1e9f-4965-5baa-9e98-43848916eae1
# ╠═7d87281d-b489-59cd-8576-61cea5b837c1
# ╠═2d51f942-58a9-5e6d-8307-8a84bcd4702e
# ╟─b97d5d76-7cdd-519c-8731-e8a0c419384a
# ╟─6048a4ae-45e4-50e5-b927-6d239398409e
# ╟─b1e867f5-9af6-539a-96c4-00f9013a928e
