### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ b0f0ddb4-dc7a-5a33-abc0-75bbf66926da
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots, ModelingToolkit, PlutoUI, WAV, Base64
    using NonLinearDynamics: prepare_audio
end

# ╔═╡ 82bdd92e-4573-5d9a-93f8-3ceeb546ee45
md"""
# Struck nonlinear oscillator and resonator

An initial velocity represents a strike; cubic stiffness shifts the source frequency as energy decays into a second resonant mode.
"""

# ╔═╡ cf08a6f6-665d-5720-a4fb-2e87a91a547f
md"""
```math
\begin{aligned}
x' &= v, & v' &= -x-\alpha x^3-2\delta v+\kappa(q-x),\\
q' &= w, & w' &= -\Omega^2q-2\zeta\Omega w+\kappa(x-q).
\end{aligned}
```

Primes denote derivatives with respect to dimensionless time ``τ=2πf_0 t``.
"""

# ╔═╡ 576759a3-4d4c-5033-84b7-d3752f69b2ed
# Variables
begin
    @independent_variables t
    @variables x(t) v(t) q(t) w(t)
end

# ╔═╡ c1669c86-73a6-5de2-9f47-3d397ddf18fb
# Parameters
@parameters α δ κ Ω ζ

# ╔═╡ 8da2a20c-086e-5ec5-843f-ac406f933104
# Time derivative
D = Differential(t)

# ╔═╡ 3fd23e85-507a-5f07-9217-6067185bcb8f
equations = [
    D(x) ~ v,
    D(v) ~ -x - α*x^3 - 2δ*v + κ*(q-x),
    D(q) ~ w,
    D(w) ~ -Ω^2*q - 2ζ*Ω*w + κ*(x-q)
]

# ╔═╡ 98d95fca-03ec-5a7e-8491-6e978003b9a4
@named system = ODESystem(equations, t)

# ╔═╡ 56e8ce57-cb80-5ed4-ad67-f217931dc497
simplified = structural_simplify(system)

# ╔═╡ 6006f626-f947-5193-a267-2381953595b7
f0 = 220.0  # Frequency scale in Hz

# ╔═╡ d2080e76-ac67-5e56-9671-befbabfa5aff
duration = 4.0  # Seconds

# ╔═╡ 3a401240-d5cb-5f59-870f-d14ea8800472
begin
    sample_rate = 48_000
    oversample = 4
end

# ╔═╡ efbbb042-46e1-59cd-85ed-8a3ee920dda3
u0 = [x => 0.0, v => 1.0, q => 0.0, w => 0.0]

# ╔═╡ 9eac8025-fe93-559f-8634-6e80382fd794
tspan = (0.0, 2pi*f0*duration)  # Dimensionless time τ

# ╔═╡ ee7ffad2-6364-56f7-bac5-3d1ed2a8f915
p = [α => 0.8, δ => 0.001, κ => 0.08, Ω => 2.4, ζ => 0.0015]

# ╔═╡ 7ca7d78b-ff5a-58f5-b7b5-d5eacf079e48
prob = ODEProblem(simplified, u0, tspan, p)

# ╔═╡ 08eec5e3-f091-5e7c-937e-0c4ab982ce66
sample_times = (0:round(Int, duration*sample_rate*oversample)-1) .* (2pi*f0/(sample_rate*oversample));

# ╔═╡ 817a84e8-49f0-5fe4-835f-99794eedd678
sol = solve(prob, Tsit5(); abstol=1e-9, reltol=1e-7,
    saveat=sample_times, save_end=false, dense=false);

# ╔═╡ 137ddf0c-ab82-594f-be92-d544f080882b
let window = 1:min(4801, length(sol.t))
    waveform = plot(sol.t[window] ./ (2pi*f0), sol[q][window];
        xlabel="time (s)", ylabel="resonator q", legend=false)
    orbit = plot(sol[x][window], sol[v][window];
        xlabel="source x", ylabel="source v", legend=false)
    plot(waveform, orbit; layout=(1,2), size=(800,300))
end

# ╔═╡ 569fca3e-cb18-53ca-92e2-f6e4ea50be5b
raw_audio = sol[q];

# ╔═╡ fe5408a9-553f-5909-be7b-87ea7997410f
let block = round(Int, 0.02*sample_rate*oversample)
    n = length(raw_audio) ÷ block
    envelope = [sqrt(sum(abs2, view(raw_audio, (k-1)*block+1:k*block))/block) for k in 1:n]
    times = ((0:n-1) .+ 0.5) .* (block/(sample_rate*oversample))
    plot(times, envelope; xlabel="time (s)", ylabel="resonator RMS", legend=false)
end

# ╔═╡ b8ceabbf-799d-545f-8c21-f2ee41cc7573
audio = prepare_audio(raw_audio, sample_rate; oversample);

# ╔═╡ 06808419-298a-5cf3-bfa7-65f38514de8a
wav_data = let buffer = IOBuffer()
    wavwrite(audio, buffer; Fs=sample_rate, nbits=16, compression=WAVE_FORMAT_PCM)
    take!(buffer)
end;

# ╔═╡ cd572ae8-b63e-5f2e-b064-d8d2fab4a524
md"""
Listen to the resonator coordinate q, with ``τ=2πf_0 t`` setting the pitch scale.
The 4× sampled solution is low-pass filtered and reduced to 48 kHz, then
peak-normalized with short fades; compare the plotted amplitudes for loudness.
"""

# ╔═╡ eed15966-4b36-55d8-aa5b-af3fc7744dc1
Resource("data:audio/wav;base64," * base64encode(wav_data), MIME("audio/wav"), ())

# ╔═╡ d3187edd-0e6a-5943-bfa9-9b18bf47ec37
DownloadButton(wav_data, "struck_resonator.wav")

# ╔═╡ Cell order:
# ╠═b0f0ddb4-dc7a-5a33-abc0-75bbf66926da
# ╟─82bdd92e-4573-5d9a-93f8-3ceeb546ee45
# ╟─cf08a6f6-665d-5720-a4fb-2e87a91a547f
# ╠═576759a3-4d4c-5033-84b7-d3752f69b2ed
# ╠═c1669c86-73a6-5de2-9f47-3d397ddf18fb
# ╠═8da2a20c-086e-5ec5-843f-ac406f933104
# ╠═3fd23e85-507a-5f07-9217-6067185bcb8f
# ╠═98d95fca-03ec-5a7e-8491-6e978003b9a4
# ╠═56e8ce57-cb80-5ed4-ad67-f217931dc497
# ╠═6006f626-f947-5193-a267-2381953595b7
# ╠═d2080e76-ac67-5e56-9671-befbabfa5aff
# ╠═3a401240-d5cb-5f59-870f-d14ea8800472
# ╠═efbbb042-46e1-59cd-85ed-8a3ee920dda3
# ╠═9eac8025-fe93-559f-8634-6e80382fd794
# ╠═ee7ffad2-6364-56f7-bac5-3d1ed2a8f915
# ╠═7ca7d78b-ff5a-58f5-b7b5-d5eacf079e48
# ╠═08eec5e3-f091-5e7c-937e-0c4ab982ce66
# ╠═817a84e8-49f0-5fe4-835f-99794eedd678
# ╠═137ddf0c-ab82-594f-be92-d544f080882b
# ╠═569fca3e-cb18-53ca-92e2-f6e4ea50be5b
# ╠═fe5408a9-553f-5909-be7b-87ea7997410f
# ╠═b8ceabbf-799d-545f-8c21-f2ee41cc7573
# ╠═06808419-298a-5cf3-bfa7-65f38514de8a
# ╟─cd572ae8-b63e-5f2e-b064-d8d2fab4a524
# ╟─eed15966-4b36-55d8-aa5b-af3fc7744dc1
# ╟─d3187edd-0e6a-5943-bfa9-9b18bf47ec37
