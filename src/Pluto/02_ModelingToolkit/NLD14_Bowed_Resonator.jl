### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 5cdebbde-2313-5da2-926d-c4c2194943c8
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots, ModelingToolkit, PlutoUI, WAV, Base64
    using NonLinearDynamics: prepare_audio
end

# ╔═╡ 8c983280-04ec-533e-8c4c-cc96da4bd634
md"""
# Bowed oscillator coupled to a body resonance

A smooth sliding-friction law drives a string mode coupled to a damped body mode; q is the sound-output proxy.
"""

# ╔═╡ c1a052cc-39d3-5da4-9167-4bfe18aed702
# Variables
begin
    @independent_variables t
    @variables x(t) v(t) q(t) w(t)
end

# ╔═╡ 81e713b3-941e-5a94-9179-4c2b38cd35b4
# Parameters
@parameters F V ε vs δ κ Ω ζ

# ╔═╡ e8914de8-fb51-5554-a17a-e51f63f23087
# Time derivative
D = Differential(t)

# ╔═╡ 52f5c023-a949-5915-901c-9539a836fc72
equations = [
    D(x) ~ v,
    D(v) ~ -x - δ*v + F*tanh((V-v)/ε)/(1+((V-v)/vs)^2) + κ*(q-x),
    D(q) ~ w,
    D(w) ~ -Ω^2*q - 2ζ*Ω*w + κ*(x-q)
]

# ╔═╡ b1034382-e29a-5e1b-ab83-58654b1583b8
@named system = ODESystem(equations, t)

# ╔═╡ 904febbc-ef15-5419-b0b5-29952ae73934
simplified = structural_simplify(system)

# ╔═╡ c7a30253-040a-5a14-8c76-319879907ef2
f0 = 196.0  # Frequency scale in Hz

# ╔═╡ 492c265f-b9d4-5e92-ab4f-d9048c71e3ce
duration = 4.0  # Seconds

# ╔═╡ ec58fb33-e62f-560f-bbf4-b29d8c707fee
begin
    sample_rate = 48_000
    oversample = 4
end

# ╔═╡ 4f35bb87-cb5f-5e5b-af8a-513fc61aef27
u0 = [x => 0.0, v => 0.0, q => 0.0, w => 0.0]

# ╔═╡ 0f9f8a1e-fbdf-5b7c-b4b3-00df2a207313
tspan = (0.0, 2pi*f0*duration)  # Dimensionless time τ

# ╔═╡ b768eae3-73e5-5e17-bd27-27b11128a6d2
p = [F => 0.7, V => 0.3, ε => 0.02, vs => 0.25, δ => 0.02, κ => 0.08, Ω => 1.4, ζ => 0.02]

# ╔═╡ ab20dcef-6a41-51f2-8a74-e81b04b558b6
prob = ODEProblem(simplified, u0, tspan, p)

# ╔═╡ 2e438808-bac1-5576-a7db-56c5280fad01
sample_times = (0:round(Int, duration*sample_rate*oversample)-1) .* (2pi*f0/(sample_rate*oversample));

# ╔═╡ 8777cf7e-23ab-5f12-a331-b5f4d60ff912
sol = solve(prob, Tsit5(); abstol=1e-9, reltol=1e-7,
    saveat=sample_times, save_end=false, dense=false);

# ╔═╡ fcb149f6-8cc1-555f-a1b5-6af42359dfa0
let window = max(1, length(sol.t)-4800):length(sol.t)
    waveform = plot(sol.t[window] ./ (2pi*f0), sol[q][window];
        xlabel="time (s)", ylabel="resonator q", legend=false)
    orbit = plot(sol[x][window], sol[v][window];
        xlabel="source x", ylabel="source v", legend=false)
    plot(waveform, orbit; layout=(1,2), size=(800,300))
end

# ╔═╡ 0aeaeaf2-6717-5e39-b4eb-9524492aec7c
raw_audio = sol[q];

# ╔═╡ a971175e-1779-52dd-8780-e42dbca83fe1
let block = round(Int, 0.02*sample_rate*oversample)
    n = length(raw_audio) ÷ block
    envelope = [sqrt(sum(abs2, view(raw_audio, (k-1)*block+1:k*block))/block) for k in 1:n]
    times = ((0:n-1) .+ 0.5) .* (block/(sample_rate*oversample))
    plot(times, envelope; xlabel="time (s)", ylabel="resonator RMS", legend=false)
end

# ╔═╡ 9e9949ec-bf0b-5ad5-ac7f-e8d53a712a11
audio = prepare_audio(raw_audio, sample_rate; oversample);

# ╔═╡ 99f76357-283b-51f0-8785-23c73182e958
wav_data = let buffer = IOBuffer()
    wavwrite(audio, buffer; Fs=sample_rate, nbits=16, compression=WAVE_FORMAT_PCM)
    take!(buffer)
end;

# ╔═╡ 58d2ab3c-0283-56aa-8da8-cb199975afca
md"""
Listen to the resonator coordinate q, with ``τ=2πf_0 t`` setting the pitch scale.
The 4× sampled solution is low-pass filtered and reduced to 48 kHz, then
peak-normalized with short fades; compare the plotted amplitudes for loudness.
"""

# ╔═╡ bf5352b0-00c5-5827-b290-5264248f5743
Resource("data:audio/wav;base64," * base64encode(wav_data), MIME("audio/wav"), ())

# ╔═╡ 4b25b9b0-ccdc-51ac-bc29-f485b81b171c
DownloadButton(wav_data, "bowed_resonator.wav")

# ╔═╡ Cell order:
# ╠═5cdebbde-2313-5da2-926d-c4c2194943c8
# ╟─8c983280-04ec-533e-8c4c-cc96da4bd634
# ╠═c1a052cc-39d3-5da4-9167-4bfe18aed702
# ╠═81e713b3-941e-5a94-9179-4c2b38cd35b4
# ╠═e8914de8-fb51-5554-a17a-e51f63f23087
# ╠═52f5c023-a949-5915-901c-9539a836fc72
# ╠═b1034382-e29a-5e1b-ab83-58654b1583b8
# ╠═904febbc-ef15-5419-b0b5-29952ae73934
# ╠═c7a30253-040a-5a14-8c76-319879907ef2
# ╠═492c265f-b9d4-5e92-ab4f-d9048c71e3ce
# ╠═ec58fb33-e62f-560f-bbf4-b29d8c707fee
# ╠═4f35bb87-cb5f-5e5b-af8a-513fc61aef27
# ╠═0f9f8a1e-fbdf-5b7c-b4b3-00df2a207313
# ╠═b768eae3-73e5-5e17-bd27-27b11128a6d2
# ╠═ab20dcef-6a41-51f2-8a74-e81b04b558b6
# ╠═2e438808-bac1-5576-a7db-56c5280fad01
# ╠═8777cf7e-23ab-5f12-a331-b5f4d60ff912
# ╠═fcb149f6-8cc1-555f-a1b5-6af42359dfa0
# ╠═0aeaeaf2-6717-5e39-b4eb-9524492aec7c
# ╠═a971175e-1779-52dd-8780-e42dbca83fe1
# ╠═9e9949ec-bf0b-5ad5-ac7f-e8d53a712a11
# ╠═99f76357-283b-51f0-8785-23c73182e958
# ╟─58d2ab3c-0283-56aa-8da8-cb199975afca
# ╟─bf5352b0-00c5-5827-b290-5264248f5743
# ╟─4b25b9b0-ccdc-51ac-bc29-f485b81b171c
