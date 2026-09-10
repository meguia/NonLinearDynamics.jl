### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ ab24b342-aca8-5c16-88b0-74f8685f8f31
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots, PlutoUI, DSP, WAV, Base64
    using NonLinearDynamics: prepare_audio
end

# ╔═╡ ffae312b-8218-51c2-86d1-d6c91d6c3176
md"""
# Nonlinear reed oscillator coupled to two resonators

Two acoustic modes interact through the reed and its nonlinear airflow, allowing a periodic tone or a two-frequency oscillation.
"""

# ╔═╡ 39cd05ad-b46f-5a91-a949-11efa6a8a40e
md"""
```math
\begin{aligned}
P &= p_1+p_2, & \Delta &= \gamma-P, & U &= \zeta\max(h,0)\frac{\Delta}{(\Delta^2+\epsilon^2)^{1/4}},\\
h' &= v, & v' &= -d_r\Omega_r v-\Omega_r^2(h-1+\Delta),\\
q_1' &= p_1, & p_1' &= -q_1-\frac{p_1}{Q_1}+F_1U,\\
q_2' &= p_2, & p_2' &= -\Omega_2^2q_2-\frac{\Omega_2}{Q_2}p_2+F_2U.
\end{aligned}
```

h is the reed opening, p₁ and p₂ are acoustic pressures, and q₁ and q₂ are their time integrals.
Constant blowing pressure γ supplies energy through the nonlinear airflow U; both resonators feed pressure back to the reed.
Primes use ``τ=2πf_0t`` and the first resonance has ``Ω_1=1``.
The small ε smooths the square-root slope at zero pressure drop. Closure stops flow; reed-contact forces are omitted.
"""

# ╔═╡ 07b1d605-8fdf-5a73-ba21-e9c9d2e4ce8d
function model!(du, u, p, t)
    h, v, q1, p1, q2, p2 = u
    γ, ζ, Ωr, dr, Ω2, Q1, Q2, F1, F2, ε = p
    Δ = γ - (p1+p2)
    flow = ζ*max(h, 0)*Δ/sqrt(sqrt(Δ^2+ε^2))
    du[1] = v
    du[2] = -dr*Ωr*v - Ωr^2*(h-1+Δ)
    du[3] = p1
    du[4] = -q1 - p1/Q1 + F1*flow
    du[5] = p2
    du[6] = -Ω2^2*q2 - Ω2*p2/Q2 + F2*flow
    nothing
end

# ╔═╡ b22bb3fc-e3a8-57af-9914-bb2b3c7aaa7a
f0 = 185.0  # Frequency scale in Hz

# ╔═╡ bda8044f-1655-5a30-946e-5865540d9c76
# Static equilibrium at the default parameters, with p₁ perturbed by 0.01.
u0 = [0.53, 0.0, 0.11582631989741293, 0.01, 0.03344769734042085, 0.0]

# ╔═╡ 59957770-3d6e-5fb2-8920-4457e68ac2b2
begin
    duration = 4.0
    sample_rate = 48_000
    oversample = 4
    transient = 5000.0  # Dimensionless time discarded before listening and analysis
end

# ╔═╡ fdf42222-3728-584c-b37c-9cf91eb374be
tspan = (0.0, transient + 2pi*f0*duration)

# ╔═╡ a3f65395-3d69-5f6d-a56e-8fb3d048342e
p = [0.47, 0.28, 2pi*2000/1161.2, 1.0, 2.5, 36.6, 41.2, 1322/1161.2, 2386/1161.2, 1e-4]
# γ, ζ, Ωr, dr, Ω2, Q1, Q2, F1, F2, ε

# ╔═╡ 2fbf8251-3a07-5ae1-a127-3eae80a398f7
prob = ODEProblem(model!, u0, tspan, p)

# ╔═╡ 9b2248b4-dcf7-58b3-ba31-a0970fd349c7
sample_times = transient .+ (0:round(Int, duration*sample_rate*oversample)-1) .* (2pi*f0/(sample_rate*oversample));

# ╔═╡ b7f6f4a9-2d16-5f73-abad-50e1112d20ad
sol = solve(prob, Tsit5(); abstol=1e-10, reltol=1e-8, maxiters=10^7,
    saveat=sample_times, save_start=false, save_end=false, dense=false);

# ╔═╡ c2630331-89ab-5912-97db-90d869aca13d
raw_audio = sol[4,:] .+ sol[6,:];

# ╔═╡ 3df9dbec-b83d-5842-afd4-e9d2279ed03b
md"""
After discarding τ < 5000, a periodic orbit gives isolated Poincaré points; a two-frequency torus gives a closed curve.
Compare the section with the spectrum: two peaks alone do not establish quasiperiodicity, and the oscillation frequencies need not equal the bare resonances.
"""

# ╔═╡ 9175ca5a-b083-5adc-a4b5-4a106ecf79e2
section = let
    keep = findall(t -> t >= 5000.0, sol.t)
    a, b, c = sol[3,:][keep], sol[5,:][keep], sol[6,:][keep]
    times = sol.t[keep]
    level = sum(a)/length(a)
    crossings = [i for i in 1:length(a)-1 if a[i] <= level < a[i+1]]
    fraction = [(level-a[i])/(a[i+1]-a[i]) for i in crossings]
    (q2 = [b[i]+s*(b[i+1]-b[i]) for (i,s) in zip(crossings,fraction)],
     p2 = [c[i]+s*(c[i+1]-c[i]) for (i,s) in zip(crossings,fraction)],
     times = [times[i]+s*(times[i+1]-times[i]) for (i,s) in zip(crossings,fraction)],
     level = level)
end;

# ╔═╡ 5e9bc9ec-90f0-53e1-a78c-12cf1fd8390e
let
    window = findall(t -> t >= last(sol.t)-80.0, sol.t)
    waveform = plot((sol.t[window] .- first(sol.t[window])) ./ (2pi*f0), raw_audio[window];
        xlabel="time (s)", ylabel="pressure p₁ + p₂", legend=false)
    section_plot = scatter(section.q2, section.p2; markersize=2, markerstrokewidth=0,
        xlabel="q₂", ylabel="p₂", title="Poincaré: q₁ = mean(q₁), p₁ > 0", legend=false)
    plot(waveform, section_plot; layout=(1,2), size=(850,320))
end

# ╔═╡ 6ce8dfeb-d603-50ae-aaa1-55564768c590
let
    spectrum = periodogram(raw_audio; fs=sample_rate*oversample, window=hanning)
    relative_power = spectrum.power ./ max(maximum(spectrum.power), eps(Float64))
    plot(spectrum.freq, 10 .* log10.(max.(relative_power, 1e-10));
        xlims=(0, 6f0), ylims=(-80, 5), xlabel="frequency (Hz)", ylabel="relative power (dB)", legend=false)
end

# ╔═╡ d0b1eb6d-1983-5742-8aa5-111bed9f8547
audio = prepare_audio(raw_audio, sample_rate; oversample);

# ╔═╡ 470a23bc-2316-5458-bd1f-2df0cdcd9265
wav_data = let buffer = IOBuffer()
    wavwrite(audio, buffer; Fs=sample_rate, nbits=16, compression=WAVE_FORMAT_PCM)
    take!(buffer)
end;

# ╔═╡ 0903b276-1fbc-52b2-8a8c-694aa5a6d73a
md"""
Listen to the sum of the two modal pressures after the transient.
The 4× sampled solution is low-pass filtered to 48 kHz, peak-normalized, and faded at both ends.
"""

# ╔═╡ e82da74d-d999-5148-baad-8a4c8f2f99b0
Resource("data:audio/wav;base64," * base64encode(wav_data), MIME("audio/wav"), ())

# ╔═╡ 99064c2b-6a05-5901-9e1b-6b9648c1b8bd
DownloadButton(wav_data, "two_resonators.wav")

# ╔═╡ b6449272-6802-58d3-8ede-a47d01379fb0
md"""
The bore parameters are based on [Doc, Vergez and Missoum (2014)](https://doi.org/10.3813/AAA.918734), with detuned resonances and stronger reed damping for this teaching example.
"""

# ╔═╡ Cell order:
# ╠═ab24b342-aca8-5c16-88b0-74f8685f8f31
# ╟─ffae312b-8218-51c2-86d1-d6c91d6c3176
# ╟─39cd05ad-b46f-5a91-a949-11efa6a8a40e
# ╠═07b1d605-8fdf-5a73-ba21-e9c9d2e4ce8d
# ╠═b22bb3fc-e3a8-57af-9914-bb2b3c7aaa7a
# ╠═bda8044f-1655-5a30-946e-5865540d9c76
# ╠═59957770-3d6e-5fb2-8920-4457e68ac2b2
# ╠═fdf42222-3728-584c-b37c-9cf91eb374be
# ╠═a3f65395-3d69-5f6d-a56e-8fb3d048342e
# ╠═2fbf8251-3a07-5ae1-a127-3eae80a398f7
# ╠═9b2248b4-dcf7-58b3-ba31-a0970fd349c7
# ╠═b7f6f4a9-2d16-5f73-abad-50e1112d20ad
# ╠═c2630331-89ab-5912-97db-90d869aca13d
# ╟─3df9dbec-b83d-5842-afd4-e9d2279ed03b
# ╠═9175ca5a-b083-5adc-a4b5-4a106ecf79e2
# ╠═5e9bc9ec-90f0-53e1-a78c-12cf1fd8390e
# ╠═6ce8dfeb-d603-50ae-aaa1-55564768c590
# ╠═d0b1eb6d-1983-5742-8aa5-111bed9f8547
# ╠═470a23bc-2316-5458-bd1f-2df0cdcd9265
# ╟─0903b276-1fbc-52b2-8a8c-694aa5a6d73a
# ╟─e82da74d-d999-5148-baad-8a4c8f2f99b0
# ╟─99064c2b-6a05-5901-9e1b-6b9648c1b8bd
# ╟─b6449272-6802-58d3-8ede-a47d01379fb0
