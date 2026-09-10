### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 23c886fb-7040-5dd4-8e40-300ae975f447
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots, PlutoUI, DSP, WAV, Base64, ModelingToolkit
    using NonLinearDynamics: prepare_audio
end

# ╔═╡ d942539b-14cd-5410-858e-7aa3e7f9d5ca
md"""
# Nonlinear reed oscillator coupled to two resonators

Two acoustic modes interact through the reed and its nonlinear airflow, allowing a periodic tone or a two-frequency oscillation.

MTK turns the symbolic variables, parameters, and equations into a numerical ODE problem.
"""

# ╔═╡ f7d292a0-31ce-565f-b9e0-5fd8ff2eb534
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

# ╔═╡ ac29f334-d6d1-522e-8312-2f40e9da6b88
# Variables
begin
    @independent_variables t
    @variables h(t) v(t) q1(t) p1(t) q2(t) p2(t)
end

# ╔═╡ 932db6f4-ed2d-53c3-b20d-1d0c2ff68f96
# Parameters
@parameters γ ζ Ωr dr Ω2 Q1 Q2 F1 F2 ε

# ╔═╡ 810da3fe-be65-54bc-a655-a1546dfd8e74
# Time derivative
D = Differential(t)

# ╔═╡ f06e1ce2-9305-54ca-9a36-90316068cfbd
equations = let
    Δ = γ - (p1+p2)
    flow = ζ*max(h, 0)*Δ/sqrt(sqrt(Δ^2+ε^2))
    [D(h) ~ v,
     D(v) ~ -dr*Ωr*v - Ωr^2*(h-1+Δ),
     D(q1) ~ p1,
     D(p1) ~ -q1 - p1/Q1 + F1*flow,
     D(q2) ~ p2,
     D(p2) ~ -Ω2^2*q2 - Ω2*p2/Q2 + F2*flow]
end

# ╔═╡ 991f67dd-7f7b-550c-b520-dbfd84c7ef9c
@named system = ODESystem(equations, t)

# ╔═╡ a04941ba-94f9-573b-bc4f-cd39e3d2f332
simplified = structural_simplify(system)

# ╔═╡ 36c15685-d666-5620-9bf0-0dbffdb0d24d
f0 = 185.0  # Frequency scale in Hz

# ╔═╡ a91d42be-0271-5065-8fdf-d7ba7d2b3b8d
# Static equilibrium at the default parameters, with p₁ perturbed by 0.01.
u0 = [h => 0.53, v => 0.0, q1 => 0.11582631989741293, p1 => 0.01, q2 => 0.03344769734042085, p2 => 0.0]

# ╔═╡ eaaa523f-8b16-5181-9319-228593255c1a
begin
    duration = 4.0
    sample_rate = 48_000
    oversample = 4
    transient = 5000.0  # Dimensionless time discarded before listening and analysis
end

# ╔═╡ ab24db88-6556-5664-8f3d-f6d3fbb6ef56
tspan = (0.0, transient + 2pi*f0*duration)

# ╔═╡ 525d089c-e895-5ca8-b338-129ab2152600
p = [γ => 0.47, ζ => 0.28, Ωr => 2pi*2000/1161.2, dr => 1.0, Ω2 => 2.5, Q1 => 36.6, Q2 => 41.2, F1 => 1322/1161.2, F2 => 2386/1161.2, ε => 1e-4]

# ╔═╡ f5d764aa-bf1a-5723-b784-80332fbc3b5f
prob = ODEProblem(simplified, u0, tspan, p)

# ╔═╡ c946d707-2f47-500e-8e03-8886feece0b6
sample_times = transient .+ (0:round(Int, duration*sample_rate*oversample)-1) .* (2pi*f0/(sample_rate*oversample));

# ╔═╡ 6c483e95-b6bd-5775-bf28-c702e7460f72
sol = solve(prob, Tsit5(); abstol=1e-10, reltol=1e-8, maxiters=10^7,
    saveat=sample_times, save_start=false, save_end=false, dense=false);

# ╔═╡ 156f9b1c-d46c-5293-a5a4-dc61881a2e4b
raw_audio = sol[p1] .+ sol[p2];

# ╔═╡ e175cbc8-fe1f-5658-983a-01d763fa4398
md"""
After discarding τ < 5000, a periodic orbit gives isolated Poincaré points; a two-frequency torus gives a closed curve.
Compare the section with the spectrum: two peaks alone do not establish quasiperiodicity, and the oscillation frequencies need not equal the bare resonances.
"""

# ╔═╡ 0ca26082-c1a2-534f-9225-977813b6918c
section = let
    keep = findall(t -> t >= 5000.0, sol.t)
    a, b, c = sol[q1][keep], sol[q2][keep], sol[p2][keep]
    times = sol.t[keep]
    level = sum(a)/length(a)
    crossings = [i for i in 1:length(a)-1 if a[i] <= level < a[i+1]]
    fraction = [(level-a[i])/(a[i+1]-a[i]) for i in crossings]
    (q2 = [b[i]+s*(b[i+1]-b[i]) for (i,s) in zip(crossings,fraction)],
     p2 = [c[i]+s*(c[i+1]-c[i]) for (i,s) in zip(crossings,fraction)],
     times = [times[i]+s*(times[i+1]-times[i]) for (i,s) in zip(crossings,fraction)],
     level = level)
end;

# ╔═╡ 107f5cc3-8587-575e-a7b5-99f943802b7b
let
    window = findall(t -> t >= last(sol.t)-80.0, sol.t)
    waveform = plot((sol.t[window] .- first(sol.t[window])) ./ (2pi*f0), raw_audio[window];
        xlabel="time (s)", ylabel="pressure p₁ + p₂", legend=false)
    section_plot = scatter(section.q2, section.p2; markersize=2, markerstrokewidth=0,
        xlabel="q₂", ylabel="p₂", title="Poincaré: q₁ = mean(q₁), p₁ > 0", legend=false)
    plot(waveform, section_plot; layout=(1,2), size=(850,320))
end

# ╔═╡ 750fb5ec-6df2-5148-9448-7346000fae31
let
    spectrum = periodogram(raw_audio; fs=sample_rate*oversample, window=hanning)
    relative_power = spectrum.power ./ max(maximum(spectrum.power), eps(Float64))
    plot(spectrum.freq, 10 .* log10.(max.(relative_power, 1e-10));
        xlims=(0, 6f0), ylims=(-80, 5), xlabel="frequency (Hz)", ylabel="relative power (dB)", legend=false)
end

# ╔═╡ 466a226d-5c57-55ca-b498-113acbc1e79e
audio = prepare_audio(raw_audio, sample_rate; oversample);

# ╔═╡ 5c64f6d1-66bd-508d-9a9a-31266e20ebc9
wav_data = let buffer = IOBuffer()
    wavwrite(audio, buffer; Fs=sample_rate, nbits=16, compression=WAVE_FORMAT_PCM)
    take!(buffer)
end;

# ╔═╡ fec11a0e-26f2-5773-943a-eb1881a08772
md"""
Listen to the sum of the two modal pressures after the transient.
The 4× sampled solution is low-pass filtered to 48 kHz, peak-normalized, and faded at both ends.
"""

# ╔═╡ 250af4ee-fe5c-59bc-9ae9-059baac357a1
Resource("data:audio/wav;base64," * base64encode(wav_data), MIME("audio/wav"), ())

# ╔═╡ 887d4166-44d0-5e6c-987a-1b29d1965c5b
DownloadButton(wav_data, "two_resonators.wav")

# ╔═╡ d448267b-608f-5801-ad67-746a8f4ea335
md"""
The bore parameters are based on [Doc, Vergez and Missoum (2014)](https://doi.org/10.3813/AAA.918734), with detuned resonances and stronger reed damping for this teaching example.
"""

# ╔═╡ Cell order:
# ╠═23c886fb-7040-5dd4-8e40-300ae975f447
# ╟─d942539b-14cd-5410-858e-7aa3e7f9d5ca
# ╟─f7d292a0-31ce-565f-b9e0-5fd8ff2eb534
# ╠═ac29f334-d6d1-522e-8312-2f40e9da6b88
# ╠═932db6f4-ed2d-53c3-b20d-1d0c2ff68f96
# ╠═810da3fe-be65-54bc-a655-a1546dfd8e74
# ╠═f06e1ce2-9305-54ca-9a36-90316068cfbd
# ╠═991f67dd-7f7b-550c-b520-dbfd84c7ef9c
# ╠═a04941ba-94f9-573b-bc4f-cd39e3d2f332
# ╠═36c15685-d666-5620-9bf0-0dbffdb0d24d
# ╠═a91d42be-0271-5065-8fdf-d7ba7d2b3b8d
# ╠═eaaa523f-8b16-5181-9319-228593255c1a
# ╠═ab24db88-6556-5664-8f3d-f6d3fbb6ef56
# ╠═525d089c-e895-5ca8-b338-129ab2152600
# ╠═f5d764aa-bf1a-5723-b784-80332fbc3b5f
# ╠═c946d707-2f47-500e-8e03-8886feece0b6
# ╠═6c483e95-b6bd-5775-bf28-c702e7460f72
# ╠═156f9b1c-d46c-5293-a5a4-dc61881a2e4b
# ╟─e175cbc8-fe1f-5658-983a-01d763fa4398
# ╠═0ca26082-c1a2-534f-9225-977813b6918c
# ╠═107f5cc3-8587-575e-a7b5-99f943802b7b
# ╠═750fb5ec-6df2-5148-9448-7346000fae31
# ╠═466a226d-5c57-55ca-b498-113acbc1e79e
# ╠═5c64f6d1-66bd-508d-9a9a-31266e20ebc9
# ╟─fec11a0e-26f2-5773-943a-eb1881a08772
# ╟─250af4ee-fe5c-59bc-9ae9-059baac357a1
# ╟─887d4166-44d0-5e6c-987a-1b29d1965c5b
# ╟─d448267b-608f-5801-ad67-746a8f4ea335
