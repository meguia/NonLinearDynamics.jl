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

# ╔═╡ c0424ad0-c07b-518d-871f-b2e8ff9881b1
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots, PlutoUI
    using RealTimeAudioDiffEq, PlutoHooks
end

# ╔═╡ 4d30f492-b596-5270-8636-0d3b2d123e46
md"""
# Bowed oscillator coupled to a body resonance: live audio

The same ODE drives `DESource`; channel mapping selects the audible state and `set_ts!` maps dimensionless time to seconds.
Choose Play to start, then move the controls; uncheck it to stop. Start with a low speaker volume.
"""

# ╔═╡ c75d8d45-5582-5741-a873-18a8d65aec0d
function model!(du, u, p, t)
    x, v, q, w = u
    F, V, ε, vs, δ, κ, Ω, ζ = p
    du[1] = v
    du[2] = -x - δ*v + F*tanh((V-v)/ε)/(1+((V-v)/vs)^2) + κ*(q-x)
    du[3] = w
    du[4] = -Ω^2*q - 2ζ*Ω*w + κ*(x-q)
    nothing
end

# ╔═╡ 41057cbe-c8bd-54ae-8f46-35521e1ad78b
u0 = [0.0, 0.0, 0.0, 0.0]

# ╔═╡ ab373653-ae42-5693-95e6-9ea0764e5e9f
p = [0.7, 0.3, 0.02, 0.25, 0.02, 0.08, 1.4, 0.02]  # F, V, ε, vs, δ, κ, Ω, ζ

# ╔═╡ cc4a7384-abe9-5ebb-bb33-2700b6417365
source = DESource(model!, copy(u0), copy(p); alg=Tsit5(), channel_map=[3, 3]);

# ╔═╡ d31f6002-3a7f-5d82-9f1a-df8d14513b7a
md"""
V $(@bind live_parameter Slider(0.1:0.05:0.6; default=0.3, show_value=true))

Pitch scale (Hz) $(@bind f0 Slider(110.0:1.0:330.0; default=196.0, show_value=true))

Gain $(@bind gain Slider(0.0:0.01:0.5; default=0.05, show_value=true))
"""

# ╔═╡ 515c4e5c-0425-5099-96df-129392b226a2
controls = begin
    set_param!(source, 2, Float64(live_parameter))
    set_ts!(source, 2pi*f0)
    set_gain!(source, Float64(gain))
    nothing
end

# ╔═╡ 2901b15f-f7d8-5592-af94-4550ae193322
# Compile a short solution before starting the audio callback.
preview = let
    controls
    prob = ODEProblem(model!, copy(u0), (0.0, 40.0), copy(get_params(source)))
    solve(prob, Tsit5(); saveat=0.02, abstol=1e-9, reltol=1e-7)
end;

# ╔═╡ 2879866e-8789-5cf2-8450-4c5bd3bd2cac
plot(preview; idxs=3, xlabel="dimensionless time", ylabel="output", legend=false)

# ╔═╡ 1764a78f-b3cc-5685-9b5e-9ef4c7617c7c
md"""Play $(@bind playing CheckBox(default=false))"""

# ╔═╡ 7781876b-7ded-5adc-8c99-6f48b0d6d33a
# Cleanup runs when Play changes, this cell is rerun, or the notebook is shut down.
@use_effect([source, playing]) do
    if playing
        device = get_default_output_device()
        device < 0 && error("No output device is available. Use the WAV examples or select an audio device with list_devices().")
        reset_state!(source)
        start_DESource(source, device; buffer_size=UInt32(2048))
    end
    return () -> begin
        isactive(source) && stop_DESource(source)
    end
end

# ╔═╡ b74dba24-6fee-59fa-b815-5b13f4869ae4
md"""Uncheck Play to stop. For the struck model, stop and play again to repeat the strike. Live output uses the selected gain; offline WAV examples also filter and normalize the signal."""

# ╔═╡ Cell order:
# ╠═c0424ad0-c07b-518d-871f-b2e8ff9881b1
# ╟─4d30f492-b596-5270-8636-0d3b2d123e46
# ╠═c75d8d45-5582-5741-a873-18a8d65aec0d
# ╠═41057cbe-c8bd-54ae-8f46-35521e1ad78b
# ╠═ab373653-ae42-5693-95e6-9ea0764e5e9f
# ╠═cc4a7384-abe9-5ebb-bb33-2700b6417365
# ╟─d31f6002-3a7f-5d82-9f1a-df8d14513b7a
# ╠═515c4e5c-0425-5099-96df-129392b226a2
# ╠═2901b15f-f7d8-5592-af94-4550ae193322
# ╠═2879866e-8789-5cf2-8450-4c5bd3bd2cac
# ╟─1764a78f-b3cc-5685-9b5e-9ef4c7617c7c
# ╠═7781876b-7ded-5adc-8c99-6f48b0d6d33a
# ╟─b74dba24-6fee-59fa-b815-5b13f4869ae4
