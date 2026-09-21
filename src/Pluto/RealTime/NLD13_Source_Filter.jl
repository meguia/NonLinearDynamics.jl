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

# ╔═╡ c9467c68-6aaa-59d0-8fde-f3b87b04a879
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots, PlutoUI
    using RealTimeAudioDiffEq, PlutoHooks
end

# ╔═╡ 4c2d0c28-f707-55b7-abf3-9b633eb0fe26
md"""
# Source–filter models

A nonlinear source exchanges energy with a damped resonator. Compare reed,
bowed, and struck excitation: the first two can sustain oscillation, while
the strike supplies its energy through the initial condition. In these
coupled models, the resonator also acts back on the source.

Each subsection has its own controls and Play switch. Stop one before
playing the next to compare the sounds.
"""

# ╔═╡ 3e839d74-f0b5-5cb6-837a-7d6b38bc1bb2
TableOfContents()

# ╔═╡ 3b98f12f-a3ee-58c6-8eb6-2833f7f6d832
md"""
## Reed oscillator coupled to a bore resonance: live audio

The same ODE drives `DESource`; channel mapping selects the audible state and `set_ts!` maps dimensionless time to seconds.
Choose Play to start, then move the controls; uncheck it to stop. Start with a low speaker volume.
"""

# ╔═╡ beaa5691-55c0-5135-b085-a8104614feee
md"""
```math
\begin{aligned}
x' &= v, & v' &= -x+\mu(1-v^2)v+\kappa(q-x),\\
q' &= w, & w' &= -\Omega^2q-2\zeta\Omega w+\kappa(x-q).
\end{aligned}
```

Primes denote derivatives with respect to dimensionless time ``τ=2πf_0 t``.
"""

# ╔═╡ 79241bbc-6ee9-5a8a-a14a-71fc0db879bc
function reed!(du, u, p, t)
    x, v, q, w = u
    μ, κ, Ω, ζ = p
    du[1] = v
    du[2] = -x + μ*(1-v^2)*v + κ*(q-x)
    du[3] = w
    du[4] = -Ω^2*q - 2ζ*Ω*w + κ*(x-q)
    nothing
end

# ╔═╡ 659ce933-3134-5402-a72f-e40c06f4f32d
reed_u0 = [0.05, 0.0, 0.0, 0.0]

# ╔═╡ 63c6d438-8adf-52d7-a6e2-c2e3cb2d5c88
reed_p = [0.3, 0.12, 1.02, 0.025]  # μ, κ, Ω, ζ

# ╔═╡ 7bcdb33c-411e-5bdb-8331-34ad56286122
reed_source = DESource(reed!, copy(reed_u0), copy(reed_p); alg=Tsit5(), channel_map=[3, 3]);

# ╔═╡ 773429e4-9f6e-5761-81c0-82c3a8de4a9c
md"""
μ $(@bind reed_live_parameter Slider(0.0:0.05:0.8; default=0.3, show_value=true))

Pitch scale (Hz) $(@bind reed_f0 Slider(110.0:1.0:330.0; default=220.0, show_value=true))

Gain $(@bind reed_gain Slider(0.0:0.01:0.5; default=0.05, show_value=true))
"""

# ╔═╡ b57d251a-6a70-5867-9eb7-40489a888a0e
reed_controls = begin
    set_param!(reed_source, 1, Float64(reed_live_parameter))
    set_ts!(reed_source, 2pi*reed_f0)
    set_gain!(reed_source, Float64(reed_gain))
    nothing
end

# ╔═╡ 9408ed8f-19c9-5c38-9f8f-56ac6b03daea
# Compile a short solution before starting the audio callback.
reed_preview = let
    reed_controls
    reed_prob = ODEProblem(reed!, copy(reed_u0), (0.0, 40.0), copy(get_params(reed_source)))
    solve(reed_prob, Tsit5(); saveat=0.02, abstol=1e-9, reltol=1e-7)
end;

# ╔═╡ a4002915-43b7-5e55-8539-14cada302fed
plot(reed_preview; idxs=3, xlabel="dimensionless time", ylabel="output", legend=false)

# ╔═╡ 9e32f370-a036-57bf-bb23-6a894c83dcd6
md"""Play $(@bind reed_playing CheckBox(default=false))"""

# ╔═╡ f68dbe74-d5bf-5394-b052-910ef4ddb106
# Cleanup runs when Play changes, this cell is rerun, or the notebook is shut down.
@use_effect([reed_source, reed_playing]) do
    if reed_playing
        device = get_default_output_device()
        device < 0 && error("No output device is available. Use the WAV examples or select an audio device with list_devices().")
        reset_state!(reed_source)
        start_DESource(reed_source, device; buffer_size=UInt32(2048))
    end
    return () -> begin
        isactive(reed_source) && stop_DESource(reed_source)
    end
end

# ╔═╡ 6c39f168-8bbf-5f8e-be63-6f9d6af262ff
md"""Uncheck Play to stop. For the struck model, stop and play again to repeat the strike. Live output uses the selected gain; offline WAV examples also filter and normalize the signal."""

# ╔═╡ 4d30f492-b596-5270-8636-0d3b2d123e46
md"""
## Bowed oscillator coupled to a body resonance: live audio

The same ODE drives `DESource`; channel mapping selects the audible state and `set_ts!` maps dimensionless time to seconds.
Choose Play to start, then move the controls; uncheck it to stop. Start with a low speaker volume.
"""

# ╔═╡ c808b4db-eac8-5b74-b094-a5628bd8c4f0
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

# ╔═╡ c75d8d45-5582-5741-a873-18a8d65aec0d
function bowed!(du, u, p, t)
    x, v, q, w = u
    F, V, ε, vs, δ, κ, Ω, ζ = p
    du[1] = v
    du[2] = -x - δ*v + F*tanh((V-v)/ε)/(1+((V-v)/vs)^2) + κ*(q-x)
    du[3] = w
    du[4] = -Ω^2*q - 2ζ*Ω*w + κ*(x-q)
    nothing
end

# ╔═╡ 41057cbe-c8bd-54ae-8f46-35521e1ad78b
bowed_u0 = [0.0, 0.0, 0.0, 0.0]

# ╔═╡ ab373653-ae42-5693-95e6-9ea0764e5e9f
bowed_p = [0.7, 0.3, 0.02, 0.25, 0.02, 0.08, 1.4, 0.02]  # F, V, ε, vs, δ, κ, Ω, ζ

# ╔═╡ cc4a7384-abe9-5ebb-bb33-2700b6417365
bowed_source = DESource(bowed!, copy(bowed_u0), copy(bowed_p); alg=Tsit5(), channel_map=[3, 3]);

# ╔═╡ d31f6002-3a7f-5d82-9f1a-df8d14513b7a
md"""
V $(@bind bowed_live_parameter Slider(0.1:0.05:0.6; default=0.3, show_value=true))

Pitch scale (Hz) $(@bind bowed_f0 Slider(110.0:1.0:330.0; default=196.0, show_value=true))

Gain $(@bind bowed_gain Slider(0.0:0.01:0.5; default=0.05, show_value=true))
"""

# ╔═╡ 515c4e5c-0425-5099-96df-129392b226a2
bowed_controls = begin
    set_param!(bowed_source, 2, Float64(bowed_live_parameter))
    set_ts!(bowed_source, 2pi*bowed_f0)
    set_gain!(bowed_source, Float64(bowed_gain))
    nothing
end

# ╔═╡ 2901b15f-f7d8-5592-af94-4550ae193322
# Compile a short solution before starting the audio callback.
bowed_preview = let
    bowed_controls
    bowed_prob = ODEProblem(bowed!, copy(bowed_u0), (0.0, 40.0), copy(get_params(bowed_source)))
    solve(bowed_prob, Tsit5(); saveat=0.02, abstol=1e-9, reltol=1e-7)
end;

# ╔═╡ 2879866e-8789-5cf2-8450-4c5bd3bd2cac
plot(bowed_preview; idxs=3, xlabel="dimensionless time", ylabel="output", legend=false)

# ╔═╡ 1764a78f-b3cc-5685-9b5e-9ef4c7617c7c
md"""Play $(@bind bowed_playing CheckBox(default=false))"""

# ╔═╡ 7781876b-7ded-5adc-8c99-6f48b0d6d33a
# Cleanup runs when Play changes, this cell is rerun, or the notebook is shut down.
@use_effect([bowed_source, bowed_playing]) do
    if bowed_playing
        device = get_default_output_device()
        device < 0 && error("No output device is available. Use the WAV examples or select an audio device with list_devices().")
        reset_state!(bowed_source)
        start_DESource(bowed_source, device; buffer_size=UInt32(2048))
    end
    return () -> begin
        isactive(bowed_source) && stop_DESource(bowed_source)
    end
end

# ╔═╡ b74dba24-6fee-59fa-b815-5b13f4869ae4
md"""Uncheck Play to stop. For the struck model, stop and play again to repeat the strike. Live output uses the selected gain; offline WAV examples also filter and normalize the signal."""

# ╔═╡ 095792a8-85ad-5f03-afaa-83d63238141a
md"""
## Struck nonlinear oscillator and resonator: live audio

The same ODE drives `DESource`; channel mapping selects the audible state and `set_ts!` maps dimensionless time to seconds.
Choose Play to start, then move the controls; uncheck it to stop. Start with a low speaker volume.
"""

# ╔═╡ caa73daf-b3a8-5911-93a4-5c2d444cb7c9
md"""
```math
\begin{aligned}
x' &= v, & v' &= -x-\alpha x^3-2\delta v+\kappa(q-x),\\
q' &= w, & w' &= -\Omega^2q-2\zeta\Omega w+\kappa(x-q).
\end{aligned}
```

Primes denote derivatives with respect to dimensionless time ``τ=2πf_0 t``.
"""

# ╔═╡ 18d44025-ceda-5d9d-bd8d-dfbf65136457
function struck!(du, u, p, t)
    x, v, q, w = u
    α, δ, κ, Ω, ζ = p
    du[1] = v
    du[2] = -x - α*x^3 - 2δ*v + κ*(q-x)
    du[3] = w
    du[4] = -Ω^2*q - 2ζ*Ω*w + κ*(x-q)
    nothing
end

# ╔═╡ 266e8c9c-380a-58f7-ab21-3d642cbf185b
struck_u0 = [0.0, 1.0, 0.0, 0.0]

# ╔═╡ fe64e557-b8c7-5d7b-9149-de36d5a38d37
struck_p = [0.8, 0.001, 0.08, 2.4, 0.0015]  # α, δ, κ, Ω, ζ

# ╔═╡ 1f29bef6-7c16-5a31-9560-f04cec1b04a6
struck_source = DESource(struck!, copy(struck_u0), copy(struck_p); alg=Tsit5(), channel_map=[3, 3]);

# ╔═╡ ff2fe376-f596-551d-94fb-f8c850c308dd
md"""
α $(@bind struck_live_parameter Slider(0.0:0.1:1.5; default=0.8, show_value=true))

Pitch scale (Hz) $(@bind struck_f0 Slider(110.0:1.0:330.0; default=220.0, show_value=true))

Gain $(@bind struck_gain Slider(0.0:0.01:0.5; default=0.05, show_value=true))
"""

# ╔═╡ 85485ee6-ac3d-538b-876d-34669f4d822d
struck_controls = begin
    set_param!(struck_source, 1, Float64(struck_live_parameter))
    set_ts!(struck_source, 2pi*struck_f0)
    set_gain!(struck_source, Float64(struck_gain))
    nothing
end

# ╔═╡ 3014026a-0f89-523f-8a9e-5b8b3e551ce9
# Compile a short solution before starting the audio callback.
struck_preview = let
    struck_controls
    struck_prob = ODEProblem(struck!, copy(struck_u0), (0.0, 40.0), copy(get_params(struck_source)))
    solve(struck_prob, Tsit5(); saveat=0.02, abstol=1e-9, reltol=1e-7)
end;

# ╔═╡ f277ef43-83f4-5e52-9499-5df02732644e
plot(struck_preview; idxs=3, xlabel="dimensionless time", ylabel="output", legend=false)

# ╔═╡ cc2a5dba-9cec-5c16-8c09-5198f465b450
md"""Play $(@bind struck_playing CheckBox(default=false))"""

# ╔═╡ 948d0bd2-863c-5909-b003-ef88febe8e44
# Cleanup runs when Play changes, this cell is rerun, or the notebook is shut down.
@use_effect([struck_source, struck_playing]) do
    if struck_playing
        device = get_default_output_device()
        device < 0 && error("No output device is available. Use the WAV examples or select an audio device with list_devices().")
        reset_state!(struck_source)
        start_DESource(struck_source, device; buffer_size=UInt32(2048))
    end
    return () -> begin
        isactive(struck_source) && stop_DESource(struck_source)
    end
end

# ╔═╡ 2abcf345-b699-5488-bbc6-c91fa154442a
md"""Uncheck Play to stop. For the struck model, stop and play again to repeat the strike. Live output uses the selected gain; offline WAV examples also filter and normalize the signal."""

# ╔═╡ Cell order:
# ╠═c9467c68-6aaa-59d0-8fde-f3b87b04a879
# ╟─4c2d0c28-f707-55b7-abf3-9b633eb0fe26
# ╟─3e839d74-f0b5-5cb6-837a-7d6b38bc1bb2
# ╟─3b98f12f-a3ee-58c6-8eb6-2833f7f6d832
# ╟─beaa5691-55c0-5135-b085-a8104614feee
# ╠═79241bbc-6ee9-5a8a-a14a-71fc0db879bc
# ╠═659ce933-3134-5402-a72f-e40c06f4f32d
# ╠═63c6d438-8adf-52d7-a6e2-c2e3cb2d5c88
# ╠═7bcdb33c-411e-5bdb-8331-34ad56286122
# ╟─773429e4-9f6e-5761-81c0-82c3a8de4a9c
# ╠═b57d251a-6a70-5867-9eb7-40489a888a0e
# ╠═9408ed8f-19c9-5c38-9f8f-56ac6b03daea
# ╠═a4002915-43b7-5e55-8539-14cada302fed
# ╟─9e32f370-a036-57bf-bb23-6a894c83dcd6
# ╠═f68dbe74-d5bf-5394-b052-910ef4ddb106
# ╟─6c39f168-8bbf-5f8e-be63-6f9d6af262ff
# ╟─4d30f492-b596-5270-8636-0d3b2d123e46
# ╟─c808b4db-eac8-5b74-b094-a5628bd8c4f0
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
# ╟─095792a8-85ad-5f03-afaa-83d63238141a
# ╟─caa73daf-b3a8-5911-93a4-5c2d444cb7c9
# ╠═18d44025-ceda-5d9d-bd8d-dfbf65136457
# ╠═266e8c9c-380a-58f7-ab21-3d642cbf185b
# ╠═fe64e557-b8c7-5d7b-9149-de36d5a38d37
# ╠═1f29bef6-7c16-5a31-9560-f04cec1b04a6
# ╟─ff2fe376-f596-551d-94fb-f8c850c308dd
# ╠═85485ee6-ac3d-538b-876d-34669f4d822d
# ╠═3014026a-0f89-523f-8a9e-5b8b3e551ce9
# ╠═f277ef43-83f4-5e52-9499-5df02732644e
# ╟─cc2a5dba-9cec-5c16-8c09-5198f465b450
# ╠═948d0bd2-863c-5909-b003-ef88febe8e44
# ╟─2abcf345-b699-5488-bbc6-c91fa154442a
