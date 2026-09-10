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

# ╔═╡ 1094b70d-1ec8-52eb-a5d5-361e23b0d1ad
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots, PlutoUI
    using RealTimeAudioDiffEq, PlutoHooks
end

# ╔═╡ 095792a8-85ad-5f03-afaa-83d63238141a
md"""
# Struck nonlinear oscillator and resonator: live audio

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
function model!(du, u, p, t)
    x, v, q, w = u
    α, δ, κ, Ω, ζ = p
    du[1] = v
    du[2] = -x - α*x^3 - 2δ*v + κ*(q-x)
    du[3] = w
    du[4] = -Ω^2*q - 2ζ*Ω*w + κ*(x-q)
    nothing
end

# ╔═╡ 266e8c9c-380a-58f7-ab21-3d642cbf185b
u0 = [0.0, 1.0, 0.0, 0.0]

# ╔═╡ fe64e557-b8c7-5d7b-9149-de36d5a38d37
p = [0.8, 0.001, 0.08, 2.4, 0.0015]  # α, δ, κ, Ω, ζ

# ╔═╡ 1f29bef6-7c16-5a31-9560-f04cec1b04a6
source = DESource(model!, copy(u0), copy(p); alg=Tsit5(), channel_map=[3, 3]);

# ╔═╡ ff2fe376-f596-551d-94fb-f8c850c308dd
md"""
α $(@bind live_parameter Slider(0.0:0.1:1.5; default=0.8, show_value=true))

Pitch scale (Hz) $(@bind f0 Slider(110.0:1.0:330.0; default=220.0, show_value=true))

Gain $(@bind gain Slider(0.0:0.01:0.5; default=0.05, show_value=true))
"""

# ╔═╡ 85485ee6-ac3d-538b-876d-34669f4d822d
controls = begin
    set_param!(source, 1, Float64(live_parameter))
    set_ts!(source, 2pi*f0)
    set_gain!(source, Float64(gain))
    nothing
end

# ╔═╡ 3014026a-0f89-523f-8a9e-5b8b3e551ce9
# Compile a short solution before starting the audio callback.
preview = let
    controls
    prob = ODEProblem(model!, copy(u0), (0.0, 40.0), copy(get_params(source)))
    solve(prob, Tsit5(); saveat=0.02, abstol=1e-9, reltol=1e-7)
end;

# ╔═╡ f277ef43-83f4-5e52-9499-5df02732644e
plot(preview; idxs=3, xlabel="dimensionless time", ylabel="output", legend=false)

# ╔═╡ cc2a5dba-9cec-5c16-8c09-5198f465b450
md"""Play $(@bind playing CheckBox(default=false))"""

# ╔═╡ 948d0bd2-863c-5909-b003-ef88febe8e44
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

# ╔═╡ 2abcf345-b699-5488-bbc6-c91fa154442a
md"""Uncheck Play to stop. For the struck model, stop and play again to repeat the strike. Live output uses the selected gain; offline WAV examples also filter and normalize the signal."""

# ╔═╡ Cell order:
# ╠═1094b70d-1ec8-52eb-a5d5-361e23b0d1ad
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
