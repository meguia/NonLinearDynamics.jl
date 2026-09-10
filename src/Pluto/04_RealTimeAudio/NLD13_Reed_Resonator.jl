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

# ╔═╡ 3b98f12f-a3ee-58c6-8eb6-2833f7f6d832
md"""
# Reed oscillator coupled to a bore resonance: live audio

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
function model!(du, u, p, t)
    x, v, q, w = u
    μ, κ, Ω, ζ = p
    du[1] = v
    du[2] = -x + μ*(1-v^2)*v + κ*(q-x)
    du[3] = w
    du[4] = -Ω^2*q - 2ζ*Ω*w + κ*(x-q)
    nothing
end

# ╔═╡ 659ce933-3134-5402-a72f-e40c06f4f32d
u0 = [0.05, 0.0, 0.0, 0.0]

# ╔═╡ 63c6d438-8adf-52d7-a6e2-c2e3cb2d5c88
p = [0.3, 0.12, 1.02, 0.025]  # μ, κ, Ω, ζ

# ╔═╡ 7bcdb33c-411e-5bdb-8331-34ad56286122
source = DESource(model!, copy(u0), copy(p); alg=Tsit5(), channel_map=[3, 3]);

# ╔═╡ 773429e4-9f6e-5761-81c0-82c3a8de4a9c
md"""
μ $(@bind live_parameter Slider(0.0:0.05:0.8; default=0.3, show_value=true))

Pitch scale (Hz) $(@bind f0 Slider(110.0:1.0:330.0; default=220.0, show_value=true))

Gain $(@bind gain Slider(0.0:0.01:0.5; default=0.05, show_value=true))
"""

# ╔═╡ b57d251a-6a70-5867-9eb7-40489a888a0e
controls = begin
    set_param!(source, 1, Float64(live_parameter))
    set_ts!(source, 2pi*f0)
    set_gain!(source, Float64(gain))
    nothing
end

# ╔═╡ 9408ed8f-19c9-5c38-9f8f-56ac6b03daea
# Compile a short solution before starting the audio callback.
preview = let
    controls
    prob = ODEProblem(model!, copy(u0), (0.0, 40.0), copy(get_params(source)))
    solve(prob, Tsit5(); saveat=0.02, abstol=1e-9, reltol=1e-7)
end;

# ╔═╡ a4002915-43b7-5e55-8539-14cada302fed
plot(preview; idxs=3, xlabel="dimensionless time", ylabel="output", legend=false)

# ╔═╡ 9e32f370-a036-57bf-bb23-6a894c83dcd6
md"""Play $(@bind playing CheckBox(default=false))"""

# ╔═╡ f68dbe74-d5bf-5394-b052-910ef4ddb106
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

# ╔═╡ 6c39f168-8bbf-5f8e-be63-6f9d6af262ff
md"""Uncheck Play to stop. For the struck model, stop and play again to repeat the strike. Live output uses the selected gain; offline WAV examples also filter and normalize the signal."""

# ╔═╡ Cell order:
# ╠═c9467c68-6aaa-59d0-8fde-f3b87b04a879
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
