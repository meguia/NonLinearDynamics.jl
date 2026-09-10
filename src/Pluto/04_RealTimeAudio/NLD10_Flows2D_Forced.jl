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

# ╔═╡ 8866f0ab-f466-5785-b9ae-bcbdae6a0314
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots, PlutoUI, RealTimeAudioDiffEq, PlutoHooks
end

# ╔═╡ 499d1360-81da-59b8-a178-05b6171d901c
md"""
# Periodically forced Duffing oscillator: live audio

The same ODE drives `DESource`; channel mapping selects the audible state and `set_ts!` maps dimensionless time to seconds.
Choose Play to start, then move the controls; uncheck it to stop. Start with a low speaker volume.
"""

# ╔═╡ c17144d0-b559-5774-9c5e-08d0680864de
md"""
```math
\begin{aligned}
\dot{x} &= v, & \dot{v} &= -\gamma v+\beta x-x^3+F\cos(\omega t).
\end{aligned}
```
"""

# ╔═╡ 15cf10e6-1168-5130-bb25-78f02cf1321e
function model!(du, u, p, t)
    x, v = u
    γ, β, F, ω = p
    du[1] = v
    du[2] = -γ*v + β*x - x^3 + F*cos(ω*t)
    nothing
end

# ╔═╡ 21b109e3-9147-51d4-891b-fc0edfb9e1c4
u0 = [0.1, 0.0]

# ╔═╡ 05b77296-0f31-539b-9abf-61f8362b1670
p = [0.2, 1.0, 0.3, 1.2]  # γ, β, F, ω

# ╔═╡ 0a3bd136-a666-5d68-b1aa-a1467fb701d7
source = DESource(model!, copy(u0), copy(p); alg=Tsit5(), channel_map=[1, 1]);

# ╔═╡ f7948483-bee8-5b6c-ac72-74900e772823
md"""
F $(@bind live_parameter Slider(0.1:0.05:0.6; default=0.3, show_value=true))

Pitch scale (Hz) $(@bind f0 Slider(110.0:1.0:330.0; default=220.0, show_value=true))

Gain $(@bind gain Slider(0.0:0.01:0.5; default=0.05, show_value=true))
"""

# ╔═╡ 2f454b68-bf8c-5d59-9665-f00b90809611
controls = begin
    set_param!(source, 3, Float64(live_parameter))
    set_ts!(source, 2pi*f0)
    set_gain!(source, Float64(gain))
    nothing
end

# ╔═╡ f9d9288f-f250-5a0d-bffc-fc4f95c614c6
# Compile a short solution before starting the audio callback.
preview = let
    controls
    prob = ODEProblem(model!, copy(u0), (0.0, 40.0), copy(get_params(source)))
    solve(prob, Tsit5(); saveat=0.02, abstol=1e-9, reltol=1e-7)
end;

# ╔═╡ e401d243-b094-55d1-9ed7-037050bf779d
plot(preview; idxs=1, xlabel="dimensionless time", ylabel="output", legend=false)

# ╔═╡ 7dd0941d-872d-5fd3-97bf-f723ed3fee2d
md"""Play $(@bind playing CheckBox(default=false))"""

# ╔═╡ cc5f029d-668f-5186-94e9-2609d321ee16
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

# ╔═╡ 88d1371e-215c-560a-9a35-66875450ad22
md"""Uncheck Play to stop. For the struck model, stop and play again to repeat the strike. Live output uses the selected gain; offline WAV examples also filter and normalize the signal."""

# ╔═╡ Cell order:
# ╠═8866f0ab-f466-5785-b9ae-bcbdae6a0314
# ╟─499d1360-81da-59b8-a178-05b6171d901c
# ╟─c17144d0-b559-5774-9c5e-08d0680864de
# ╠═15cf10e6-1168-5130-bb25-78f02cf1321e
# ╠═21b109e3-9147-51d4-891b-fc0edfb9e1c4
# ╠═05b77296-0f31-539b-9abf-61f8362b1670
# ╠═0a3bd136-a666-5d68-b1aa-a1467fb701d7
# ╟─f7948483-bee8-5b6c-ac72-74900e772823
# ╠═2f454b68-bf8c-5d59-9665-f00b90809611
# ╠═f9d9288f-f250-5a0d-bffc-fc4f95c614c6
# ╠═e401d243-b094-55d1-9ed7-037050bf779d
# ╟─7dd0941d-872d-5fd3-97bf-f723ed3fee2d
# ╠═cc5f029d-668f-5186-94e9-2609d321ee16
# ╟─88d1371e-215c-560a-9a35-66875450ad22
