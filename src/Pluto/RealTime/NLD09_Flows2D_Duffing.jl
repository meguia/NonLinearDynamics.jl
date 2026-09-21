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

# ╔═╡ 1b3706f3-b912-5d78-b278-7e422e3ca568
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots, PlutoUI, RealTimeAudioDiffEq, PlutoHooks
end

# ╔═╡ 362aa14d-d1df-5b69-b6d1-c1785a9bbcdd
md"""
# The Duffing oscillator: live audio

The same double-well oscillator as in ODE drives the audio. Listen to
velocity so that a resting displacement does not produce DC. Damping makes
the tone decay; set γ = 0 for conservative oscillation.
"""

# ╔═╡ ab526a23-527c-5168-a1b6-a95bd1c39fc9
md"""
```math
\begin{aligned}
\dot{x} &= v, & \dot{v} &= -\gamma v+\beta x-x^3.
\end{aligned}
```
"""

# ╔═╡ 03027286-e51a-5db3-9391-3c3995a9837e
function model!(du, u, p, t)
    x, v = u
    γ, β = p
    du[1] = v
    du[2] = -γ*v + β*x - x^3
    nothing
end

# ╔═╡ b7d6b9f1-b249-55da-aaed-0dad2553f0ff
u0 = [0.1, 0.7]

# ╔═╡ 1fff0d46-e540-5af6-b6fc-52c4aed4abe6
p = [0.15, 1.0]  # γ, β

# ╔═╡ 8073e0d7-dcfe-5349-9182-62a737d8f23d
source = let
    audio_source = DESource(model!, copy(u0), copy(p);
        alg=Tsit5(), channel_map=[2,2])
    audio_source.data.problem = remake(audio_source.data.problem;
        abstol=1e-9, reltol=1e-7)
    audio_source
end;

# ╔═╡ cff76f41-4a0d-535b-bf1e-03cf3049ceb1
md"""
Damping γ $(@bind damping Slider(0.0:0.01:0.3; default=0.15, show_value=true))

Pitch scale (Hz) $(@bind f0 Slider(80.0:5.0:400.0; default=220.0, show_value=true))

Gain $(@bind gain Slider(0.0:0.01:0.2; default=0.05, show_value=true))
"""

# ╔═╡ 2af58a1a-ec00-58f8-9e4e-dc3f42db044f
controls = begin
    set_param!(source, 1, Float64(damping))
    set_ts!(source, Float64(2pi*f0))
    set_gain!(source, Float64(gain))
    nothing
end

# ╔═╡ ad7fb115-3a4a-5484-aca8-ffb9dc46d775
preview = let
    controls
    prob = remake(source.data.problem; u0=copy(u0), p=copy(get_params(source)), tspan=(0.0, 60.0))
    solve(prob, Tsit5(); saveat=0.02)
end;

# ╔═╡ 627372e8-362b-5125-910f-f4215b5b0a17
plot(preview; idxs=[1, 2], xlabel="simulation time", ylabel="state",
    margin=5 * Plots.mm)

# ╔═╡ 947aee4c-95a3-5ccd-8e23-ddb3877b10e1
md"""
Play $(@bind playing CheckBox(default=false))
"""

# ╔═╡ ade44a0e-6266-52b8-b840-3c2e6a354a5e
@use_effect([source, playing]) do
    if playing
        device = get_default_output_device()
        device < 0 && error("No output device is available. Select an audio device with list_devices().")
        reset_state!(source)

        start_DESource(source, device; buffer_size=UInt32(2048))
    end
    return () -> begin
        isactive(source) && stop_DESource(source)
    end
end

# ╔═╡ bf5e1fa6-41cd-509e-929d-c634df88f6d7
md"""
Uncheck Play to stop. Stop and play again to restart from the initial condition.
Live output uses the selected gain; closing the notebook in Pluto stops its source.
"""

# ╔═╡ Cell order:
# ╠═1b3706f3-b912-5d78-b278-7e422e3ca568
# ╟─362aa14d-d1df-5b69-b6d1-c1785a9bbcdd
# ╟─ab526a23-527c-5168-a1b6-a95bd1c39fc9
# ╠═03027286-e51a-5db3-9391-3c3995a9837e
# ╠═b7d6b9f1-b249-55da-aaed-0dad2553f0ff
# ╠═1fff0d46-e540-5af6-b6fc-52c4aed4abe6
# ╠═8073e0d7-dcfe-5349-9182-62a737d8f23d
# ╟─cff76f41-4a0d-535b-bf1e-03cf3049ceb1
# ╠═2af58a1a-ec00-58f8-9e4e-dc3f42db044f
# ╠═ad7fb115-3a4a-5484-aca8-ffb9dc46d775
# ╠═627372e8-362b-5125-910f-f4215b5b0a17
# ╟─947aee4c-95a3-5ccd-8e23-ddb3877b10e1
# ╠═ade44a0e-6266-52b8-b840-3c2e6a354a5e
# ╟─bf5e1fa6-41cd-509e-929d-c634df88f6d7
