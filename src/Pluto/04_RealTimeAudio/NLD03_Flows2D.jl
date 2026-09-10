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

# ╔═╡ d47730f1-8075-5431-b1dc-2caaea9e94af
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots, PlutoUI, RealTimeAudioDiffEq, PlutoHooks
end

# ╔═╡ 97c5f81a-77b3-5a77-ac9e-35b592e5550e
md"""
# Flows in 2D: damped oscillator: live audio

The same ODE drives `DESource`; channel mapping selects the audible state and `set_ts!` maps dimensionless time to seconds.
Choose Play to start, then move the controls; uncheck it to stop. Start with a low speaker volume.
"""

# ╔═╡ 929a6e41-527a-502a-9f12-a172521a3369
function model!(du, u, p, t)
    x, v = u
    k, γ = p
    du[1] = v
    du[2] = -k*x - γ*v
    nothing
end

# ╔═╡ b48824df-388a-5d36-9fcf-e69127c5c971
u0 = [1.0, 0.0]

# ╔═╡ be36ccd4-f6e9-5e25-8ed0-6132b4484151
p = [1.0, 0.0]  # k, γ; zero damping sustains this introductory tone

# ╔═╡ 26335d90-9716-5cbf-9ae9-faddc9a995dc
source = DESource(model!, copy(u0), copy(p); alg=Tsit5(), channel_map=[1, 1]);

# ╔═╡ 3c02880f-0626-5ba8-a480-e8ef33b69864
md"""
k $(@bind live_parameter Slider(0.5:0.05:2.0; default=1.0, show_value=true))

Pitch scale (Hz) $(@bind f0 Slider(110.0:1.0:330.0; default=220.0, show_value=true))

Gain $(@bind gain Slider(0.0:0.01:0.5; default=0.05, show_value=true))
"""

# ╔═╡ b1434dc0-a3cd-589c-a33f-a5aaf6b0a3c9
controls = begin
    set_param!(source, 1, Float64(live_parameter))
    set_ts!(source, 2pi*f0)
    set_gain!(source, Float64(gain))
    nothing
end

# ╔═╡ 8a80c215-70c2-5873-91db-a9fddf7ae2b2
# Compile a short solution before starting the audio callback.
preview = let
    controls
    prob = ODEProblem(model!, copy(u0), (0.0, 40.0), copy(get_params(source)))
    solve(prob, Tsit5(); saveat=0.02, abstol=1e-9, reltol=1e-7)
end;

# ╔═╡ f800d62c-5841-5de5-b3cd-1410127743ff
plot(preview; idxs=1, xlabel="dimensionless time", ylabel="output", legend=false)

# ╔═╡ 8dd195ed-c4c8-569c-9008-5aa7c422e22a
md"""Play $(@bind playing CheckBox(default=false))"""

# ╔═╡ 5da11c70-99c1-5b27-9e12-a2309d5faea0
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

# ╔═╡ 2ce158d7-6d5b-5378-955b-700fd51e5954
md"""Uncheck Play to stop. For the struck model, stop and play again to repeat the strike. Live output uses the selected gain; offline WAV examples also filter and normalize the signal."""

# ╔═╡ Cell order:
# ╠═d47730f1-8075-5431-b1dc-2caaea9e94af
# ╟─97c5f81a-77b3-5a77-ac9e-35b592e5550e
# ╠═929a6e41-527a-502a-9f12-a172521a3369
# ╠═b48824df-388a-5d36-9fcf-e69127c5c971
# ╠═be36ccd4-f6e9-5e25-8ed0-6132b4484151
# ╠═26335d90-9716-5cbf-9ae9-faddc9a995dc
# ╟─3c02880f-0626-5ba8-a480-e8ef33b69864
# ╠═b1434dc0-a3cd-589c-a33f-a5aaf6b0a3c9
# ╠═8a80c215-70c2-5873-91db-a9fddf7ae2b2
# ╠═f800d62c-5841-5de5-b3cd-1410127743ff
# ╟─8dd195ed-c4c8-569c-9008-5aa7c422e22a
# ╠═5da11c70-99c1-5b27-9e12-a2309d5faea0
# ╟─2ce158d7-6d5b-5378-955b-700fd51e5954
