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

# ╔═╡ 701cd6c8-ec63-534b-9a14-e1a60093fbd6
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots, PlutoUI, RealTimeAudioDiffEq, PlutoHooks
end

# ╔═╡ 10ffec0e-e17b-5d51-8f10-dd2164113596
md"""
# The Lorenz system: live audio

Accelerate the Lorenz trajectory and listen to x/25. The fixed scaling
keeps the plotted state and the sound amplitude distinct. The parameter ρ
changes the flow; the pitch scale changes playback speed.
"""

# ╔═╡ ac7e0cc9-3e63-548e-963d-2d8df818540c
md"""
```math
\begin{aligned}
\dot{x} &= \sigma(y-x),\\
\dot{y} &= x(\rho-z)-y,\\
\dot{z} &= xy-\beta z.
\end{aligned}
```
"""

# ╔═╡ a1370801-68fb-5182-9f42-864f1d99b2bf
function model!(du, u, p, t)
    x, y, z = u
    σ, ρ, β = p
    du[1] = σ*(y-x)
    du[2] = x*(ρ-z)-y
    du[3] = x*y-β*z
    nothing
end

# ╔═╡ 2759e95f-90c5-59a6-b171-c84b564f8e96
u0 = [1.0, 0.0, 0.0]

# ╔═╡ f08907b5-4e63-553f-bf69-64de0aa92d15
p = [10.0, 28.0, 8/3]  # σ, ρ, β

# ╔═╡ b29b6c2a-697f-59f8-b153-0a589dce8c01
source = let
    audio_source = DESource(model!, copy(u0), copy(p);
        alg=Tsit5(), channel_map=[1/25 1/25; 0.0 0.0; 0.0 0.0])
    audio_source.data.problem = remake(audio_source.data.problem;
        abstol=1e-9, reltol=1e-7)
    audio_source
end;

# ╔═╡ 31cf593c-e13b-5468-886c-51f538b58458
md"""
ρ $(@bind rho Slider(1.0:1.0:40.0; default=28.0, show_value=true))

Pitch scale (Hz) $(@bind f0 Slider(1.0:1.0:30.0; default=10.0, show_value=true))

Gain $(@bind gain Slider(0.0:0.01:0.2; default=0.05, show_value=true))
"""

# ╔═╡ 349eb1da-8fbd-5491-bc67-5209945fcac0
controls = begin
    set_param!(source, 2, Float64(rho))
    set_ts!(source, Float64(2pi*f0))
    set_gain!(source, Float64(gain))
    nothing
end

# ╔═╡ e2967db0-aece-51d0-8ef6-76893b7dada9
preview = let
    controls
    prob = remake(source.data.problem; u0=copy(u0), p=copy(get_params(source)), tspan=(0.0, 60.0))
    solve(prob, Tsit5(); saveat=0.01)
end;

# ╔═╡ e19fbe23-95e6-514c-9cbf-3149a92570db
plot(preview; idxs=[1,2,3], xlabel="simulation time", ylabel="state",
    margin=5 * Plots.mm)

# ╔═╡ d6d2edd8-b2f2-5bd5-9370-9e00f8c6c687
md"""
Play $(@bind playing CheckBox(default=false))
"""

# ╔═╡ ab81bcbd-3c1e-5d3a-952b-2c8b64f70783
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

# ╔═╡ 7f52cc4d-50b8-5059-ac54-747fe7639367
md"""
Uncheck Play to stop. Stop and play again to restart from the initial condition.
Live output uses the selected gain; closing the notebook in Pluto stops its source.
"""

# ╔═╡ Cell order:
# ╠═701cd6c8-ec63-534b-9a14-e1a60093fbd6
# ╟─10ffec0e-e17b-5d51-8f10-dd2164113596
# ╟─ac7e0cc9-3e63-548e-963d-2d8df818540c
# ╠═a1370801-68fb-5182-9f42-864f1d99b2bf
# ╠═2759e95f-90c5-59a6-b171-c84b564f8e96
# ╠═f08907b5-4e63-553f-bf69-64de0aa92d15
# ╠═b29b6c2a-697f-59f8-b153-0a589dce8c01
# ╟─31cf593c-e13b-5468-886c-51f538b58458
# ╠═349eb1da-8fbd-5491-bc67-5209945fcac0
# ╠═e2967db0-aece-51d0-8ef6-76893b7dada9
# ╠═e19fbe23-95e6-514c-9cbf-3149a92570db
# ╟─d6d2edd8-b2f2-5bd5-9370-9e00f8c6c687
# ╠═ab81bcbd-3c1e-5d3a-952b-2c8b64f70783
# ╟─7f52cc4d-50b8-5059-ac54-747fe7639367
