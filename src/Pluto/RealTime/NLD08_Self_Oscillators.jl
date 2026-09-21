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

# ╔═╡ e73a2692-e982-5d9a-82d9-bb238110aa6d
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, Plots, PlutoUI, RealTimeAudioDiffEq, PlutoHooks
end

# ╔═╡ 702d0973-ea3d-5adf-96b7-b6ddef104835
md"""
# Self-oscillators: live audio

Compare both examples from ODE NLD08. Negative damping sustains the Rayleigh
reed; a smooth cubic slip law supplies the rubbed oscillator's excitation.
The latter illustrates self-oscillation without imposing exact sticking.
"""

# ╔═╡ 00992aef-ef8a-51f7-99cf-8a4bd9ef51d4
md"""
The first-order equations are

```math
\begin{aligned}
\dot{x} &= v,\\
\dot{v} &= -Kx-\gamma(v)v=-Kx+(\mu-s_c^2v^2)v.
\end{aligned}
```

The cubic term $-s_c^2v^3$ makes the system nonlinear. The function below writes
these two derivatives into `du`, with `u = [x, v]` and `p = [μ, K, sc]`.
"""

# ╔═╡ 5a57db75-c11b-5e82-a561-0de066e27013
function model!(du, u, p, t)
    x, v = u
    μ, K, sc = p
    du[1] = v
    du[2] = -K*x + (μ-sc^2*v^2)*v
    nothing
end

# ╔═╡ d523d95e-12c7-5de1-afba-37dbafaca77a
md"""
For belt speed $V$ and slip $s=v-V$, use $C(s)=-s(1-s^2)$:

```math
x'=v,\qquad v'=-x-\mu C(v-V).
```

The first state is displacement and the second is velocity. Both examples use
the pitch scale $\tau=2\pi f_0t$.
"""

# ╔═╡ 93f2ab27-4c99-5e0a-b892-1600243e4eca
friction(s) = -s*(1-s^2)

# ╔═╡ 8bfc12f6-34c4-53b8-9158-153f3e584418
function bow!(du, u, p, t)
    x, v = u
    μ, V = p
    du[1] = v
    du[2] = -x - μ * friction(v - V)
    nothing
end

# ╔═╡ cf23f192-89e8-519c-b3b3-788b56ada0d8
md"""
Example $(@bind example Select(["Rayleigh reed", "Rubbed oscillator"]; default="Rayleigh reed"))
"""

# ╔═╡ bdc5b1c7-d94d-5701-add6-df3ba2c8d3bd
u0 = example == "Rayleigh reed" ? [0.1,0.0] : [0.7,0.0]

# ╔═╡ c6e6394f-a8dd-5b96-b749-63862ccafba3
p = example == "Rayleigh reed" ? [0.3,1.0,1.0] : [0.02,0.02]

# ╔═╡ e610e593-85b7-50f9-bcb3-86aadca75484
chosen_model! = example == "Rayleigh reed" ? model! : bow!

# ╔═╡ 8fe2e0a2-d5ae-56f6-b758-962463f268c7
source = let
    audio_source = DESource(chosen_model!, copy(u0), copy(p);
        alg=Tsit5(), channel_map=[1,1])
    audio_source.data.problem = remake(audio_source.data.problem;
        abstol=1e-9, reltol=1e-7)
    audio_source
end;

# ╔═╡ 5cd9bfab-f22c-52c1-bd9c-56c6623753d0
md"""
Excitation μ $(@bind drive Slider(0.0:0.01:0.5; default=p[1], show_value=true))

Pitch scale (Hz) $(@bind f0 Slider(80.0:5.0:400.0; default=220.0, show_value=true))

Gain $(@bind gain Slider(0.0:0.01:0.2; default=0.05, show_value=true))

The other parameters remain visible in p. Choose the example before adjusting
its excitation; switching examples restarts the corresponding source.
"""

# ╔═╡ d90ec5ca-9b20-5f9e-81ef-b94613775fb9
controls = begin
    set_param!(source, 1, Float64(drive))
    set_ts!(source, 2pi*f0)
    set_gain!(source, Float64(gain))
    nothing
end

# ╔═╡ f800327c-136f-549a-ab01-687c4aca1dba
preview = let
    controls
    prob = remake(source.data.problem; tspan=(0.0,250.0), p=copy(get_params(source)))
    solve(prob, Tsit5(); saveat=0.02)
end;

# ╔═╡ e722356c-8218-5e24-839b-7af760d2b9cb
let
    waveform = plot(preview; idxs=1, xlabel="simulation time", ylabel="x", legend=false)
    orbit = plot(preview; idxs=(1,2), xlabel="x", ylabel="v", legend=false)
    plot(waveform, orbit; layout=(1,2), size=(850,320), margin=5 * Plots.mm)
end

# ╔═╡ e1d37e7e-a0e5-54ee-9748-38bac0428680
md"""
Play $(@bind playing CheckBox(default=false))
"""

# ╔═╡ 0bd376b4-8165-5139-ab18-7e26b7324d09
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

# ╔═╡ 5d42e7e0-ffa3-5d3e-8ef0-096fbb09e603
md"""
Uncheck Play to stop. Stop and play again to restart from the initial condition.
Live output uses the selected gain; closing the notebook in Pluto stops its source.
"""

# ╔═╡ Cell order:
# ╠═e73a2692-e982-5d9a-82d9-bb238110aa6d
# ╟─702d0973-ea3d-5adf-96b7-b6ddef104835
# ╟─00992aef-ef8a-51f7-99cf-8a4bd9ef51d4
# ╠═5a57db75-c11b-5e82-a561-0de066e27013
# ╟─d523d95e-12c7-5de1-afba-37dbafaca77a
# ╠═93f2ab27-4c99-5e0a-b892-1600243e4eca
# ╠═8bfc12f6-34c4-53b8-9158-153f3e584418
# ╟─cf23f192-89e8-519c-b3b3-788b56ada0d8
# ╠═bdc5b1c7-d94d-5701-add6-df3ba2c8d3bd
# ╠═c6e6394f-a8dd-5b96-b749-63862ccafba3
# ╠═e610e593-85b7-50f9-bcb3-86aadca75484
# ╠═8fe2e0a2-d5ae-56f6-b758-962463f268c7
# ╟─5cd9bfab-f22c-52c1-bd9c-56c6623753d0
# ╠═d90ec5ca-9b20-5f9e-81ef-b94613775fb9
# ╠═f800327c-136f-549a-ab01-687c4aca1dba
# ╠═e722356c-8218-5e24-839b-7af760d2b9cb
# ╟─e1d37e7e-a0e5-54ee-9748-38bac0428680
# ╠═0bd376b4-8165-5139-ab18-7e26b7324d09
# ╟─5d42e7e0-ffa3-5d3e-8ef0-096fbb09e603
