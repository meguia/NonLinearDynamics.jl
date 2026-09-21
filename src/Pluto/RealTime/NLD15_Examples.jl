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

# ╔═╡ 64a40c39-28e8-5115-977e-d4c329c5badc
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
    Pkg.instantiate()
    using DifferentialEquations, LinearAlgebra, Plots, PlutoUI, RealTimeAudioDiffEq, PlutoHooks
    include(joinpath(@__DIR__, "..", "..", "realtime_delays.jl"))
    using .RealtimeDelays: buffered_dde_source, warmup_delay_source!
end

# ╔═╡ 20db5510-3672-55d0-b7cb-e35c23e731f2
md"""
# Examples: vocal folds and birdsong

Two models connect nonlinear oscillation to vocal sound. Begin with phonation
onset and saturation in a reduced vocal-fold oscillator, then follow a
syringeal source through a delayed vocal tract and an acoustic filter.

Each subsection has its own controls and Play switch. Stop one before
playing the next to compare the sounds.
"""

# ╔═╡ 4de90467-b82e-5c5b-9ea7-966ab8a026f2
TableOfContents()

# ╔═╡ 4de0c122-664f-5f68-aae3-8097e02f400e
md"""
## Vocal folds: live audio

A reduced vocal-fold oscillator combines drive, stiffness, a constant bias,
and nonlinear damping. These are the same equations and initial conditions
as the vocal-fold section in ODE NLD15; listen to displacement while changing drive and stiffness.
"""

# ╔═╡ e02194ae-76a5-5b91-92d4-a19fc9878d6f
md"""
With dimensionless time $\tau=2\pi f_0t$, define the effective stiffness
$K_{\mathrm{eff}}=k+0.1\mu(k+1)$. The equations are

```math
\begin{aligned}
x' &= v,\\
v' &= -K_{\mathrm{eff}}x+(\mu-x^2)v-F.
\end{aligned}
```

The state is `u = [x, v]` and the parameters are `p = [μ, k, F]`.
For $K_{\mathrm{eff}}>0$, the equilibrium is $x_*=-F/K_{\mathrm{eff}}$, $v_*=0$.
Its linearized damping changes sign at $\mu-x_*^2=0$; a bias can therefore
shift the onset of oscillation. The scale $f_0$ converts simulation time to
seconds; the resulting pitch also depends on stiffness and drive.
"""

# ╔═╡ 061cf701-efa1-5a0b-8db4-b6b7464c9c9d
function vocal!(du, u, p, t)
    x, v = u
    μ, k, F = p
    effective_stiffness = k + 0.1μ*(k + 1)
    du[1] = v
    du[2] = -effective_stiffness*x + (μ - x^2)*v - F
    nothing
end

# ╔═╡ 11013c08-a8ea-52b6-bbc9-a8c0685b4648
vocal_u0 = [0.1, 0.1]  # x(0), v(0)

# ╔═╡ 4c786b09-a539-552a-9efa-0c377d0ba2d9
vocal_p = [0.3, 1.0, 0.0]  # μ, k, F

# ╔═╡ c6dd9cbc-f9ed-5a48-b77b-8d496dc9472a
vocal_source = let
    audio_source = DESource(vocal!, copy(vocal_u0), copy(vocal_p);
        alg=Tsit5(), channel_map=[1,1])
    audio_source.data.problem = remake(audio_source.data.problem;
        abstol=1e-9, reltol=1e-7)
    audio_source
end;

# ╔═╡ 6d585394-ebb7-5bfd-a8e2-62ba32fd4b43
md"""
Drive μ $(@bind vocal_drive Slider(-0.2:0.01:0.8; default=0.3, show_value=true))

Stiffness k $(@bind vocal_stiffness Slider(0.2:0.05:2.0; default=1.0, show_value=true))

Bias F $(@bind vocal_bias Slider(-0.3:0.01:0.3; default=0.0, show_value=true))

Pitch scale (Hz) $(@bind vocal_f0 Slider(80.0:5.0:320.0; default=180.0, show_value=true))

Gain $(@bind vocal_gain Slider(0.0:0.01:0.2; default=0.05, show_value=true))
"""

# ╔═╡ a94ced60-7e23-5fb0-8ade-eec64e19406a
vocal_controls = begin
    set_param!(vocal_source, 1, Float64(vocal_drive))
    set_param!(vocal_source, 2, Float64(vocal_stiffness))
    set_param!(vocal_source, 3, Float64(vocal_bias))
    set_ts!(vocal_source, Float64(2pi*vocal_f0))
    set_gain!(vocal_source, Float64(vocal_gain))
    nothing
end

# ╔═╡ c1a39015-9e36-5c93-8ec4-1a642b4e685b
vocal_preview = let
    vocal_controls
    vocal_prob = remake(vocal_source.data.problem; u0=copy(vocal_u0), p=copy(get_params(vocal_source)), tspan=(0.0, 100.0))
    solve(vocal_prob, Tsit5(); saveat=0.02)
end;

# ╔═╡ 1eaf43bd-23b7-5963-aa74-4474ee30b817
plot(vocal_preview; idxs=[1, 2], xlabel="simulation time", ylabel="state",
    margin=5 * Plots.mm)

# ╔═╡ 70269b4d-9556-5244-8ef0-0f570a334e6e
md"""
Play $(@bind vocal_playing CheckBox(default=false))
"""

# ╔═╡ bee2ac94-fc5f-5bbd-9502-98d2c3a141fa
@use_effect([vocal_source, vocal_playing]) do
    if vocal_playing
        device = get_default_output_device()
        device < 0 && error("No output device is available. Select an audio device with list_devices().")
        reset_state!(vocal_source)

        start_DESource(vocal_source, device; buffer_size=UInt32(2048))
    end
    return () -> begin
        isactive(vocal_source) && stop_DESource(vocal_source)
    end
end

# ╔═╡ 72427869-673d-5a26-bd8a-f4ef9e610de1
md"""
Uncheck Play to stop. Stop and play again to restart from the initial condition.
Live output uses the selected gain; closing the notebook in Pluto stops its source.
"""

# ╔═╡ f8149ecc-541e-5b70-b048-7e8ea67945b8
md"""
## Birdsong: syrinx and delayed vocal tract

The syringeal labia are represented by a nonlinear oscillator with displacement
$x$ and velocity $y$. The controls $\alpha$ and $\beta$ change the pressure-like
drive and stiffness, while $\gamma$ sets the source time scale. The source
feeds a travelling pressure wave into a tract of length $L$; reflections return
after the round-trip time $T=2L/c$, with assumed sound speed $c=343$ m/s.
A linear acoustic circuit then filters the transmitted pressure.

This preserves both the source and the delayed tract from `DNL_Birdsong.jl` in
`NonlinearDynamicsPluto`. In this formulation the source drives the tract in
one direction: acoustic pressure does not feed back into the labial equation.
Time $t$ is in seconds.
"""

# ╔═╡ 9750f0c0-940c-5a7d-8649-70ec5cc058d8
md"""
The source and the two travelling pressures satisfy

```math
\begin{aligned}
\dot{x}&=y,\\
\dot{y}&=\gamma^2(-\alpha-\beta x-x^3+x^2)-\gamma x(x+1)y,\\
P_i(t)&=\nu y(t)-rP_i(t-T),\\
P_0(t)&=(1-r)P_i(t-T/2).
\end{aligned}
```

Here $\nu$ converts source velocity to pressure, $r$ is the reflection
coefficient, and $P_0$ drives the acoustic circuit. Its coordinates are
$i_1$, $\Omega=\dot{i}_1$, and $i_3$; $i_3$ is the output-flow proxy.
Using the original circuit parameters, define

```math
A=\frac{1}{C_hM_G},\quad
B=R_h\left(\frac{1}{M_B}+\frac{1}{M_G}\right),\quad
C=A+\frac{R_hR_B}{M_GM_B},\quad
D=\frac{R_hR_B}{M_GM_B}.
```

The circuit equations are

```math
\begin{aligned}
\dot{i}_1&=\Omega,\\
\dot{\Omega}&=-Ai_1-B\Omega+Ci_3+DP_0+\frac{\dot{P}_0}{M_G},\\
\dot{i}_3&=-\frac{M_G}{M_B}\Omega-\frac{R_h}{M_B}i_3+\frac{P_0}{M_B}.
\end{aligned}
```

$C_h$ is the compliance, $M_G,M_B$ are inertance parameters, and $R_h,R_B$
set the circuit's resistive terms.
"""

# ╔═╡ bee546d4-ffd4-58ab-a516-0d597a1daaae
md"""
To avoid differentiating a delayed algebraic pressure, set
$w=\Omega-P_0/M_G$. The same circuit becomes

```math
\begin{aligned}
\dot{i}_1&=w+P_0/M_G,\\
\dot{w}&=-Ai_1-Bw+Ci_3+(D-B/M_G)P_0,\\
\dot{i}_3&=-(M_G/M_B)w-(R_h/M_B)i_3.
\end{aligned}
```

The state is `u = [x, y, i1, w, i3, Pi, P0]`. The last two equations are
algebraic constraints, so their mass-matrix entries are zero. The other five
entries are one: `M*u′ = f(u, history, p, t)`.
This keeps the original pressure relations and both propagation delays.
"""

# ╔═╡ d0df2fb0-5b32-5c10-aa25-aa5d9c7877dc
function bird!(du, u, history, p, t)
    x, y, i1, w, i3, Pi, P0 = u
    α, β, γ, Ch, MG, MB, RB, Rh, r, ν, T = p
    A = 1/(Ch*MG)
    B = Rh*(1/MB + 1/MG)
    D = Rh*RB/(MG*MB)
    C = A + D
    du[1] = y
    du[2] = γ^2*(-α - β*x - x^3 + x^2) - γ*x*(x+1)*y
    du[3] = w + P0/MG
    du[4] = -A*i1 - B*w + C*i3 + (D - B/MG)*P0
    du[5] = -(MG/MB)*w - (Rh/MB)*i3
    du[6] = ν*y - r*history(p, t-T; idxs=6) - Pi
    du[7] = (1-r)*history(p, t-T/2; idxs=6) - P0
    nothing
end

# ╔═╡ cabb3097-9b2d-52d4-8773-2e65f9e2a7f9
bird_u0 = [0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

# ╔═╡ cb767325-70c0-5168-bf33-11bd3a81e043
bird_p = [0.1, 0.1, 23500.0, 1.43e-10, 0.001, 1e4, 5e6, 2.4e4, 0.65, 0.001, 2*0.025/343]
# α, β, γ, Ch, MG, MB, RB, Rh, r, ν, T

# ╔═╡ 95b6d13d-3ec0-53d4-ad6d-416ae36c2cd0
bird_mass_matrix = Diagonal([1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0])

# ╔═╡ e9bb27ee-448e-5444-b4ed-79be8d15673a
bird_system = DDEFunction(bird!; mass_matrix=bird_mass_matrix)

# ╔═╡ dafff062-3d82-5dc5-80b0-b9ff859a8694
bird_lags = [bird_p[11]/2, bird_p[11]]  # One-way and round-trip times

# ╔═╡ a37b752b-3782-5f1b-a937-a789ba539219
md"""
The live version retains the full delayed tract. A small adapter stores the
solver's dense history between audio buffers; restarting clears both the state
and that history. The first Play compiles two silent buffers before starting.
The model already uses physical seconds, so playback uses
`set_ts!(bird_source, 1.0)`. The gain multiplies `10⁷ i₃` because this circuit's flow
coordinate is small.
"""

# ╔═╡ 401a5ffd-7936-5455-beea-e124cf381251
bird_source, bird_delay_history = buffered_dde_source(bird_system, bird_u0, bird_p;
    constant_lags=bird_lags, channel_map=[5,5]);

# ╔═╡ c064d89e-e707-57f1-968f-1e12e8aa907e
md"""
Pressure control α $(@bind bird_alpha Slider(0.0:0.01:0.2; default=0.1, show_value=true))

Stiffness control β $(@bind bird_beta Slider(0.0:0.01:0.3; default=0.1, show_value=true))

Source rate γ (s⁻¹) $(@bind bird_rate Slider(10000.0:500.0:35000.0; default=23500.0, show_value=true))

Gain $(@bind bird_gain Slider(0.0:0.01:0.2; default=0.05, show_value=true))

Change the source controls while it runs. Edit the fixed tract parameters
above with playback stopped to rebuild the source with new delays.
"""

# ╔═╡ 2c6f9a8c-1077-599d-9695-d0572efba672
bird_controls = begin
    set_param!(bird_source, 1, Float64(bird_alpha))
    set_param!(bird_source, 2, Float64(bird_beta))
    set_param!(bird_source, 3, Float64(bird_rate))
    set_ts!(bird_source, 1.0)
    set_gain!(bird_source, 1e7*Float64(bird_gain))
    nothing
end

# ╔═╡ 879b9afe-0433-5510-afa2-2155b8e8856e
bird_preview = let
    bird_controls
    bird_history(bird_p,t;idxs=nothing) = isnothing(idxs) ? zeros(7) : zeros(7)[idxs]
    bird_prob = DDEProblem(bird_system, copy(bird_u0), bird_history, (0.0,0.05), copy(get_params(bird_source));
        constant_lags=bird_lags)
    solve(bird_prob, MethodOfSteps(Rodas5P()); abstol=1e-10, reltol=1e-6, dt=1e-6)
end;

# ╔═╡ ee1384d0-f864-5ef6-b4a9-f499fc00bc06
let
    times = range(max(0.0, last(bird_preview.t)-0.01), last(bird_preview.t); length=2000)
    values = bird_preview(times)
    source_plot = plot(times, values[1,:]; xlabel="time (s)", ylabel="labial x", legend=false)
    output_plot = plot(times, values[5,:]; xlabel="time (s)", ylabel="output flow i₃", legend=false)
    orbit = plot(values[1,:], values[2,:] ./ bird_rate;
        xlabel="x", ylabel="y/γ", legend=false)
    pressure_plot = plot(times, [values[6,:] values[7,:]];
        xlabel="time (s)", ylabel="pressure", label=["Pᵢ" "P₀"])
    plot(source_plot, output_plot, orbit, pressure_plot;
        layout=(2,2), size=(850,550), margin=5 * Plots.mm)
end

# ╔═╡ 9cc27184-7009-5af3-9dc7-ecc0bf2ac0cd
md"""
Play $(@bind bird_playing CheckBox(default=false))
"""

# ╔═╡ a973ea1b-e137-5390-8a57-c26ccc8cf676
@use_effect([bird_source, bird_playing]) do
    if bird_playing
        device = get_default_output_device()
        device < 0 && error("No output device is available. Select an audio device with list_devices().")
        warmup_delay_source!(bird_source, bird_delay_history)
        start_DESource(bird_source, device; sample_rate=48000.0, buffer_size=UInt32(2048))
    end
    return () -> begin
        isactive(bird_source) && stop_DESource(bird_source)
    end
end

# ╔═╡ e2e62851-2558-53d6-9e17-3354477a8b3d
md"""
Uncheck Play to stop. Stop and play again to restart from the initial condition.
Live output uses the selected gain; closing the notebook in Pluto stops its source.
"""

# ╔═╡ Cell order:
# ╠═64a40c39-28e8-5115-977e-d4c329c5badc
# ╟─20db5510-3672-55d0-b7cb-e35c23e731f2
# ╟─4de90467-b82e-5c5b-9ea7-966ab8a026f2
# ╟─4de0c122-664f-5f68-aae3-8097e02f400e
# ╟─e02194ae-76a5-5b91-92d4-a19fc9878d6f
# ╠═061cf701-efa1-5a0b-8db4-b6b7464c9c9d
# ╠═11013c08-a8ea-52b6-bbc9-a8c0685b4648
# ╠═4c786b09-a539-552a-9efa-0c377d0ba2d9
# ╠═c6dd9cbc-f9ed-5a48-b77b-8d496dc9472a
# ╟─6d585394-ebb7-5bfd-a8e2-62ba32fd4b43
# ╠═a94ced60-7e23-5fb0-8ade-eec64e19406a
# ╠═c1a39015-9e36-5c93-8ec4-1a642b4e685b
# ╠═1eaf43bd-23b7-5963-aa74-4474ee30b817
# ╟─70269b4d-9556-5244-8ef0-0f570a334e6e
# ╠═bee2ac94-fc5f-5bbd-9502-98d2c3a141fa
# ╟─72427869-673d-5a26-bd8a-f4ef9e610de1
# ╟─f8149ecc-541e-5b70-b048-7e8ea67945b8
# ╟─9750f0c0-940c-5a7d-8649-70ec5cc058d8
# ╟─bee546d4-ffd4-58ab-a516-0d597a1daaae
# ╠═d0df2fb0-5b32-5c10-aa25-aa5d9c7877dc
# ╠═cabb3097-9b2d-52d4-8773-2e65f9e2a7f9
# ╠═cb767325-70c0-5168-bf33-11bd3a81e043
# ╠═95b6d13d-3ec0-53d4-ad6d-416ae36c2cd0
# ╠═e9bb27ee-448e-5444-b4ed-79be8d15673a
# ╠═dafff062-3d82-5dc5-80b0-b9ff859a8694
# ╟─a37b752b-3782-5f1b-a937-a789ba539219
# ╠═401a5ffd-7936-5455-beea-e124cf381251
# ╟─c064d89e-e707-57f1-968f-1e12e8aa907e
# ╠═2c6f9a8c-1077-599d-9695-d0572efba672
# ╠═879b9afe-0433-5510-afa2-2155b8e8856e
# ╠═ee1384d0-f864-5ef6-b4a9-f499fc00bc06
# ╟─9cc27184-7009-5af3-9dc7-ecc0bf2ac0cd
# ╠═a973ea1b-e137-5390-8a57-c26ccc8cf676
# ╟─e2e62851-2558-53d6-9e17-3354477a8b3d
