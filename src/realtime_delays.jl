module RealtimeDelays

using DifferentialEquations
using RealTimeAudioDiffEq
const PA = RealTimeAudioDiffEq.LibPortAudio

function segment_value(segment, t, idxs)
    first(segment.t)-1e-12 <= t <= last(segment.t)+1e-12 || return nothing
    segment(clamp(t, first(segment.t), last(segment.t)); idxs)
end

"""
    buffered_dde_source(f, u0, p; constant_lags, channel_map, ...)

Construct a live delay source without opening an audio device. Each completed
audio buffer retains the accepted-step history needed by the next buffer.
Returns (source, segments); empty segments when resetting the source.

This adapter uses the problem constructor in the pinned RealTimeAudioDiffEq 0.1
and the internal dense history of DelayDiffEq 5. Its callback still produces the
regular sample grid expected by RealTimeAudioDiffEq.
"""
function buffered_dde_source(f, u0, p; constant_lags, channel_map=[1,1],
        algorithm=MethodOfSteps(Rodas5P()), abstol=1e-9, reltol=1e-5)
    all(>(0), constant_lags) || throw(ArgumentError("Delays must be positive"))
    segments = Any[]
    initial_history = zeros(length(u0))
    function history(parameters, t; idxs=nothing)
        if t <= 0
            return isnothing(idxs) ? copy(initial_history) : initial_history[idxs]
        end
        for segment in Iterators.reverse(segments)
            value = segment_value(segment, t, idxs)
            if !isnothing(value)
                # The function barrier keeps interpolation specialized despite
                # storing solutions from separately initialized audio buffers.
                return idxs isa Integer ? oftype(t, value) : value
            end
        end
        error("Missing delay history at t=$t; reset the source and its history together.")
    end
    function retain_history(callback, u, t, integrator)
        t == last(integrator.sol.prob.tspan) || return
        # The public solution contains audio samples; the internal one retains
        # dense interpolation between the delay solver's accepted steps.
        push!(segments, integrator.integrator.sol)
        cutoff = t - maximum(constant_lags)
        while length(segments) > 1 && last(first(segments).t) < cutoff
            popfirst!(segments)
        end
        nothing
    end
    retain = DiscreteCallback((u,t,integrator) -> false, integrator -> nothing;
        finalize=retain_history, save_positions=(false,false))
    problem = DDEProblem(f, copy(u0), history, (0.0,1.0), copy(p);
        constant_lags=collect(constant_lags), callback=retain,
        abstol, reltol, dt=min(minimum(constant_lags)/10, 1e-6), dense=false)
    source = RealTimeAudioDiffEq._DESource(problem, algorithm, channel_map)
    source, segments
end

"""
    warmup_delay_source!(source, segments; sample_rate=48000.0, frames=2048)

Compile two consecutive callbacks into an in-memory buffer before opening the
audio device. Afterwards the state and delay history are reset for playback.
"""
function warmup_delay_source!(source, segments; sample_rate=48000.0, frames=2048)
    isactive(source) && throw(ArgumentError("Stop playback before warming up"))
    stream = source.data.stream_data
    previous = (stream.sample_rate, stream.buffer_size, stream.n_channels)
    mapping = get_channelmap(source)
    channels = mapping isa AbstractMatrix ? size(mapping,2) : length(mapping)
    stream.sample_rate, stream.buffer_size, stream.n_channels = sample_rate, frames, channels
    buffer = fill(Float32(NaN), frames*channels)
    data = Ref(source.data)
    try
        reset_state!(source)
        empty!(segments)
        GC.@preserve source buffer data begin
            for _ in 1:2
                fill!(buffer, Float32(NaN))
                before = source.data.state.t
                result = ccall(Base.unsafe_convert(Ptr{Cvoid}, source.callback), PA.PaStreamCallbackResult,
                    (Ptr{Cvoid}, Ptr{Cvoid}, Culong, Ptr{PA.PaStreamCallbackTimeInfo},
                     PA.PaStreamCallbackFlags, Ptr{Cvoid}),
                    C_NULL, pointer(buffer), frames, C_NULL, 0, data)
                result == PA.paContinue && all(isfinite,buffer) && source.data.state.t > before ||
                    error("The delay solver did not complete its audio warmup")
            end
        end
    finally
        reset_state!(source)
        empty!(segments)
        stream.sample_rate, stream.buffer_size, stream.n_channels = previous
    end
    nothing
end

end
