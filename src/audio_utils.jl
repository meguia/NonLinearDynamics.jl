using DSP: resample

"""
    prepare_audio(samples, sample_rate; oversample=4, peak=0.7, fade_seconds=0.01)

Low-pass and downsample an oversampled mono signal, remove its DC offset,
normalize its peak, and apply short cosine fades for WAV playback.
The input length must be a multiple of `oversample`; output has exactly
`length(samples) ÷ oversample` samples. Normalization changes loudness, so use
the original solution when comparing model amplitudes.
"""
function prepare_audio(samples::AbstractVector, sample_rate::Integer;
        oversample::Integer=4, peak::Real=0.7, fade_seconds::Real=0.01)
    sample_rate > 0 || throw(ArgumentError("sample_rate must be positive"))
    oversample > 0 || throw(ArgumentError("oversample must be positive"))
    0 < peak <= 1 || throw(ArgumentError("peak must be in (0, 1]"))
    isfinite(fade_seconds) && fade_seconds >= 0 || throw(ArgumentError("invalid fade duration"))
    isempty(samples) && throw(ArgumentError("samples must not be empty"))
    length(samples) % oversample == 0 || throw(ArgumentError("incomplete oversampling frame"))
    all(isfinite, samples) || throw(ArgumentError("audio contains non-finite values"))
    n = length(samples) ÷ oversample
    centered = Float64.(samples)
    centered .-= sum(centered) / length(centered)
    # Remove DC before filtering so a constant state does not excite filter edges.
    maximum(abs, centered) <= 1e-12 && return zeros(n)
    signal = oversample == 1 ? centered : resample(centered, 1 // oversample)[1:n]
    signal .-= sum(signal) / n
    amplitude = maximum(abs, signal)
    amplitude <= 1e-12 && return zeros(n)
    signal .*= peak / amplitude
    fade = min(round(Int, fade_seconds * sample_rate), n ÷ 2)
    if fade >= 2
        ramp = (1 .- cos.(range(0, pi; length=fade))) ./ 2
        signal[1:fade] .*= ramp
        signal[end-fade+1:end] .*= reverse(ramp)
    end
    signal
end
