# Run from the repository: julia --project=. scripts/render_audio.jl [output_folder]
# Uses the actual ODE notebooks, including their sample rate and export settings.
ENV["GKSwstype"] = "100"
empty!(LOAD_PATH)
append!(LOAD_PATH, ["@", "@stdlib"])
using WAV

const ROOT = dirname(@__DIR__)
destination = isempty(ARGS) ? joinpath(ROOT, "audio") : abspath(only(ARGS))
mkpath(destination)

for (notebook, filename) in (
        ("NLD13_Reed_Resonator.jl", "reed_resonator.wav"),
        ("NLD14_Bowed_Resonator.jl", "bowed_resonator.wav"),
        ("NLD15_Struck_Resonator.jl", "struck_resonator.wav"))
    lesson = Module(gensym(:AudioLesson))
    Base.include(lesson, joinpath(ROOT, "src", "Pluto", "03_ODE", notebook))
    Base.invokelatest() do
        path = joinpath(destination, filename)
        write(path, lesson.wav_data)
        samples, rate = wavread(path)
        @assert rate == lesson.sample_rate
        @assert size(samples, 1) == round(Int, lesson.duration * rate)
        @assert all(isfinite, samples) && 0 < maximum(abs, samples) <= 0.701
        println(filename, ": ", rate, " Hz, ", size(samples, 1) / rate, " s, mono PCM16")
    end
end
