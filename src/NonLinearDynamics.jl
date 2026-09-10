module NonLinearDynamics

    using Plots
    using DifferentialEquations
    using ForwardDiff
    using IntervalRootFinding
    using StaticArrays 

    export  flow1D,
            potential1D,
            flow2d_vectorfield,
            flow2d_nullclines,
            flow2d_animated,
            flow2d_manifolds,
            flow2d_forced,
            classification_linear,
            phase_portrait,
            attractor_basin,
            poincare_forced,
            poincare_forced!,
            poincare_forced_zoom,
            saddle_orbit2D,
            saddle_manifolds_forced,
            butterfly,
            flow3d,
            prepare_audio

    include("NLD_utils.jl")
    include("audio_utils.jl")

end
