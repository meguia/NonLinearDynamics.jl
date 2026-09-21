module NonLinearDynamics

export flow1D,
       flow2d_vectorfield,
       flow2d_nullclines,
       classification_linear,
       flow2d_manifolds,
       phase_portrait,
       attractor_basin,
       flow2d_forced,
       poincare_forced,
       poincare_forced!,
       poincare_forced_zoom,
       saddle_orbit2D,
       saddle_manifolds_forced,
       prepare_audio

include("NLD_utils.jl")
include("audio_utils.jl")

end
