module NonLinearDynamics

    using Plots
    using LinearAlgebra
    using DifferentialEquations
    using ForwardDiff
    using IntervalRootFinding
    using StaticArrays 

    export  flux1D,
            potential1D,
            flux2d_vectorfield,
            flux2d_nullclines,
            flux2d_animated,
            flux2d_manifolds,
            flux2d_forced,
            classification_linear,
            phase_portrait,
            attractor_basin,
            poincare_forced,
            poincare_forced_zoom,
            recurrence_plot,
            saddle_orbit2D,
            saddle_manifolds_forced,
            butterfly

    include("DNL_utils.jl")

end
