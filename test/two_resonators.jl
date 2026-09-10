# Run: julia --project=. --startup-file=no test/two_resonators.jl
# Read the actual teaching function without rendering the notebook or opening audio.
using Test, DifferentialEquations, ForwardDiff, LinearAlgebra, Statistics

const LESSON_PATH = joinpath(@__DIR__, "..", "src", "Pluto", "03_ODE", "NLD16_Two_Resonators.jl")
const LESSON = Module(:TwoResonatorModel)
const LESSON_TEXT = read(LESSON_PATH, String)
for pattern in (r"(?ms)^function model!.*?^end", r"(?m)^p = .*", r"(?m)^u0 = .*")
    Base.include_string(LESSON, match(pattern, LESSON_TEXT).match, LESSON_PATH)
end

function check_two_resonators()
    p = copy(LESSON.p)
    u0 = copy(LESSON.u0)
    model! = LESSON.model!
    function section_solution(parameters, initial; stop=10000.0, discard=5000.0, tolerance=1e-8, algorithm=Tsit5())
        points = Vector{Float64}[]
        times = Float64[]
        # A fixed surface and direction, using continuous root finding rather than audio samples.
        condition(u,t,integrator) = u[3]-u0[3]
        function record!(integrator)
            if integrator.t >= discard
                push!(points,copy(integrator.u))
                push!(times,integrator.t)
            end
            u_modified!(integrator,false)
        end
        callback = ContinuousCallback(condition,record!,nothing;save_positions=(false,false))
        sol = solve(ODEProblem(model!,initial,(0.0,stop),parameters),algorithm;
            callback,abstol=tolerance/100,reltol=tolerance,maxiters=10^7,save_everystep=false)
        @test sol.retcode == ReturnCode.Success
        @test !isempty(points)
        (sol=sol,points=points,times=times)
    end
    spread(run,k) = std(getindex.(run.points,k))

    @testset "Sustained two-resonator quasiperiodicity" begin
        torus = section_solution(p,u0;stop=30000.0)
        @test length(torus.times) > 3500
        early = findall(<(10000.0),torus.times)
        late = findall(>=(25000.0),torus.times)
        for k in (5,6)
            first_spread = std(getindex.(torus.points[early],k))
            last_spread = std(getindex.(torus.points[late],k))
            @test last_spread > 0.005
            @test last_spread ≈ first_spread rtol=0.02
        end

        periodic_parameters = copy(p)
        periodic_parameters[1] = 0.45
        periodic = section_solution(periodic_parameters,u0)
        @test spread(periodic,5) < 1e-5
        @test spread(periodic,6) < 1e-5

        # Changing the blowing preset from the periodic state must also reach the torus.
        switched = section_solution(p,last(periodic.sol.u))
        @test spread(switched,5) > 0.005
        @test spread(switched,6) > 0.005

        # Check a different integration method and a smaller airflow smoothing scale.
        refined = section_solution(p,u0;algorithm=Vern7(),tolerance=1e-9)
        smoother_parameters = copy(p)
        smoother_parameters[end] /= 10
        smaller_smoothing = section_solution(smoother_parameters,u0)
        for run in (refined,smaller_smoothing), k in (5,6)
            @test spread(run,k) ≈ spread(torus,k) rtol=0.02
        end

        # Three tangent directions distinguish a stable torus from a limit cycle or chaos.
        function tangent!(dY,Y,parameters,t)
            model!(view(dY,1:6),view(Y,1:6),parameters,t)
            J = ForwardDiff.jacobian(view(Y,1:6)) do u
                du = similar(u)
                model!(du,u,parameters,t)
                du
            end
            mul!(reshape(view(dY,7:24),6,3),J,reshape(view(Y,7:24),6,3))
        end
        Y = vcat(last(torus.sol.u),vec(Matrix{Float64}(I,6,6)[:,1:3]))
        integrator = init(ODEProblem(tangent!,Y,(0.0,21000.0),p),Tsit5();
            abstol=1e-10,reltol=1e-8,maxiters=10^8,save_everystep=false)
        growth = zeros(3)
        for j in 1:4200
            step!(integrator,5.0,true)
            # Abort on a failed integration; elapsed requested time cannot replace actual time.
            isapprox(integrator.t,5j;atol=1e-7) || error("Tangent integration stopped early")
            factor = qr(reshape(view(integrator.u,7:24),6,3))
            j > 200 && (growth .+= log.(abs.(diag(factor.R))))
            integrator.u[7:24] .= vec(Matrix(factor.Q))
            u_modified!(integrator,true)
        end
        exponents = growth/20000.0
        println("Largest Lyapunov exponents per dimensionless time: ",exponents)
        @test all(abs.(exponents[1:2]) .< 5e-4)
        @test exponents[3] < -0.005
    end

    @testset "Two-resonator control endpoints" begin
        for blowing in (0.45,0.47), ratio in (2.0,2.6)
            parameters = copy(p)
            parameters[1],parameters[5] = blowing,ratio
            sol = solve(ODEProblem(model!,u0,(0.0,6000.0),parameters),Tsit5();
                abstol=1e-9,reltol=1e-7,maxiters=10^7,saveat=1.0)
            @test sol.retcode == ReturnCode.Success
            @test all(u -> all(isfinite,u),sol.u)
            @test maximum(norm,sol.u) < 10
        end
    end
end

Base.invokelatest(check_two_resonators)
