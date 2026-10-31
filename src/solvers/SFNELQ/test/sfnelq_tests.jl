using LinearAlgebra
using Test
using DifferentialGamesBase
using DifferentialGamesBaseSolvers

@testset "SFNELQ" begin
    A = [1.0 1.0; 0.0 1.0]
    B = [reshape([0.0, 1.0], 2, 1)]
    Q = diagm([1.0, 1.0]); R = reshape([0.1], 1, 1); Qf = diagm([10.0, 10.0])
    g = LQGameProblem(A, B, [Q], [R], [Qf], [1.0, 0.0], 2.0; dt = 1.0)

    Σ = diagm([0.01, 0.01])
    sgame = StochasticGameProblem(g, ProcessNoiseModel(Σ))

    @testset "certainty equivalence: gains match deterministic FNELQ" begin
        det_sol = solve(g, FNELQ())
        sto_sol = solve(sgame, SFNELQ())

        @test sto_sol.solver_info[:certainty_equivalent]
        # Feedback gains (strategy.gains) must be identical under additive-only noise.
        for i in 1:1, k in 1:length(det_sol.strategy.gains[i])
            @test det_sol.strategy.gains[i][k] ≈ sto_sol.strategy.gains[i][k]
        end
        # Mean trajectory is identical too (noise is zero-mean).
        @test det_sol.state_trajectory ≈ sto_sol.state_trajectory
    end

    @testset "expected cost increases by the noise correction" begin
        det_sol = solve(g, FNELQ())
        sto_sol = solve(sgame, SFNELQ())
        @test all(sto_sol.solver_info[:noise_cost] .>= 0)
        for i in 1:1
            @test sto_sol.trajectories[i].cost ≈
                  det_sol.trajectories[i].cost + sto_sol.solver_info[:noise_cost][i]
        end
    end

    @testset "rejects multiplicative-noise StochasticGameProblem" begin
        drift = LinearDynamics(A, B)
        sde = ItoSDEDynamics(drift, (x,u,p,t) -> [0.1 0.0; 0.0 0.1], 2)
        sgame_mult = StochasticGameProblem(g, sde)
        @test_throws ErrorException solve(sgame_mult, SFNELQ())
    end
end
