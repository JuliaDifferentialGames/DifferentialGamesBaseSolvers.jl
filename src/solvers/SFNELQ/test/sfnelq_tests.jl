using LinearAlgebra
using Random
using Statistics
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

    @testset "expected cost matches Monte Carlo (affine terms, singular noise)" begin
        A2 = [1.0 0.1; 0.0 1.0]; B2 = reshape([0.0, 0.1], 2, 1)
        Q2 = [1.0 0.0; 0.0 2.0]; R2 = reshape([0.5], 1, 1); Qf2 = [3.0 0.5; 0.5 1.0]
        qv = [0.3, -0.2]; rv = [0.4]; x0 = [1.0, -1.0]
        g2 = LQGameProblem(A2, [B2], [Q2], [R2], [Qf2], x0, 3.0; dt = 1.0, q = [qv], r = [rv])
        W = [0.09 0.0; 0.0 0.0]   # rank-deficient
        sg = StochasticGameProblem(g2, ProcessNoiseModel(W))
        sol = solve(sg, SFNELQ())
        P = sol.strategy.gains[1]; ff = sol.strategy.feedforward[1]
        N = length(P)
        rng = MersenneTwister(7); nmc = 100_000
        costs = zeros(nmc)
        for s in 1:nmc
            x = copy(x0); c = 0.0
            for k in 1:N
                u = -P[k] * x - ff[k]
                c += 0.5 * (x' * Q2 * x + u' * R2 * u) + qv' * x + rv' * u
                x = A2 * x + B2 * u + [0.3 * randn(rng), 0.0]
            end
            costs[s] = c + 0.5 * x' * Qf2 * x
        end
        se = Statistics.std(costs) / sqrt(nmc)
        @test abs(Statistics.mean(costs) - sol.trajectories[1].cost) < 4 * se
    end

    @testset "zero noise reproduces the deterministic cost" begin
        sg0 = StochasticGameProblem(g, ProcessNoiseModel(zeros(2, 2)))
        s0 = solve(sg0, SFNELQ()); sd = solve(g, FNELQ())
        @test s0.trajectories[1].cost ≈ sd.trajectories[1].cost
        @test all(s0.solver_info[:noise_cost] .== 0)
    end

    @testset "noise cost respects objective scaling" begin
        obj = g.objectives[1]
        scaled_obj = PlayerObjective(obj.player_id, obj.stage_cost, obj.terminal_cost, 2.5)
        gs = GameProblem{Float64}(g.n_players, [scaled_obj], g.dynamics, g.initial_state,
            g.private_constraints, g.shared_constraints, g.time_horizon, g.metadata)
        sgs = StochasticGameProblem(gs, ProcessNoiseModel(Σ))
        base = solve(sgame, SFNELQ())
        sc = solve(sgs, SFNELQ())
        @test sc.trajectories[1].cost ≈ 2.5 * base.trajectories[1].cost
    end

    @testset "rejects multiplicative-noise StochasticGameProblem" begin
        drift = LinearDynamics(A, B)
        sde = ItoSDEDynamics(drift, (x,u,p,t) -> [0.1 0.0; 0.0 0.1], 2)
        sgame_mult = StochasticGameProblem(g, sde)
        @test_throws ErrorException solve(sgame_mult, SFNELQ())
    end
end
