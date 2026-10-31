using LinearAlgebra
using Random
using Statistics
using Test
using DifferentialGamesBase
using DifferentialGamesBaseSolvers

function _two_player_lq(; N = 6)
    A = [1.0 0.1; -0.05 0.97]
    B = [reshape([0.0, 0.1], 2, 1), reshape([0.1, 0.0], 2, 1)]
    Q = [[2.0 0.2; 0.2 1.0], [1.0 0.0; 0.0 3.0]]
    R = [reshape([0.5], 1, 1), reshape([0.8], 1, 1)]
    Qf = [Matrix(3.0I, 2, 2), Matrix(2.0I, 2, 2)]
    LQGameProblem(A, B, Q, R, Qf, [1.0, -0.5], Float64(N); dt = 1.0)
end

@testset "SFNELQ: moments formula and shared observations" begin
    g = _two_player_lq()
    W = [0.03 0.01; 0.01 0.02]
    Σ0 = [0.05 0.0; 0.0 0.02]

    @testset "moment-based expected cost reproduces the perfect-information value" begin
        sg = StochasticGameProblem(g, ProcessNoiseModel(W); initial_covariance = Σ0)
        sol = solve(sg, SFNELQ())
        det = solve(g, FNELQ())
        ec = expected_costs_from_moments(sg, det, sol.solver_info[:state_covariances],
            sol.solver_info[:control_covariances], sol.solver_info[:state_control_covariances])
        @test ec ≈ sol.solver_info[:expected_costs] rtol = 1e-10
    end

    @testset "near-perfect measurements recover perfect-information moments" begin
        perfect = solve(StochasticGameProblem(g, ProcessNoiseModel(W); initial_covariance = Σ0), SFNELQ())
        obs = linear_observation(Matrix(1.0I, 2, 2), Matrix(1e-12I, 2, 2))
        sh = solve(StochasticGameProblem(g, ProcessNoiseModel(W); initial_covariance = Σ0, information = obs), SFNELQ())
        @test sh.solver_info[:separation_theorem_used]
        @test sh.solver_info[:expected_costs] ≈ perfect.solver_info[:expected_costs] rtol = 1e-3
        @test sh.solver_info[:state_covariances][end] ≈ perfect.solver_info[:state_covariances][end] rtol = 1e-3
        @test norm(sh.solver_info[:posterior_covariances][end]) < 1e-9
    end

    @testset "partial noisy measurements: filter + policy vs Monte Carlo" begin
        H = [1.0 0.0]; V = reshape([0.05], 1, 1)
        obs = linear_observation(H, V)
        sg = StochasticGameProblem(g, ProcessNoiseModel(W); initial_covariance = Σ0, information = obs)
        sol = solve(sg, SFNELQ())
        info = sol.solver_info
        N = length(info[:control_covariances])
        # the physical-state covariance is NOT the posterior estimation covariance
        @test !isapprox(info[:state_covariances][end], info[:posterior_covariances][end]; rtol = 0.1)
        @test all(k -> isposdef(Symmetric(info[:state_covariances][k] + 1e-12I)), 1:N+1)

        dyn = g.dynamics
        P = [vcat((sol.strategy.gains[i][k] for i in 1:2)...) for k in 1:N]
        α = [vcat((sol.strategy.feedforward[i][k] for i in 1:2)...) for k in 1:N]
        A = get_A(dyn, 1); Bj = get_B_concatenated(dyn, 1)
        rng = MersenneTwister(21); nmc = 60_000
        Lw = cholesky(W).L; L0 = cholesky(Σ0).L
        costs = zeros(2, nmc); xend = zeros(2, nmc); x̂end = zeros(2, nmc)
        objs = [get_objective(g, i) for i in 1:2]
        ranges = [1:1, 2:2]
        for s in 1:nmc
            x = g.initial_state + L0 * randn(rng, 2); c = zeros(2)
            y1 = H * x + sqrt(V[1]) * randn(rng, 1)
            x̂ = g.initial_state + info[:initial_kalman_gain] * (y1 - H * g.initial_state)
            for k in 1:N
                u = -P[k] * x̂ - α[k]
                for i in 1:2
                    c[i] += evaluate_stage_cost(objs[i].stage_cost, x, u[ranges[i]], nothing, k)
                end
                x = A * x + Bj * u + Lw * randn(rng, 2)
                y = H * x + sqrt(V[1]) * randn(rng, 1)
                x̂⁻ = A * x̂ + Bj * u
                x̂ = x̂⁻ + info[:kalman_gains][k] * (y - H * x̂⁻)
            end
            for i in 1:2
                c[i] += evaluate_terminal_cost(objs[i].terminal_cost, x, nothing)
            end
            costs[:, s] = c; xend[:, s] = x; x̂end[:, s] = x̂
        end
        for i in 1:2
            se = Statistics.std(costs[i, :]) / sqrt(nmc)
            @test abs(Statistics.mean(costs[i, :]) - info[:expected_costs][i]) < 4 * se
        end
        @test isapprox(info[:state_covariances][end], Statistics.cov(xend'); atol = 5e-3, rtol = 0.05)
        @test isapprox(info[:estimate_covariances][end], Statistics.cov(x̂end'); atol = 5e-3, rtol = 0.05)
        @test isapprox(info[:posterior_covariances][end], Statistics.cov((xend - x̂end)'); atol = 5e-3, rtol = 0.05)
    end

    @testset "capability errors for unsupported information" begin
        nonlinear = SharedObservation((x, v, p, t) -> x[1:1] .^ 2 .+ v, 1, reshape([0.1], 1, 1))
        @test_throws ErrorException solve(StochasticGameProblem(g, ProcessNoiseModel(W); information = nonlinear), SFNELQ())
        @test_throws ErrorException solve(StochasticGameProblem(g, ProcessNoiseModel(W); information = OpenLoopInformation()), SFNELQ())
    end
end
