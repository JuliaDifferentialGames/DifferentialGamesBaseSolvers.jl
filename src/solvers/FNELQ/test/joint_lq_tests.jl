using LinearAlgebra
using Random
using Statistics
using Test
using DifferentialGamesBase
using DifferentialGamesBaseSolvers

# Unilateral-deviation check: a feedback Nash equilibrium is Nash in the space of
# affine feedback laws u^i_k = -P^i_k x_k - α^i_k, so every player's cost must be
# stationary w.r.t. its own (P^i, α^i) with the other players' laws held fixed.
@testset "FNELQ: joint-control costs are Nash (unilateral deviation)" begin
    rng = MersenneTwister(3)
    n, N = 2, 3
    A = [1.0 0.1; -0.05 0.95]
    B = [reshape([0.0, 0.1], 2, 1), reshape([0.1, 0.05], 2, 1)]
    m = 2
    sym(X) = (X + X') / 2
    Qs = [sym(randn(rng, n, n)) + 4I for _ in 1:2] .|> Matrix
    Rs = [[2.0 0.3; 0.3 0.4], [0.5 -0.2; -0.2 1.5]]   # full 2×2 joint-control weights
    Ms = [0.2 * randn(rng, n, m) for _ in 1:2]
    qs = [randn(rng, n) for _ in 1:2]
    rs = [randn(rng, m) for _ in 1:2]
    Qfs = [Matrix(sym(randn(rng, n, n)) + 4I) for _ in 1:2]
    qfs = [randn(rng, n) for _ in 1:2]
    x0 = [1.0, -0.5]
    g = JointLQGameProblem(A, B, Qs, Rs, Qfs, x0, Float64(N); dt = 1.0, M = Ms, q = qs, r = rs, qf = qfs)
    sol = solve(g, FNELQ())
    @test sol.solver_info[:second_order_ok]

    P = [copy(sol.strategy.gains[i]) for i in 1:2]
    α = [copy(sol.strategy.feedforward[i]) for i in 1:2]

    function cost(i, P, α)
        x = copy(x0); c = 0.0
        for k in 1:N
            u = [only(-P[j][k] * x - α[j][k]) for j in 1:2]
            c += 0.5 * (x' * Qs[i] * x + u' * Rs[i] * u) + x' * Ms[i] * u + qs[i]' * x + rs[i]' * u
            x = A * x + vec(B[1]) * u[1] + vec(B[2]) * u[2]
        end
        return c + 0.5 * x' * Qfs[i] * x + qfs[i]' * x
    end

    @test cost(1, P, α) ≈ sol.trajectories[1].cost
    @test cost(2, P, α) ≈ sol.trajectories[2].cost

    h = 1e-5
    for i in 1:2, k in 1:N
        for idx in eachindex(P[i][k])
            Pp = deepcopy(P); Pm = deepcopy(P)
            Pp[i][k][idx] += h; Pm[i][k][idx] -= h
            @test abs((cost(i, Pp, α) - cost(i, Pm, α)) / 2h) < 1e-6
        end
        αp = deepcopy(α); αm = deepcopy(α)
        αp[i][k][1] += h; αm[i][k][1] -= h
        @test abs((cost(i, P, αp) - cost(i, P, αm)) / 2h) < 1e-6
    end
end

@testset "FNELQ: own-control cross term M together with affine terms" begin
    # Regression: the linear value coefficient needs a -M α term when M, q and r are all active.
    A = [1.0 0.1; 0.0 1.0]; B = [reshape([0.0, 0.1], 2, 1)]
    Q = [1.0 0.0; 0.0 2.0]; R = reshape([0.5], 1, 1); Qf = [3.0 0.5; 0.5 1.0]
    M = reshape([0.4, -0.3], 2, 1); q = [0.3, -0.2]; r = [0.4]
    x0 = [1.0, -1.0]
    g = LQGameProblem(A, [B[1]], [Q], [R], [Qf], x0, 3.0; dt = 1.0, M = [M], q = [q], r = [r])
    sol = solve(g, FNELQ())
    u = vec(sol.trajectories[1].controls)
    function J(u)
        x = copy(x0); c = 0.0
        for k in 1:3
            c += 0.5 * (x' * Q * x + R[1] * u[k]^2) + x' * M[:, 1] * u[k] + q' * x + r[1] * u[k]
            x = A * x + vec(B[1]) * u[k]
        end
        return c + 0.5 * x' * Qf * x
    end
    h = 1e-5
    for k in 1:3
        e = zeros(3); e[k] = h
        @test abs((J(u + e) - J(u - e)) / 2h) < 1e-6   # single player: feedback = open-loop optimum
    end
end

@testset "SFNELQ: LTV noise and uncertain initial state" begin
    A = [1.0 0.1; 0.0 1.0]; B = [reshape([0.0, 0.1], 2, 1)]
    Q = [1.0 0.0; 0.0 2.0]; R = reshape([0.5], 1, 1); Qf = [3.0 0.5; 0.5 1.0]
    x0 = [1.0, -1.0]
    g = LQGameProblem(A, B, [Q], [R], [Qf], x0, 3.0; dt = 1.0)
    Ws = [[0.04 0.0; 0.0 0.0], [0.0 0.0; 0.0 0.09], zeros(2, 2)]
    Σ0 = [0.05 0.01; 0.01 0.02]
    sg = StochasticGameProblem(g, ProcessNoiseModel(Ws); initial_covariance = Σ0)
    sol = solve(sg, SFNELQ())
    P = sol.strategy.gains[1]; ff = sol.strategy.feedforward[1]

    rng = MersenneTwister(11); nmc = 100_000
    L0 = cholesky(Σ0).L
    Ls = [[0.2, 0.0], [0.0, 0.3], [0.0, 0.0]]
    costs = zeros(nmc); xT = zeros(2, nmc)
    for s in 1:nmc
        x = x0 + L0 * randn(rng, 2); c = 0.0
        for k in 1:3
            u = -P[k] * x - ff[k]
            c += 0.5 * (x' * Q * x + u' * R * u)
            x = A * x + vec(B[1]) * only(u) + Ls[k] .* randn(rng, 2)
        end
        costs[s] = c + 0.5 * x' * Qf * x
        xT[:, s] = x
    end
    se = Statistics.std(costs) / sqrt(nmc)
    @test abs(Statistics.mean(costs) - sol.trajectories[1].cost) < 4 * se

    # Analytic covariance trajectory vs sample covariance of the terminal state.
    Σ_end = get_covariance_trajectory(sol)[end]
    @test isapprox(Σ_end, Statistics.cov(xT'); atol = 0.01, rtol = 0.05)
    @test get_covariance_trajectory(sol)[1] ≈ Σ0
    # Cov(u) = P Σ Pᵀ and Cov(x,u) = -Σ Pᵀ at k = 1
    @test sol.solver_info[:control_covariances][1] ≈ P[1] * Σ0 * P[1]'
    @test sol.solver_info[:state_control_covariances][1] ≈ -Σ0 * P[1]'
end

@testset "FNELQ/SFNELQ: failure statuses instead of false equilibria" begin
    one(x) = reshape([float(x)], 1, 1)
    mk(R1, R2) = JointLQGameProblem(one(1), [one(1), one(1)], [one(1), one(1)], [R1, R2],
                                    [one(0), one(0)], [1.0], 3.0; dt = 1.0)

    @testset "singular stacked Nash system" begin
        g = mk([1.0 1.0; 1.0 1.0], [1.0 1.0; 1.0 1.0])
        sol = @test_logs (:warn,) match_mode = :any solve(g, FNELQ())
        @test !sol.converged
        @test sol.solver_info[:status] === :singular_nash_system
        sg = StochasticGameProblem(g, ProcessNoiseModel(one(0.1)))
        ssol = @test_logs (:warn,) match_mode = :any solve(sg, SFNELQ())
        @test !ssol.converged
        @test ssol.solver_info[:status] === :singular_nash_system
    end

    @testset "non-convex own-control block is not reported as Nash" begin
        g = mk([-1.0 0.0; 0.0 1.0], [1.0 0.0; 0.0 1.0])
        sol = @test_logs (:warn,) match_mode = :any solve(g, FNELQ())
        @test !sol.converged
        @test sol.solver_info[:status] === :second_order_violation
        @test sol.solver_info[:min_own_hessian_eig] < 0
    end

    @testset "well-posed game reports :solved" begin
        g = mk([1.0 0.2; 0.2 1.0], [1.0 -0.1; -0.1 1.0])
        sol = solve(g, FNELQ())
        @test sol.converged && sol.solver_info[:status] === :solved
    end
end

@testset "SFNELQ: capability guards" begin
    A = [1.0 0.1; 0.0 1.0]; B = [reshape([0.0, 0.1], 2, 1)]
    g = LQGameProblem(A, B, [Matrix(1.0I, 2, 2)], [reshape([1.0], 1, 1)],
                      [Matrix(1.0I, 2, 2)], [0.0, 0.0], 3.0; dt = 1.0)
    W = ProcessNoiseModel(Matrix(0.01I, 2, 2))
    nonlinear = SharedObservation((x, v, p, t) -> x[1:1] .^ 2 .+ v, 1, reshape([0.1], 1, 1))
    @test_throws ErrorException solve(StochasticGameProblem(g, W; information = nonlinear), SFNELQ())
    @test_throws ErrorException solve(StochasticGameProblem(g, W; information = OpenLoopInformation()), SFNELQ())
    cc = ChanceConstraint(StateBounds(1; state_offset = 0, state_dim = 2,
                                      lower = [-5.0, -5.0], upper = [1.0, 5.0]), 0.05)
    @test_throws ErrorException solve(StochasticGameProblem(g, W; chance_constraints = [cc]), SFNELQ())
    @test solve(StochasticGameProblem(g, W), SFNELQ()).converged
end
