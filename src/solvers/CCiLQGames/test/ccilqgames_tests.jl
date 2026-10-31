using LinearAlgebra
using SparseArrays
using Random
using Test
using DifferentialGamesBase
using DifferentialGamesBaseSolvers

# Nonlinear stochastic game builder: joint state, per-player stage/terminal costs of (x, u_i).
function _cc_game(f, n, cd, stage, term, x0, N, dt; private = [], shared = [])
    np = length(cd)
    offs = [0; cumsum(cd)[1:end-1]]
    objs = [PlayerObjective(i, NonlinearStageCost(stage[i]), NonlinearTerminalCost(term[i])) for i in 1:np]
    meta = GameMetadata([n], cd, [0], offs, CouplingGraph(sparse(trues(np, np)), Vector{Int}[], nothing), false, nothing)
    GameProblem{Float64}(np, objs, CoupledNonlinearDynamics(f, n, sum(cd)), x0,
        AbstractPrivateConstraint[private...], AbstractSharedConstraint[shared...], DiscreteTime(N * dt, dt), meta)
end

@testset "CCiLQGames" begin
    @testset "unconstrained, noise-free LQ game reduces to FNELQ" begin
        A = [1.0 0.1; 0.0 1.0]; B = [reshape([0.0, 0.1], 2, 1)]
        g = LQGameProblem(A, B, [Matrix(1.0I, 2, 2)], [reshape([1.0], 1, 1)],
                          [Matrix(2.0I, 2, 2)], [1.0, 0.0], 1.0; dt = 0.1)
        sg = StochasticGameProblem(g, ProcessNoiseModel(zeros(2, 2)))
        det = solve(g, FNELQ())
        sol = solve(sg, CCiLQGames(mc_samples = 0))
        @test sol.converged && sol.solver_info[:status] === :locally_stationary_surrogate
        @test sol.trajectories[1].cost ≈ det.trajectories[1].cost rtol = 1e-6
        @test sol.state_trajectory ≈ det.state_trajectory atol = 1e-6
    end

    @testset "affine Gaussian chance constraint on an LQ game is met with probability 1-ε" begin
        dt, N = 0.1, 10
        A = Matrix(1.0I, 2, 2); B = [reshape([dt, 0.0], 2, 1), reshape([0.0, dt], 2, 1)]
        Q1 = [4.0 0.0; 0.0 0.0]; Q2 = [0.0 0.0; 0.0 4.0]
        t1 = [1.0, 0.0]; t2 = [0.0, 1.0]
        g = LQGameProblem(A, B, [Q1, Q2], [reshape([0.1], 1, 1), reshape([0.1], 1, 1)], [Q1, Q2],
                          [0.0, 0.0], N * dt; dt = dt, q = [-(Q1 * t1), -(Q2 * t2)])
        W = Matrix(0.01I, 2, 2)
        ε = 0.1
        cc = ChanceConstraint(StateBounds(1; state_offset = 0, state_dim = 1, lower = [-5.0], upper = [0.6]), ε)
        sg = StochasticGameProblem(g, ProcessNoiseModel(W); chance_constraints = [cc])
        sol = solve(sg, CCiLQGames(mc_samples = 20000, seed = 1))
        info = sol.solver_info
        @test info[:status] in (:locally_stationary_surrogate, :surrogate_feasible)
        @test info[:r_feas] <= 1e-4
        @test all(==(:exact_affine_gaussian), info[:guarantee_classes])
        # The upper-bound row (component 2) must bind somewhere, and then P(violation) = ε exactly.
        ub = filter(c -> c.component == 2, info[:chance_checks])
        @test all(c -> c.residual <= 1e-4, info[:chance_checks])
        # a two-row bound splits ε across its rows: the upper row carries ε/2
        @test maximum(c.violation_probability for c in ub) ≈ ε / 2 atol = 5e-3
        # Monte Carlo of the raw constraint agrees with the analytic per-event probabilities.
        mc = info[:monte_carlo]
        rows = [2 + 2 * (k - 1) for k in 1:N]          # component 2 at each step
        for (k, r) in enumerate(rows)
            p = ub[k].violation_probability
            lo, hi = mc.per_event_interval[r]
            @test lo - 2e-3 <= p <= hi + 2e-3
        end
        # Per-event ε does NOT give a joint guarantee: any-violation exceeds ε but obeys the union bound.
        @test mc.any_violation >= maximum(mc.per_event)
        @test mc.any_violation <= info[:analytic_joint_violation_bound] + 5e-3
        # Without the constraint the unconstrained equilibrium violates the bound far more often.
        free = solve(StochasticGameProblem(g, ProcessNoiseModel(W)), CCiLQGames(mc_samples = 0))
        @test maximum(free.state_trajectory[1, :]) > 0.6
        @test maximum(sol.state_trajectory[1, 1:N]) < 0.6
    end

    @testset "nonlinear 2-player game with a private chance constraint" begin
        c = 0.5; dt, N = 0.1, 20
        f = (x, u, p, t) -> [x[2]; u[1] - c * abs(x[2]) * x[2]; x[4]; u[2] - c * abs(x[4]) * x[4]]
        stage = [(x, u, p, t) -> 2.0 * (x[1] - 1.0)^2 + 0.5 * (x[1] - x[3])^2 + 0.05 * x[2]^2 + 0.1 * u[1]^2,
                 (x, u, p, t) -> 2.0 * (x[3] - 0.5)^2 + 0.5 * (x[1] - x[3])^2 + 0.05 * x[4]^2 + 0.1 * u[1]^2]
        term = [(x, p) -> 5.0 * (x[1] - 1.0)^2 + x[2]^2, (x, p) -> 5.0 * (x[3] - 0.5)^2 + x[4]^2]
        cc = ChanceConstraint(PrivateInequality(1; func = (x, u, p, t) -> [x[1] - 0.8], dim = 1), 0.05)
        g = _cc_game(f, 4, [1, 1], stage, term, zeros(4), N, dt)
        G = [0.0 0.0; 0.05 0.0; 0.0 0.0; 0.0 0.05]
        sde = ItoSDEDynamics(g.dynamics, constant_diffusion(G), 2)
        sg = StochasticGameProblem(g, sde; chance_constraints = [cc])
        sol = solve(sg, CCiLQGames(mc_samples = 4000, seed = 2, max_iter = 30))
        info = sol.solver_info
        @test info[:status] in (:locally_stationary_surrogate, :surrogate_feasible)
        @test info[:r_feas] <= 1e-4
        @test info[:guarantee_classes] == [:linearized_approximation]
        @test all(isfinite, sol.state_trajectory)
        # the deterministic mean path stays strictly inside the bound by the risk margin
        @test maximum(sol.state_trajectory[1, 1:N]) < 0.8
        @test info[:monte_carlo].per_event[argmax(info[:monte_carlo].per_event)] < 0.15
        # reproducible
        sol2 = solve(sg, CCiLQGames(mc_samples = 4000, seed = 2, max_iter = 30))
        @test sol2.state_trajectory ≈ sol.state_trajectory
    end

    @testset "3-player game with a shared chance constraint" begin
        dt, N = 0.1, 20
        f = (x, u, p, t) -> [x[2]; u[1]; x[4]; u[2]; x[6]; u[3]]
        pos = [1, 3, 5]
        stage = [(x, u, p, t) -> 2.0 * (x[1] - 1.0)^2 + 0.1 * x[2]^2 + 0.1 * u[1]^2,
                 (x, u, p, t) -> 2.0 * (x[3] - 0.0)^2 + 0.1 * x[4]^2 + 0.1 * u[1]^2,
                 (x, u, p, t) -> 1.0 * (x[5] - 0.5)^2 + 0.1 * x[6]^2 + 0.1 * u[1]^2]
        term = [(x, p) -> 5.0 * (x[1] - 1.0)^2 + x[2]^2, (x, p) -> 5.0 * x[3]^2 + x[4]^2,
                (x, p) -> 2.0 * (x[5] - 0.5)^2 + x[6]^2]
        # players 1 and 2 want to swap sides but must keep p₂ - p₁ ≥ 0.2 (shared)
        sc = ChanceConstraint(SharedInequality([1, 2, 3]; func = (x, u, p, t) -> [0.2 - (x[3] - x[1])], dim = 1), 0.05)
        x0 = [0.0, 0.0, 1.0, 0.0, 0.5, 0.0]
        g = _cc_game(f, 6, [1, 1, 1], stage, term, x0, N, dt)
        G = zeros(6, 3); G[2, 1] = 0.05; G[4, 2] = 0.05; G[6, 3] = 0.05
        sg = StochasticGameProblem(g, ItoSDEDynamics(g.dynamics, constant_diffusion(G), 3); chance_constraints = [sc])
        sol = solve(sg, CCiLQGames(mc_samples = 4000, seed = 3, max_iter = 30))
        info = sol.solver_info
        @test info[:status] in (:locally_stationary_surrogate, :surrogate_feasible)
        @test info[:r_feas] <= 1e-4
        @test minimum(sol.state_trajectory[3, 1:N] - sol.state_trajectory[1, 1:N]) >= 0.2
        # shared constraint: a single multiplier array is shared by all three players
        @test length(info[:multipliers]) == 1
    end

    @testset "informative failures and guards" begin
        A = [1.0 0.1; 0.0 1.0]; B = [reshape([0.0, 0.1], 2, 1)]
        g = LQGameProblem(A, B, [Matrix(1.0I, 2, 2)], [reshape([1.0], 1, 1)],
                          [Matrix(1.0I, 2, 2)], [1.0, 0.0], 1.0; dt = 0.1)
        sg = StochasticGameProblem(g, ProcessNoiseModel(zeros(2, 2)))
        @test_throws ErrorException solve(sg, CCiLQGames(covariance_mode = :consistent_covariance))
        shared = SharedObservation((x, v, p, t) -> x[1:1] .+ v, 1, reshape([0.1], 1, 1))
        @test_throws ErrorException solve(StochasticGameProblem(g, ProcessNoiseModel(zeros(2, 2)); information = shared), CCiLQGames())
        @test_throws ArgumentError CCiLQGames(covariance_mode = :bogus)
        @test solver_capabilities(CCiLQGames) ⊆ DifferentialGamesBase.CANONICAL_CAPABILITIES
    end
end
