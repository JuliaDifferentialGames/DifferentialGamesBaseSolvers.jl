using Test
using LinearAlgebra
using DifferentialGamesBase
using DifferentialGamesBaseSolvers

# ============================================================================
# Shared test fixtures
# ============================================================================

"""
Scalar 1-player 1-state 1-control LQ game with known closed-form solution.

  min  Σ_{k=0}^{N-1} (½q·x² + ½r·u²) + ½qf·xf²
  s.t. x(k+1) = a·x(k) + b·u(k)

Optimal gain: P(k) = qf, then backwards:
  Z(k) = ... standard Riccati.
For scalar, closed-form at N=1:
  S  = r + b²·qf
  P  = b·qf·a / S
  u* = -P·x₀
"""
function scalar_lq_game(;
    a=1.0, b=1.0, q=1.0, r=1.0, qf=10.0, x0=1.0, N=10, dt=0.1
)
    A  = reshape([a], 1, 1)
    B  = [reshape([b], 1, 1)]
    Q  = [reshape([q], 1, 1)]
    R  = [reshape([r], 1, 1)]
    Qf_arr = [reshape([qf], 1, 1)]
    game = LQGameProblem(A, B, Q, R, Qf_arr, [x0], Float64(N*dt); dt=dt)
    return game
end

"""
Two-player double integrator game (shared state, n=4, m=[2,2]).

Player 1: regulate state to zero + control cost.
Player 2: drive x₃ to 1 + control cost.
"""
function two_player_double_integrator(; N=20, dt=0.1)
    n = 4
    A = [1.0 dt 0.0 0.0;
         0.0 1.0 0.0 0.0;
         0.0 0.0 1.0 dt;
         0.0 0.0 0.0 1.0]
    B1 = [0.0; dt; 0.0; 0.0;; ]         # n×1
    B2 = [0.0; 0.0; 0.0; dt;; ]         # n×1

    Q1  = diagm([1.0, 0.1, 0.0, 0.0])
    Q2  = diagm([0.0, 0.0, 1.0, 0.1])
    R1  = 0.1 * Matrix{Float64}(I, 1, 1)
    R2  = 0.1 * Matrix{Float64}(I, 1, 1)
    Qf1 = diagm([10.0, 1.0, 0.0, 0.0])
    Qf2 = diagm([0.0, 0.0, 10.0, 1.0])
    x0  = [2.0, 0.0, -1.0, 0.0]

    game = LQGameProblem(A, [B1, B2], [Q1, Q2], [R1, R2], [Qf1, Qf2], x0, Float64(N*dt); dt=dt)
    return game
end

"""
Two-player nonlinear game: double integrator dynamics with nonlinear costs.
Used to test iLQGames against FNELQ (iLQGames should recover FNELQ solution
on the LQ version of this problem up to linearisation tolerance).
"""
function nonlinear_two_player(; N=20, dt=0.1)
    # Dynamics: two independent 2D double integrators (4 states total)
    f1 = (x, u, p, t) -> [x[2]; u[1]]
    f2 = (x, u, p, t) -> [x[2]; u[1]]
    dyn = SeparableDynamics(
        [ContinuousPlayerDynamics(f1, 2, 1),
         ContinuousPlayerDynamics(f2, 2, 1)],
        [2, 2], [1, 1]
    )

    # Quadratic costs (so iLQGames = FNELQ on converged LQ approximation)
    Q1  = diagm([1.0, 0.1])
    Q2  = diagm([1.0, 0.1])
    Qf1 = diagm([10.0, 1.0])
    Qf2 = diagm([10.0, 1.0])
    R1  = reshape([0.1], 1, 1)
    R2  = reshape([0.1], 1, 1)

    o1 = PlayerObjective(1, LQStageCost(Q1, R1), LQTerminalCost(Qf1))
    o2 = PlayerObjective(2, LQStageCost(Q2, R2), LQTerminalCost(Qf2))

    p1 = Player{Float64}(1, 2, 1, [1.0, 0.0], f1, o1, AbstractPrivateConstraint[])
    p2 = Player{Float64}(2, 2, 1, [-1.0, 0.0], f2, o2, AbstractPrivateConstraint[])

    game = DifferentialGame([p1, p2], Float64(N*dt), Float64(dt))
    return game
end

# ============================================================================
# FNELQ tests
# ============================================================================

@testset "FNELQ: Scalar LQ — analytical ground truth" begin
    # Single player, scalar state/control, N=2 steps (dt=0.5, tf=1.0).
    # DiscreteTime requires dt < tf, so N=1/dt=1 is forbidden.
    #
    # The backward pass runs k=2 then k=1. At k=N=2, Z(N)=Qf exactly:
    #   S₂  = r + b²·Qf          ← coupled system at last step
    #   P₂  = (b·Qf·a) / S₂     ← gain at k=2
    #   u*(2) = -P₂ · x(2)
    a=1.0; b=1.0; q=1.0; r=1.0; qf=10.0; x0=2.0; N=2

    game   = scalar_lq_game(a=a, b=b, q=q, r=r, qf=qf, x0=x0, N=N, dt=0.5)
    solver = FNELQ()
    sol    = solve(game, solver; check_compatibility=false)

    # Closed-form for last backward step (k=N)
    S_N  = r + b^2 * qf
    P_N  = (b * qf * a) / S_N

    @testset "Solution status" begin
        @test sol.converged
        @test sol.equilibrium_type == :FeedbackNash
        @test has_strategy(sol)
        @test sol.strategy isa FeedbackStrategy
    end

    @testset "Gain at k=N matches closed form" begin
        P_num = get_gain(sol.strategy, 1, N)[1, 1]
        @test P_num ≈ P_N atol=1e-10
    end

    @testset "Control at k=N satisfies optimality" begin
        x_N   = sol.state_trajectory[1, N]
        u_N   = sol.trajectories[1].controls[1, N]
        @test u_N ≈ -P_N * x_N atol=1e-10
    end

    @testset "State trajectory consistent with dynamics" begin
        X = sol.state_trajectory
        U = sol.trajectories[1].controls
        for k in 1:N
            @test X[1, k+1] ≈ a * X[1, k] + b * U[1, k] atol=1e-10
        end
    end

    @testset "Cost is non-negative" begin
        @test sol.trajectories[1].cost >= 0.0
    end
end

@testset "FNELQ: Two-player double integrator" begin
    game   = two_player_double_integrator()
    solver = FNELQ()
    sol    = solve(game, solver; check_compatibility=false)

    @testset "Solution structure" begin
        @test sol.converged
        @test sol.equilibrium_type == :FeedbackNash
        @test length(sol.trajectories) == 2
        @test has_strategy(sol)
        @test sol.strategy isa FeedbackStrategy
        @test sol.strategy.n_players == 2
    end

    @testset "Trajectory dimensions" begin
        N  = n_steps(game)
        @test size(sol.state_trajectory) == (4, N+1)
        @test size(sol.trajectories[1].controls) == (1, N)
        @test size(sol.trajectories[2].controls) == (1, N)
    end

    @testset "Initial state matches" begin
        @test sol.state_trajectory[:, 1] ≈ game.initial_state
    end

    @testset "Trajectory is dynamically feasible" begin
        # x(k+1) = A·x(k) + B₁·u₁(k) + B₂·u₂(k)
        N   = n_steps(game)
        dyn = game.dynamics
        A   = get_A(dyn, 1)
        B1  = get_B(dyn, 1, 1)
        B2  = get_B(dyn, 2, 1)
        X   = sol.state_trajectory
        U1  = sol.trajectories[1].controls
        U2  = sol.trajectories[2].controls
        for k in 1:N
            x_next = A * X[:, k] + B1 * U1[:, k] + B2 * U2[:, k]
            @test X[:, k+1] ≈ x_next atol=1e-10
        end
    end

    @testset "Nash condition: gains are finite and non-trivial" begin
        strat = sol.strategy
        N     = n_steps(game)
        # At least some gains should be non-zero (non-trivial Nash)
        P1_norms = [norm(get_gain(strat, 1, k)) for k in 1:N]
        P2_norms = [norm(get_gain(strat, 2, k)) for k in 1:N]
        @test all(isfinite, P1_norms)
        @test all(isfinite, P2_norms)
        @test maximum(P1_norms) > 1e-6
        @test maximum(P2_norms) > 1e-6
    end

    @testset "Cost values are finite and positive" begin
        @test all(c -> isfinite(c) && c > 0, get_costs(sol))
    end
end

@testset "FNELQ: LTV game — costs decrease vs zero-gain" begin
    # For an LQ game, optimal control should give lower cost than zero control.
    game   = two_player_double_integrator(N=15, dt=0.1)
    solver = FNELQ()
    sol    = solve(game, solver; check_compatibility=false)

    N   = n_steps(game)
    dyn = game.dynamics
    A   = get_A(dyn, 1)
    B1  = get_B(dyn, 1, 1)
    B2  = get_B(dyn, 2, 1)

    # Compute cost under zero control
    X_zero = zeros(4, N+1)
    X_zero[:, 1] = game.initial_state
    for k in 1:N
        X_zero[:, k+1] = A * X_zero[:, k]
    end
    costs_zero = [begin
        obj = get_objective(game, i)
        total_cost(obj, [X_zero[:, k] for k in 1:N+1],
                   [zeros(1) for _ in 1:N], nothing)
    end for i in 1:2]

    for i in 1:2
        @test sol.trajectories[i].cost < costs_zero[i]
    end
end

@testset "FNELQ: Affine costs (cross-term M)" begin
    # Build a game with non-zero cross-term M and verify M is used.
    n = 2; m = 1; dt = 0.1; N = 5
    A  = [1.0 dt; 0.0 1.0]
    B  = [reshape([0.0; dt], 2, 1)]
    Q  = [Matrix{Float64}(I, 2, 2)]
    R  = [reshape([1.0], 1, 1)]
    M  = [reshape([0.1; 0.2], 2, 1)]   # cross-term coupling state and control
    Qf = [diagm([5.0, 1.0])]
    x0 = [1.0, 0.0]

    game_M  = LQGameProblem(A, B, Q, R, Qf, x0, Float64(N*dt); dt=dt, M=[M[1]])
    game_0  = LQGameProblem(A, B, Q, R, Qf, x0, Float64(N*dt); dt=dt)

    sol_M = solve(game_M, FNELQ(); check_compatibility=false)
    sol_0 = solve(game_0, FNELQ(); check_compatibility=false)

    @testset "M≠0 gives different gains than M=0" begin
        for k in 1:N
            P_M = get_gain(sol_M.strategy, 1, k)
            P_0 = get_gain(sol_0.strategy, 1, k)
            @test !(P_M ≈ P_0)
        end
    end

    @testset "Both solutions converge" begin
        @test sol_M.converged
        @test sol_0.converged
    end
end

@testset "FNELQ: Warmstart does not affect solution" begin
    # FNELQ is an exact backward pass — warmstart is ignored.
    # Two independent solves should give identical results.
    game = two_player_double_integrator()
    sol1 = solve(game, FNELQ(); check_compatibility=false)
    sol2 = solve(game, FNELQ(); warmstart=sol1, check_compatibility=false)
    @test get_costs(sol1) ≈ get_costs(sol2)
end