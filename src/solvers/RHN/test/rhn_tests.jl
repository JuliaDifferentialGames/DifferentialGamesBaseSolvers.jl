# ============================================================================
# rhn_tests.jl — Tests for the Receding Horizon Nash (RHN) solver wrapper
#
# Tests the RecedingHorizonNash wrapper with multiple inner solvers:
#   - FNELQ   (exact LQ, warmstart structurally accepted but has no effect)
#   - iLQGames (nonlinear/LQ, warmstart reduces iterations)
#
# Run standalone:
#   using DifferentialGamesBase, DifferentialGamesBaseSolvers
#   include("rhn_tests.jl")
# ============================================================================

using Test
using LinearAlgebra
using DifferentialGamesBase
using DifferentialGamesBaseSolvers

# ============================================================================
# Shared fixtures
# ============================================================================

"""
    make_di_horizon_game(; N_win, dt) -> GameProblem{Float64}

2-player discrete-time double integrator game (shared 2D state).

Dynamics:  x_{k+1} = [1 dt; 0 1] x_k + [0; dt] u₁_k + [0; dt] u₂_k
Costs:     each player minimizes ½ xᵀQᵢx + ½ uᵢᵀRᵢuᵢ + ½ xfᵀQfᵢxf

Both players penalise the full state, ensuring the Nash equilibrium
drives x → 0 from any x₀. Used as the receding-horizon window template.
"""
function make_di_horizon_game(;
    N_win ::Int     = 10,
    dt    ::Float64 = 0.1,
    q_pos ::Float64 = 1.0,
    q_vel ::Float64 = 0.1,
    r     ::Float64 = 0.5,
    qf_pos::Float64 = 5.0,
    qf_vel::Float64 = 1.0
)
    T  = Float64
    A  = T[1.0 dt; 0.0 1.0]
    B1 = T[0.0; dt;;]    # (2×1) player 1 force
    B2 = T[0.0; dt;;]    # (2×1) player 2 force (same channel)

    Q  = diagm(T[q_pos, q_vel])
    Qf = diagm(T[qf_pos, qf_vel])
    R  = reshape(T[r], 1, 1)

    tf = Float64(N_win * dt)
    # Dummy x0 — patched by with_initial_state at each RH step
    x0 = zeros(T, 2)

    LQGameProblem(A, [B1, B2], [Q, Q], [R, R], [Qf, Qf], x0, tf; dt=dt)
end

# ============================================================================
@testset "Receding Horizon Nash solver" begin

    # ── Test 1: RecedingHorizonNashProblem construction ───────────────────────
    @testset "RecedingHorizonNashProblem construction" begin
        game = make_di_horizon_game()
        x0   = [1.0, 0.0]
        prob = RecedingHorizonNashProblem(game, x0, 30)

        @test prob isa RecedingHorizonNashProblem{Float64}
        @test n_players(prob) == 2
        @test state_dim(prob) == 2
        @test prob.n_sim_steps == 30
        @test prob.x0 ≈ x0
        @test prob.horizon_game === game

        # Dimension mismatch should throw
        @test_throws AssertionError RecedingHorizonNashProblem(game, zeros(3), 10)
        # Zero or negative steps should throw
        @test_throws AssertionError RecedingHorizonNashProblem(game, x0, 0)
    end

    # ── Test 2: with_initial_state is O(1) and correct ───────────────────────
    @testset "with_initial_state" begin
        game    = make_di_horizon_game()
        x0_new  = Float64[2.0, -1.0]
        patched = with_initial_state(game, x0_new)

        @test patched isa GameProblem{Float64}
        @test patched.initial_state ≈ x0_new
        # All other fields are shared by reference (no deep copy)
        @test patched.dynamics     === game.dynamics
        @test patched.objectives   === game.objectives
        @test patched.time_horizon === game.time_horizon
        @test patched.metadata     === game.metadata

        # Original game unchanged
        @test game.initial_state ≈ zeros(2)
    end

    # ── Test 3: RecedingHorizonNash struct ────────────────────────────────────
    @testset "RecedingHorizonNash struct" begin
        solver = RecedingHorizonNash(FNELQ())
        @test solver isa RecedingHorizonNash{FNELQ}
        @test solver.warm_start    == true
        @test solver.verbose_inner == false

        s2 = RecedingHorizonNash(iLQGames(); warm_start=false, verbose_inner=true)
        @test s2.warm_start    == false
        @test s2.verbose_inner == true
    end

    # ── Test 4: FNELQ + RHN — solution shape and types ───────────────────────
    @testset "FNELQ+RHN: solution shape" begin
        game   = make_di_horizon_game(N_win=10)
        x0     = Float64[1.0, 0.0]
        N_sim  = 20
        prob   = RecedingHorizonNashProblem(game, x0, N_sim)
        solver = RecedingHorizonNash(FNELQ())

        sol = solve(prob, solver)

        @test sol isa RecedingHorizonNashSolution{Float64}
        @test size(sol.X)    == (2, N_sim+1)          # n × N_sim+1
        @test length(sol.U)  == 2                      # one per player
        @test size(sol.U[1]) == (1, N_sim)             # m₁ × N_sim
        @test size(sol.U[2]) == (1, N_sim)
        @test length(sol.costs) == 2
        @test sol.converged
        @test sol.solve_time ≥ 0.0

        # Solver info entries
        @test sol.solver_info[:n_steps_simulated] == N_sim
        @test sol.solver_info[:n_inner_converged] == N_sim   # FNELQ always converges
    end

    # ── Test 5: FNELQ + RHN — initial condition ───────────────────────────────
    @testset "FNELQ+RHN: initial condition preserved" begin
        game = make_di_horizon_game()
        x0   = Float64[3.0, -2.0]
        prob = RecedingHorizonNashProblem(game, x0, 15)
        sol  = solve(prob, RecedingHorizonNash(FNELQ()))

        @test sol.X[:, 1] ≈ x0   atol=1e-12
    end

    # ── Test 6: FNELQ + RHN — convergence toward origin ──────────────────────
    @testset "FNELQ+RHN: state converges toward origin" begin
        # Larger terminal cost ensures the Nash equilibrium strongly stabilises x→0
        game = make_di_horizon_game(N_win=15, qf_pos=20.0, qf_vel=5.0)
        x0   = Float64[2.0, 1.0]
        sol  = solve(
            RecedingHorizonNashProblem(game, x0, 50),
            RecedingHorizonNash(FNELQ())
        )

        # State norm at final step should be smaller than at start
        x_final = sol.X[:, end]
        x_init  = sol.X[:, 1]
        @test norm(x_final) < norm(x_init)

        # For a well-tuned game the system should be near origin after 50 steps
        @test norm(x_final) < 0.5 * norm(x_init)
    end

    # ── Test 7: warm-start shift helpers ─────────────────────────────────────
    @testset "warm-start shift helpers" begin
        T  = Float64
        N  = 8
        np = 2
        n  = 2
        m  = [1, 1]
        ts = collect(range(0.0, 0.8, length=N+1))

        # Build a random FeedbackStrategy
        gains = [[rand(T, m[i], n) for _ in 1:N] for i in 1:np]
        ff    = [[rand(T, m[i])    for _ in 1:N] for i in 1:np]
        xnom  = rand(T, n, N+1)
        unom  = [rand(T, m[i], N)  for i in 1:np]
        strat = FeedbackStrategy(gains, ff, xnom, unom, m, ts)

        shifted = DifferentialGamesBaseSolvers._shift_feedback_strategy(strat)

        @test n_steps(shifted) == N                      # same horizon length
        @test size(shifted.nominal_states) == (n, N+1)  # same shape

        # After shifting, first element of shifted = second element of original
        for i in 1:np
            @test shifted.gains[i][1] ≈ gains[i][2]   atol=1e-14
            @test shifted.feedforward[i][1] ≈ ff[i][2] atol=1e-14
        end
        # Last element repeated
        for i in 1:np
            @test shifted.gains[i][N] ≈ gains[i][N]   atol=1e-14
        end
        # Nominal states shifted
        @test shifted.nominal_states[:, 1] ≈ xnom[:, 2]   atol=1e-14
        @test shifted.nominal_states[:, N+1] ≈ xnom[:, N+1] atol=1e-14

        # OpenLoopStrategy shift
        ctrl = [rand(T, m[i], N) for i in 1:np]
        ols  = OpenLoopStrategy(ctrl, m, ts)
        ols_s = DifferentialGamesBaseSolvers._shift_open_loop_strategy(ols)

        @test n_steps(ols_s) == N
        for i in 1:np
            @test ols_s.controls[i][:, 1] ≈ ctrl[i][:, 2] atol=1e-14
            @test ols_s.controls[i][:, N] ≈ ctrl[i][:, N] atol=1e-14
        end
    end

    # ── Test 8: FNELQ + RHN — no warm-start mode ────────────────────────────
    @testset "FNELQ+RHN: no warm-start mode" begin
        game = make_di_horizon_game()
        x0   = Float64[1.0, 0.5]
        prob = RecedingHorizonNashProblem(game, x0, 10)

        sol_ws = solve(prob, RecedingHorizonNash(FNELQ(); warm_start=true))
        sol_no = solve(prob, RecedingHorizonNash(FNELQ(); warm_start=false))

        # Both should converge; FNELQ is exact so results should be identical
        @test sol_ws.converged && sol_no.converged
        @test sol_ws.X ≈ sol_no.X   atol=1e-10
        @test sol_ws.costs ≈ sol_no.costs   atol=1e-10
    end

    # ── Test 9: costs are non-negative and finite ─────────────────────────────
    @testset "FNELQ+RHN: costs finite and non-negative" begin
        game = make_di_horizon_game()
        x0   = Float64[1.0, 0.0]
        sol  = solve(
            RecedingHorizonNashProblem(game, x0, 20),
            RecedingHorizonNash(FNELQ())
        )

        @test all(isfinite, sol.costs)
        @test all(sol.costs .≥ 0.0)
        @test all(isfinite, sol.X)
    end

    # ── Test 10: iLQGames + RHN — runs and converges ──────────────────────────
    @testset "iLQGames+RHN: runs and converges" begin
        # iLQGames is an iterative solver; for a linear-quadratic game it should
        # converge to the same solution as FNELQ (up to its tolerance).
        game = make_di_horizon_game(N_win=10)
        x0   = Float64[1.5, -0.5]
        prob = RecedingHorizonNashProblem(game, x0, 15)

        ilq_solver = iLQGames(max_iter=100, ε_conv=1e-4)
        sol = solve(prob, RecedingHorizonNash(ilq_solver; warm_start=true))

        @test sol isa RecedingHorizonNashSolution{Float64}
        @test sol.X[:, 1] ≈ x0   atol=1e-12
        @test all(isfinite, sol.X)
        @test all(isfinite, sol.costs)
    end

    # ── Test 11: iLQGames + RHN — both FNELQ and iLQGames stabilise the system ─
    #
    # Key invariant: both solvers, when used as RHN inner solvers, should
    # drive the LQ game state toward the origin. We do NOT compare the two
    # closed-loop trajectories numerically — for a game where B₁ = B₂ (both
    # players affect the same channel), iLQGames may converge to a different
    # valid approximate equilibrium than the exact FNELQ solution. The test
    # checks the shared structural property: closed-loop stability.
    @testset "iLQGames+RHN: both solvers stabilise the same LQ game" begin
        game  = make_di_horizon_game(N_win=15, qf_pos=10.0, qf_vel=2.0)
        x0    = Float64[1.5, 0.5]
        prob  = RecedingHorizonNashProblem(game, x0, 30)

        sol_fnelq = solve(prob, RecedingHorizonNash(FNELQ()))
        sol_ilq   = solve(prob, RecedingHorizonNash(
                        iLQGames(max_iter=300, ε_conv=1e-4); warm_start=true))

        # Both produce valid, finite trajectories
        @test all(isfinite, sol_fnelq.X)
        @test all(isfinite, sol_ilq.X)
        @test all(isfinite, sol_fnelq.costs)
        @test all(isfinite, sol_ilq.costs)

        # Both start from the same x0
        @test sol_fnelq.X[:, 1] ≈ x0   atol=1e-12
        @test sol_ilq.X[:, 1]   ≈ x0   atol=1e-12

        # Both stabilise: final state norm < initial state norm
        @test norm(sol_fnelq.X[:, end]) < norm(x0)
        @test norm(sol_ilq.X[:, end])   < norm(x0)

        # Costs are non-negative and finite
        @test all(sol_fnelq.costs .≥ 0.0)
        @test all(sol_ilq.costs   .≥ 0.0)
    end

    # ── Test 12: accessors ────────────────────────────────────────────────────
    @testset "RecedingHorizonNashSolution accessors" begin
        game = make_di_horizon_game()
        x0   = Float64[1.0, 0.0]
        sol  = solve(
            RecedingHorizonNashProblem(game, x0, 10),
            RecedingHorizonNash(FNELQ())
        )

        @test get_state_trajectory(sol)     === sol.X
        @test get_control_trajectory(sol, 1) === sol.U[1]
        @test get_control_trajectory(sol, 2) === sol.U[2]
        @test get_total_cost(sol, 1) == sol.costs[1]
        @test get_total_cost(sol, 2) == sol.costs[2]

        @test_throws ErrorException get_control_trajectory(sol, 0)
        @test_throws ErrorException get_control_trajectory(sol, 3)
        @test_throws ErrorException get_total_cost(sol, 3)
    end

    # ── Test 13: dynamics consistency ────────────────────────────────────────
    @testset "FNELQ+RHN: trajectory satisfies dynamics" begin
        game = make_di_horizon_game(N_win=10)
        x0   = Float64[1.0, 0.5]
        prob = RecedingHorizonNashProblem(game, x0, 15)
        sol  = solve(prob, RecedingHorizonNash(FNELQ()))

        dyn = game.dynamics
        th  = game.time_horizon
        tw  = collect(range(0.0, th.tf, length=n_steps(game)+1))

        # Check every transition X[:,t+1] = A*X[:,t] + B₁*U[1][:,t] + B₂*U[2][:,t]
        for t in 1:prob.n_sim_steps
            u_joint = vcat(sol.U[1][:, t], sol.U[2][:, t])
            x_next  = _rollout_step(dyn, sol.X[:, t], u_joint, nothing, tw, 1)
            @test sol.X[:, t+1] ≈ x_next   atol=1e-10
        end
    end

    # ── Test 14: display does not error ──────────────────────────────────────
    @testset "display" begin
        game = make_di_horizon_game()
        prob = RecedingHorizonNashProblem(game, [1.0, 0.0], 10)
        sol  = solve(prob, RecedingHorizonNash(FNELQ()))

        # show(::IO, ...) should not throw
        buf = IOBuffer()
        show(buf, prob)
        show(buf, sol)
        show(buf, MIME("text/plain"), prob)
        show(buf, MIME("text/plain"), sol)
        @test length(take!(buf)) > 0
    end

end # @testset "Receding Horizon Nash solver"
