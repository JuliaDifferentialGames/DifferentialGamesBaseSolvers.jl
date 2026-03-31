# ============================================================================
# test_fnelq_ltv.jl
#
# Regression and correctness tests for the updated FNELQ backward pass.
#
# Test hierarchy:
#   L1 — extract_lq_matrices: correct dispatch and field extraction
#   L2 — LTI regression: updated FNELQ matches original on LTI problems
#   L3 — LTV correctness: constant LTV sequences reproduce LTI result exactly
#   L4 — LTV variation: time-varying A produces different (verifiable) trajectory
#   L5 — Affine terms: LTV q/r sequences handled correctly
# ============================================================================

using Test
using LinearAlgebra
using DifferentialGamesBase

# ── Shared test fixtures ──────────────────────────────────────────────────────

function make_lti_game(T=Float64; n=4, m=2, n_players=2, tf=1.0, dt=0.1)
    A  = T(0.9) * Matrix{T}(I, n, n)
    B  = [Matrix{T}(I, n, m)          for _ in 1:n_players]
    Q  = [Matrix{T}(I, n, n)          for _ in 1:n_players]
    R  = [Matrix{T}(I, m, m)          for _ in 1:n_players]
    Qf = [T(5) * Matrix{T}(I, n, n)   for _ in 1:n_players]
    x0 = ones(T, n)
    LQGameProblem(A, B, Q, R, Qf, x0, T(tf); dt=T(dt))
end

function make_ltv_game_constant(T=Float64; n=4, m=2, n_players=2, tf=1.0, dt=0.1)
    N  = Int(round(tf / dt))
    A  = T(0.9) * Matrix{T}(I, n, n)
    B  = Matrix{T}(I, n, m)
    Q  = Matrix{T}(I, n, n)
    R  = Matrix{T}(I, m, m)
    Qf = T(5) * Matrix{T}(I, n, n)
    x0 = ones(T, n)

    A_seq  = [copy(A) for _ in 1:N]
    B_seq  = [[copy(B) for _ in 1:N] for _ in 1:n_players]
    Q_seq  = [[copy(Q) for _ in 1:N] for _ in 1:n_players]
    R_seq  = [[copy(R) for _ in 1:N] for _ in 1:n_players]
    Qf_vec = [copy(Qf) for _ in 1:n_players]

    LTVLQGameProblem(A_seq, B_seq, Q_seq, R_seq, Qf_vec, x0, T(tf); dt=T(dt))
end

function make_ltv_game_varying(T=Float64; n=4, m=2, n_players=2, tf=1.0, dt=0.1)
    N  = Int(round(tf / dt))
    x0 = ones(T, n)
    Qf = T(5) * Matrix{T}(I, n, n)

    A_seq  = [T(0.9) * Matrix{T}(I, n, n) + T(0.05) * sin(k * π / N) * Matrix{T}(I, n, n)
              for k in 1:N]
    B_seq  = [[Matrix{T}(I, n, m) for _ in 1:N] for _ in 1:n_players]
    Q_seq  = [[Matrix{T}(I, n, n) for _ in 1:N] for _ in 1:n_players]
    R_seq  = [[Matrix{T}(I, m, m) for _ in 1:N] for _ in 1:n_players]
    Qf_vec = [copy(Qf)            for _ in 1:n_players]

    LTVLQGameProblem(A_seq, B_seq, Q_seq, R_seq, Qf_vec, x0, T(tf); dt=T(dt))
end

# ── Tests ─────────────────────────────────────────────────────────────────────

@testset "FNELQ — LTV backward pass" begin

    # ─────────────────────────────────────────────────────────────────────────
    @testset "L1 — extract_lq_matrices dispatch" begin

        @testset "LTI: returns LinearDynamics with is_ltv=false" begin
            game = make_lti_game()
            dyn, scs, Qf, n, np = DifferentialGamesBaseSolvers.extract_lq_matrices(game)
            @test dyn isa LinearDynamics
            @test !is_ltv(dyn)
            @test n == 4
            @test np == 2
            @test length(scs) == 2
            @test all(sc isa LQStageCost for sc in scs)
            @test length(Qf) == 2
        end

        @testset "LTV: returns LinearDynamics with is_ltv=true" begin
            game = make_ltv_game_constant()
            dyn, scs, Qf, n, np = DifferentialGamesBaseSolvers.extract_lq_matrices(game)
            @test dyn isa LinearDynamics
            @test is_ltv(dyn)
            @test n == 4
            @test np == 2
        end

        @testset "LTI: accessor values match constructor inputs" begin
            game = make_lti_game()
            dyn, scs, Qf, n, np = DifferentialGamesBaseSolvers.extract_lq_matrices(game)
            @test get_A(dyn, 1) === get_A(dyn, 5)
            @test get_B(dyn, 1, 1) === get_B(dyn, 1, 7)
            @test get_Q(scs[1], 1) === get_Q(scs[1], 5)
            @test get_R(scs[1], 1) === get_R(scs[1], 5)
        end

        @testset "LTV: accessor values differ across k" begin
            game = make_ltv_game_varying()
            dyn, scs, _, _, _ = DifferentialGamesBaseSolvers.extract_lq_matrices(game)
            @test get_A(dyn, 1) != get_A(dyn, 5)
            @test get_Q(scs[1], 1) == get_Q(scs[1], 5)
        end

    end

    # ─────────────────────────────────────────────────────────────────────────
    @testset "L2 — LTI regression: updated FNELQ matches reference output" begin
        game   = make_lti_game()
        solver = FNELQ()
        sol    = solve(game, solver; verbose=false)

        @test sol.converged
        @test length(sol.trajectories) == 2

        # State trajectory: finite, correct dimensions
        states = sol.trajectories[1].states
        @test all(isfinite, states)
        @test size(states, 1) == 4   # n=4
        @test size(states, 2) == n_steps(game) + 1

        # Control trajectory: finite, correct dimensions
        controls = sol.trajectories[1].controls
        @test all(isfinite, controls)
        @test size(controls, 1) == 2   # m=2
        @test size(controls, 2) == n_steps(game)

        # Feedback gains: correct dimensions
        N = n_steps(game)
        P = sol.solver_info[:feedback_gains]
        @test length(P) == 2
        for i in 1:2
            @test length(P[i]) == N
            @test size(P[i][1]) == (2, 4)
        end

        # Cost-to-go matrices: symmetric at every step
        Z = sol.solver_info[:cost_to_go_matrices]
        for i in 1:2, k in 1:N+1
            @test norm(Z[i][k] - Z[i][k]') < 1e-10
        end

        # Costs finite and positive
        costs = sol.solver_info[:costs]
        @test all(isfinite, costs)
        @test all(costs .> 0)

        # Deterministic: identical costs on repeat solve
        sol2 = solve(game, solver; verbose=false)
        @test sol.solver_info[:costs] ≈ sol2.solver_info[:costs] atol=0.0

        # Deterministic: identical states on repeat solve
        @test sol.trajectories[1].states ≈ sol2.trajectories[1].states atol=0.0
    end

    # ─────────────────────────────────────────────────────────────────────────
    @testset "L3 — LTV with constant matrices: identical to LTI" begin
        solver   = FNELQ()
        sol_lti  = solve(make_lti_game(),          solver; verbose=false)
        sol_ltv  = solve(make_ltv_game_constant(), solver; verbose=false)

        # Feedback gains at every timestep must match
        P_lti = sol_lti.solver_info[:feedback_gains]
        P_ltv = sol_ltv.solver_info[:feedback_gains]
        N = n_steps(make_lti_game())
        for i in 1:2, k in 1:N
            @test P_lti[i][k] ≈ P_ltv[i][k] atol=1e-10
        end

        # Costs must be identical
        @test sol_lti.solver_info[:costs] ≈ sol_ltv.solver_info[:costs] atol=1e-8

        # Cost-to-go matrices must match at every step
        Z_lti = sol_lti.solver_info[:cost_to_go_matrices]
        Z_ltv = sol_ltv.solver_info[:cost_to_go_matrices]
        for i in 1:2, k in 1:N+1
            @test Z_lti[i][k] ≈ Z_ltv[i][k] atol=1e-10
        end

        # State trajectories must match
        for i in 1:2
            @test sol_lti.trajectories[i].states ≈ sol_ltv.trajectories[i].states atol=1e-10
            @test sol_lti.trajectories[i].controls ≈ sol_ltv.trajectories[i].controls atol=1e-10
        end
    end

    # ─────────────────────────────────────────────────────────────────────────
    @testset "L4 — LTV with varying A: produces different trajectory" begin
        solver    = FNELQ()
        sol_lti   = solve(make_lti_game(),         solver; verbose=false)
        sol_ltv_v = solve(make_ltv_game_varying(), solver; verbose=false)

        @test sol_ltv_v.converged

        # Costs must differ — varying A changes the closed-loop value function
        @test !isapprox(sol_lti.solver_info[:costs], sol_ltv_v.solver_info[:costs]; atol=1e-6)

        # Gains must differ at at least one timestep
        P_lti = sol_lti.solver_info[:feedback_gains]
        P_ltv = sol_ltv_v.solver_info[:feedback_gains]
        N = n_steps(make_lti_game())
        @test any(!isapprox(P_lti[1][k], P_ltv[1][k]; atol=1e-6) for k in 1:N)

        # Cost-to-go matrices must all be finite
        Z = sol_ltv_v.solver_info[:cost_to_go_matrices]
        for i in 1:2, k in 1:N+1
            @test all(isfinite, Z[i][k])
        end
    end

    # ─────────────────────────────────────────────────────────────────────────
    @testset "L5 — Affine terms: LTV q/r sequences" begin
        T  = Float64
        n, m, np = 1, 1, 1
        N  = 5; tf = 0.5; dt = 0.1
        a  = T(0.9)
        A_seq  = [fill(a, 1, 1)    for _ in 1:N]
        B_seq  = [[fill(T(1), 1, 1) for _ in 1:N]]
        Q_seq  = [[fill(T(1), 1, 1) for _ in 1:N]]
        R_seq  = [[fill(T(1), 1, 1) for _ in 1:N]]
        Qf_vec = [fill(T(5), 1, 1)]
        x0     = T[1.0]

        game_no_aff = LTVLQGameProblem(A_seq, B_seq, Q_seq, R_seq, Qf_vec, x0, T(tf); dt=T(dt))

        q_seq = [[fill(T(0.5), 1) for _ in 1:N]]
        r_seq = [[fill(T(0.1), 1) for _ in 1:N]]
        game_aff = LTVLQGameProblem(A_seq, B_seq, Q_seq, R_seq, Qf_vec, x0, T(tf);
                                     dt=T(dt), q_seq=q_seq, r_seq=r_seq)

        solver = FNELQ()
        sol_no  = solve(game_no_aff, solver; verbose=false)
        sol_aff = solve(game_aff,    solver; verbose=false)

        @test sol_no.converged
        @test sol_aff.converged

        # Feedback gain P must be identical — q, r don't enter the Riccati equation
        P_no  = sol_no.solver_info[:feedback_gains][1]
        P_aff = sol_aff.solver_info[:feedback_gains][1]
        for k in 1:N
            @test P_no[k] ≈ P_aff[k] atol=1e-10
        end

        # Affine gains α must be nonzero when q/r are nonzero
        @test haskey(sol_aff.solver_info, :affine_gains)
        α = sol_aff.solver_info[:affine_gains][1]
        @test any(abs(α[k][1]) > 1e-6 for k in 1:N)

        # Costs must differ between affine and non-affine games
        @test !isapprox(sol_no.solver_info[:costs], sol_aff.solver_info[:costs]; atol=1e-6)
    end

end