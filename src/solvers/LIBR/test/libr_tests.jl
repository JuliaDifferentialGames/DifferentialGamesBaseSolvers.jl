# ============================================================================
# libr_tests.jl — Tests for the L-IBR solver
#
# Run as part of the package test suite from runtests.jl, or standalone via:
#   using DifferentialGamesBaseSolvers
#   include("libr_tests.jl")
# ============================================================================

using Test
using LinearAlgebra
using DifferentialGamesBase
using DifferentialGamesBaseSolvers

# ── Problem factory: 2-agent double-integrator lexicographic game ─────────────
#
# State: [px, py, vx, vy], Control: [ax, ay] (double integrator, Euler)
# Agent 1: start (1,0), goal (-1,0)
# Agent 2: start (-1,0), goal (1,0)
# z_i = (X_i :: 4×(N+1), U_i :: 2×N)

function _build_2agent_lex_game(; d_safe=0.3, tf=1.0, dt=0.1)
    T   = Float64
    dyn = (x, u, p, t) -> [x[3], x[4], u[1], u[2]]
    Q   = zeros(T, 4, 4)
    R   = Matrix{T}(I, 2, 2) .* 1e-6   # tiny placeholder — LIBR ignores these
    Qf  = zeros(T, 4, 4)

    function _player(id, x0)
        stage    = LQStageCost(Q, R)
        terminal = LQTerminalCost(Qf)
        obj      = PlayerObjective(id, stage, terminal)
        PlayerSpec{T}(id, 4, 2, x0, dyn, obj, T[])
    end

    p1 = _player(1, [1.0, 0.0, 0.0, 0.0])
    p2 = _player(2, [-1.0, 0.0, 0.0, 0.0])

    d2  = d_safe^2
    f_col = (zi, zj) -> begin
        Xi, _ = zi;  Xj, _ = zj
        val = zero(eltype(Xi))
        for k in 1:size(Xi, 2)
            Δx = Xi[1, k] - Xj[1, k];  Δy = Xi[2, k] - Xj[2, k]
            viol = max(zero(eltype(Xi)), d2 - Δx^2 - Δy^2)
            val += viol^2
        end
        val
    end

    col_pairs = Function[((zi, zj) -> zero(eltype(zi[1])))  f_col;
                          f_col  ((zi, zj) -> zero(eltype(zi[1])))]

    goal1 = [-1.0, 0.0];  goal2 = [1.0, 0.0]
    g1 = zi -> begin
        Xi, Ui = zi
        sum((Xi[1, k] - goal1[1])^2 + (Xi[2, k] - goal1[2])^2 for k in 1:size(Xi, 2)) +
        0.1 * sum(Ui[j, k]^2 for j in 1:size(Ui, 1), k in 1:size(Ui, 2))
    end
    g2 = zi -> begin
        Xi, Ui = zi
        sum((Xi[1, k] - goal2[1])^2 + (Xi[2, k] - goal2[2])^2 for k in 1:size(Xi, 2)) +
        0.1 * sum(Ui[j, k]^2 for j in 1:size(Ui, 1), k in 1:size(Ui, 2))
    end

    return LexicographicGameProblem([p1, p2], col_pairs, Function[g1, g2], T(tf), T(dt))
end

# ============================================================================
@testset "LIBR solver" begin

    # ── Test 1: solver_capabilities ──────────────────────────────────────────
    @testset "solver_capabilities" begin
        caps = solver_capabilities(LIBR)
        @test :LexicographicGame in caps
        @test :OpenLoopPolicies   in caps
        @test :SeparableDynamics  in caps
        @test :DiscreteTime       in caps
    end

    # ── Test 2: LIBR struct defaults and keyword construction ────────────────
    @testset "LIBR struct" begin
        s = LIBR()
        @test s.max_iter   == 100
        @test s.inner_iter == 200
        @test s.col_penalty > 0
        @test s.col_slack   >= 0

        s2 = LIBR(max_iter=10, inner_iter=5, step_size=0.001)
        @test s2.max_iter   == 10
        @test s2.inner_iter == 5
        @test s2.step_size  ≈ 0.001
    end

    # ── Test 3: solve returns GNEPSolution with correct structure ─────────────
    @testset "solve returns GNEPSolution" begin
        game   = _build_2agent_lex_game(tf=1.0, dt=0.2)
        solver = LIBR(max_iter=2, inner_iter=3, step_size=0.001)
        sol    = solve(game, solver; verbose=false)

        @test sol isa GNEPSolution
        @test length(sol.trajectories) == 2
        @test sol.trajectories[1].player_id == 1
        @test sol.trajectories[2].player_id == 2
        @test sol.strategy isa OpenLoopStrategy
        @test haskey(sol.solver_info, :J_col)
        @test haskey(sol.solver_info, :J_per)
        @test haskey(sol.solver_info, :Δ_hist)
    end

    # ── Test 4: trajectory dimensions ────────────────────────────────────────
    @testset "trajectory dimensions" begin
        N      = 10
        game   = _build_2agent_lex_game(tf=N*0.1, dt=0.1)
        solver = LIBR(max_iter=2, inner_iter=2, step_size=0.001)
        sol    = solve(game, solver)

        for i in 1:2
            tr = sol.trajectories[i]
            @test size(tr.states,   1) == 4      # private state dim
            @test size(tr.states,   2) == N + 1  # N+1 columns
            @test size(tr.controls, 1) == 2      # private control dim
            @test size(tr.controls, 2) == N      # N columns
            @test length(tr.times)     == N + 1
        end
    end

    # ── Test 5: collision costs are non-negative ──────────────────────────────
    @testset "J_col non-negative" begin
        game   = _build_2agent_lex_game()
        solver = LIBR(max_iter=3, inner_iter=5, step_size=0.005)
        sol    = solve(game, solver)

        J_col = sol.solver_info[:J_col]
        @test length(J_col) == 2
        @test all(J_col .>= -1e-10)
    end

    # ── Test 6: convergence flag and iteration count ──────────────────────────
    @testset "convergence flag and history" begin
        game   = _build_2agent_lex_game(tf=1.0, dt=0.2)
        solver = LIBR(max_iter=4, inner_iter=3, step_size=0.001)
        sol    = solve(game, solver)

        @test sol.iterations ≤ 4
        Δ = sol.solver_info[:Δ_hist]
        @test length(Δ) == sol.iterations
        @test all(Δ .>= 0)
    end

    # ── Test 7: determinism — two identical solves give identical results ─────
    @testset "determinism" begin
        game   = _build_2agent_lex_game(tf=1.0, dt=0.2)
        solver = LIBR(max_iter=3, inner_iter=4, step_size=0.001)
        sol1   = solve(game, solver)
        sol2   = solve(game, solver)

        for i in 1:2
            @test sol1.trajectories[i].states   ≈ sol2.trajectories[i].states
            @test sol1.trajectories[i].controls ≈ sol2.trajectories[i].controls
        end
    end

    # ── Test 8: zero-init trajectory satisfies dynamics ───────────────────────
    #   With zero controls and double-integrator dynamics, states should be const.
    @testset "zero-control initial state" begin
        game   = _build_2agent_lex_game(tf=1.0, dt=0.1)
        solver = LIBR(max_iter=1, inner_iter=0, step_size=0.0)  # no gradient steps
        sol    = solve(game, solver)

        # With zero controls and zero initial velocity, positions should be const
        X1 = sol.trajectories[1].states
        @test all(X1[1, :] .≈ 1.0)   # px = 1.0 throughout
        @test all(X1[2, :] .≈ 0.0)   # py = 0.0 throughout
        @test all(X1[3, :] .≈ 0.0)   # vx = 0.0 throughout
        @test all(X1[4, :] .≈ 0.0)   # vy = 0.0 throughout
    end

    # ── Test 9: Miller & Mitra 4-agent smoke test ─────────────────────────────
    @testset "Miller & Mitra 4-agent smoke" begin
        include(joinpath(@__DIR__, "..", "examples", "miller_mitra_example.jl"))
        sol = solve_miller_mitra(max_iter=2, inner_iter=3, step_size=0.001)

        @test sol isa GNEPSolution
        @test length(sol.trajectories) == 4
        J_col = sol.solver_info[:J_col]
        @test length(J_col) == 4
        @test all(J_col .>= -1e-10)
    end

end # @testset "LIBR solver"
