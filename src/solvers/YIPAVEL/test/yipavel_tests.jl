# ============================================================================
# yipavel_tests.jl — Tests for the YiPavel distributed primal-dual GNE solver
#
# Run as part of the package test suite from runtests.jl, or standalone via:
#   using DifferentialGamesBaseSolvers
#   include("yipavel_tests.jl")
# ============================================================================

using Test
using LinearAlgebra
using DifferentialGamesBase
using DifferentialGamesBaseSolvers

# ── Helper: build a minimal 1-step ConvexGameProblem ─────────────────────────
#
# The "1-step static game" trick:
#   - N_steps = 1, dt = 1.0, tf = 1.0
#   - Trivial dynamics: ẋ = 0 (state is irrelevant; decision = control u)
#   - Placeholder LQ objective (Q = 0, R ≈ 0, Qf = 0)
#   - Real costs provided via cost_fns to solve()
#
# Player i's control at the single time step = the static decision variable xᵢ.

function _build_static_game(np::Int, m_i::Vector{Int};
                             shared_constraints=AbstractSharedConstraint[])
    T  = Float64
    dt = T(1.0);  tf = T(1.2)   # N = round(1.2/1.0) = 1; dt < tf ✓
    dyn = (x, u, p, t) -> zeros(T, length(x))   # ẋ = 0

    players = map(1:np) do i
        mi = m_i[i]
        Q  = zeros(T, mi, mi)
        R  = Matrix{T}(I, mi, mi) * 1e-8
        Qf = zeros(T, mi, mi)
        PlayerSpec{T}(i, mi, mi, zeros(T, mi),
                      dyn, PlayerObjective(i, LQStageCost(Q, R), LQTerminalCost(Qf)), T[])
    end

    return ConvexGameProblem(players, shared_constraints, tf, dt; assume_convex=true)
end

# ── 2-player unconstrained quadratic game ─────────────────────────────────────
#
# f₁(u₁, u₂) = ½(u₁ − 2)² + ½c(u₁ + u₂)²
# f₂(u₁, u₂) = ½(u₂ − 1)² + ½c(u₁ + u₂)²
#
# Nash equilibrium (∇ᵤᵢ fᵢ = 0):
#   (1+c)u₁ + c·u₂ = 2     [c = 0.1]
#    c·u₁ + (1+c)u₂ = 1
#
# Solution: u₁* = (1.1·2 − 0.1·1) / (1.1² − 0.1²) = 2.1/1.2 = 1.75
#           u₂* = (1.1·1 − 0.1·2) / (1.1² − 0.1²) = 0.9/1.2 = 0.75

const _C = 0.1

function _unconstrained_costs_2p()
    f1 = (s1, s2) -> 0.5*(s1[1] - 2.0)^2 + 0.5*_C*(s1[1] + s2[1])^2
    f2 = (s1, s2) -> 0.5*(s2[1] - 1.0)^2 + 0.5*_C*(s1[1] + s2[1])^2
    return Function[f1, f2]
end

# Exact unconstrained NE
const _U1_STAR_UNCONSTRAINED = 1.75
const _U2_STAR_UNCONSTRAINED = 0.75

# ── 2-player constrained quadratic game ──────────────────────────────────────
#
# Same game, with coupling: u₁ + u₂ ≤ 2.0
# At NE, u₁* + u₂* = 2.5 > 2.0, so constraint is ACTIVE.
#
# Constrained variational GNE (from KKT conditions, λ ≥ 0):
#   Player 1: (1+c)u₁ + c·u₂ − 2 + λ = 0
#   Player 2: c·u₁ + (1+c)u₂ − 1 + λ = 0
#   With u₁ + u₂ = 2 (active constraint):
#     u₁* = 1.8 − λ  AND  u₁* = 1.2 + λ  →  λ* = 0.3
#     u₁* = 1.5,  u₂* = 0.5

const _CAP_2P = 2.0
const _U1_STAR_CONSTRAINED = 1.5
const _U2_STAR_CONSTRAINED = 0.5
const _LAM_STAR_CONSTRAINED = 0.3

# ============================================================================
@testset "YiPavel solver" begin

    # ── Test 1: solver_capabilities ──────────────────────────────────────────
    @testset "solver_capabilities" begin
        caps = solver_capabilities(YiPavel)
        @test :ConvexGame       in caps
        @test :GeneralizedNash  in caps
        @test :VariationalNash  in caps
        @test :SharedConstraints in caps
    end

    # ── Test 2: struct defaults ───────────────────────────────────────────────
    @testset "YiPavel struct defaults" begin
        s = YiPavel()
        @test s.max_iter == 5000
        @test s.tol      ≈ 1e-6
        @test s.τ        > 0
        @test s.ν        > 0
        @test s.σ        > 0

        s2 = YiPavel(max_iter=100, τ=0.01, ν=0.005, σ=0.01)
        @test s2.max_iter == 100
        @test s2.τ        ≈ 0.01
        @test s2.ν        ≈ 0.005
        @test s2.σ        ≈ 0.01
    end

    # ── Test 3: solve returns GNEPSolution ────────────────────────────────────
    @testset "solve returns GNEPSolution" begin
        game   = _build_static_game(2, [1, 1])
        solver = YiPavel(max_iter=10, τ=0.05, ν=0.02, σ=0.05)
        sol    = solve(game, solver; cost_fns=_unconstrained_costs_2p(), verbose=false)

        @test sol isa GNEPSolution
        @test length(sol.trajectories) == 2
        @test sol.trajectories[1].player_id == 1
        @test sol.trajectories[2].player_id == 2
        @test sol.strategy isa OpenLoopStrategy
        @test haskey(sol.solver_info, :λ_final)
        @test haskey(sol.solver_info, :Δ_hist)
    end

    # ── Test 4: trajectory dimensions ────────────────────────────────────────
    @testset "trajectory dimensions" begin
        m = 3   # 3D decision per player
        game   = _build_static_game(2, [m, m])
        solver = YiPavel(max_iter=5)

        cost1 = (s1, s2) -> sum(s1 .^ 2) + 0.1*sum((s1 .+ s2) .^ 2)
        cost2 = (s1, s2) -> sum(s2 .^ 2) + 0.1*sum((s1 .+ s2) .^ 2)
        sol = solve(game, solver; cost_fns=Function[cost1, cost2])

        for i in 1:2
            tr = sol.trajectories[i]
            @test size(tr.states,   1) == m   # private state dim = m
            @test size(tr.states,   2) == 2   # N+1 = 2 (1-step game)
            @test size(tr.controls, 1) == m   # control dim = m
            @test size(tr.controls, 2) == 1   # N = 1 step
            @test length(tr.times)     == 2   # N+1 = 2
        end
    end

    # ── Test 5: convergence to unconstrained Nash equilibrium ─────────────────
    @testset "unconstrained NE convergence (2-player quadratic)" begin
        game   = _build_static_game(2, [1, 1])
        solver = YiPavel(max_iter=5000, tol=1e-7, τ=0.1, ν=0.05, σ=0.1)
        sol    = solve(game, solver;
                       cost_fns = _unconstrained_costs_2p(),
                       verbose  = false)

        u1 = sol.trajectories[1].controls[1, 1]
        u2 = sol.trajectories[2].controls[1, 1]

        @test sol.converged
        @test u1 ≈ _U1_STAR_UNCONSTRAINED  atol=1e-4
        @test u2 ≈ _U2_STAR_UNCONSTRAINED  atol=1e-4
    end

    # ── Test 6: convergence to constrained Nash equilibrium ───────────────────
    @testset "constrained NE convergence (u₁+u₂ ≤ 2)" begin
        T      = Float64
        game   = _build_static_game(2, [1, 1])
        solver = YiPavel(max_iter=5000, tol=1e-7, τ=0.05, ν=0.02, σ=0.05)

        # Coupling: u₁ + u₂ ≤ 2  →  coupling_A = [1 1] (1×2),  coupling_b = [2]
        coup_A = T[1.0 1.0]   # 1 × 2 matrix
        coup_b = T[_CAP_2P]

        sol = solve(game, solver;
                    cost_fns     = _unconstrained_costs_2p(),
                    coupling_A   = coup_A,
                    coupling_b   = coup_b,
                    coupling_leq = true,
                    verbose      = false)

        u1  = sol.trajectories[1].controls[1, 1]
        u2  = sol.trajectories[2].controls[1, 1]
        λ1  = sol.solver_info[:λ_final][1]
        λ2  = sol.solver_info[:λ_final][2]

        @test sol.converged
        @test u1 ≈ _U1_STAR_CONSTRAINED  atol=5e-3
        @test u2 ≈ _U2_STAR_CONSTRAINED  atol=5e-3
        # Constraint satisfied
        @test u1 + u2 ≤ _CAP_2P + 1e-3
        # Dual consensus: λ₁ ≈ λ₂ ≈ λ*  (each player's local multiplier agrees)
        @test λ1[1] ≈ _LAM_STAR_CONSTRAINED  atol=0.05
        @test λ2[1] ≈ _LAM_STAR_CONSTRAINED  atol=0.05
        @test λ1[1] ≈ λ2[1]                  atol=1e-3
    end

    # ── Test 7: box constraints respected ────────────────────────────────────
    @testset "box constraints respected" begin
        T    = Float64
        game = _build_static_game(2, [2, 2])
        solver = YiPavel(max_iter=500, τ=0.1, ν=0.05, σ=0.1)

        # Push the unconstrained solution outside [0, 1.5]^2
        f1 = (s1, s2) -> -dot([5.0, 5.0], s1)   # wants s1 as large as possible
        f2 = (s1, s2) -> -dot([5.0, 5.0], s2)

        lb = [zeros(T, 2) for _ in 1:2]
        ub = [ones(T, 2) .* 1.5 for _ in 1:2]

        sol = solve(game, solver;
                    cost_fns = Function[f1, f2],
                    lb = lb, ub = ub, verbose = false)

        for i in 1:2
            u = sol.trajectories[i].controls[:, 1]
            @test all(u .>= -1e-8)
            @test all(u .<= 1.5 + 1e-8)
        end
    end

    # ── Test 8: Δ_hist is non-negative and decreasing (on average) ───────────
    @testset "convergence history" begin
        game   = _build_static_game(2, [1, 1])
        solver = YiPavel(max_iter=200, tol=1e-7, τ=0.1, ν=0.05, σ=0.1)
        sol    = solve(game, solver; cost_fns=_unconstrained_costs_2p())

        Δ = sol.solver_info[:Δ_hist]
        @test length(Δ) == sol.iterations
        @test all(Δ .>= 0)
        # Last Δ should be at or below tolerance (if converged)
        sol.converged && @test Δ[end] ≤ solver.tol * 10
    end

    # ── Test 9: determinism ───────────────────────────────────────────────────
    @testset "determinism" begin
        game   = _build_static_game(3, [2, 2, 2])
        solver = YiPavel(max_iter=50, τ=0.05, ν=0.02, σ=0.05)

        # 3-player decoupled quadratic (easy, fast)
        f1 = (s1, s2, s3) -> sum((s1 .- [1.0, 2.0]) .^ 2)
        f2 = (s1, s2, s3) -> sum((s2 .- [3.0, 1.0]) .^ 2)
        f3 = (s1, s2, s3) -> sum((s3 .- [0.5, 1.5]) .^ 2)
        fns = Function[f1, f2, f3]

        sol1 = solve(game, solver; cost_fns=fns)
        sol2 = solve(game, solver; cost_fns=fns)

        for i in 1:3
            @test sol1.trajectories[i].controls ≈ sol2.trajectories[i].controls
        end
    end

    # ── Test 10: Cournot smoke test (Yi & Pavel Section VI) ──────────────────
    @testset "Cournot smoke test (N=5, m=2)" begin
        include(joinpath(@__DIR__, "..", "examples", "cournot_example.jl"))
        sol = solve_cournot(max_iter=3000, tol=1e-5, verbose=false)

        @test sol isa GNEPSolution
        @test length(sol.trajectories) == 5

        # Market capacity approximately satisfied per market
        total_supply = [
            sum(sol.trajectories[i].controls[k, 1] for i in 1:5)
            for k in 1:CN_NMARKETS
        ]
        @test all(total_supply .<= CN_R .+ 0.05)

        # Each company's production is non-negative and within bounds
        for i in 1:5
            xi = sol.trajectories[i].controls[:, 1]
            @test all(xi .>= -1e-6)
            @test all(xi .<= CN_THETA[i] .+ 1e-4)
        end

        # Multipliers are non-negative (dual feasibility)
        λ_all = sol.solver_info[:λ_final]
        @test all(all(λi .>= -1e-8) for λi in λ_all)
    end

end # @testset "YiPavel solver"
