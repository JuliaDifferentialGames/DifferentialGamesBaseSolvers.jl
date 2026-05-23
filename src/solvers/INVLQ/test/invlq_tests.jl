# ============================================================================
# invlq_tests.jl — Tests for the InverseLQGames solver
#
# Tests the algebraic inverse LQ game algorithm of Inga et al. (2019).
# Reference: IEEE Control Systems Letters, 3(4), 871–876.
#
# Run standalone:
#   using DifferentialGamesBase, DifferentialGamesBaseSolvers
#   include("invlq_tests.jl")
# ============================================================================

using Test
using LinearAlgebra
using DifferentialGamesBase
using DifferentialGamesBaseSolvers

# Include the canonical example (defines constants and helpers)
include(joinpath(@__DIR__, "..", "examples", "double_integrator_example.jl"))

# ── Helpers ───────────────────────────────────────────────────────────────────

"""
Check whether a vector `θ` lies in the column space of `V`
(i.e., whether it is a linear combination of V's columns).
Returns the relative projection error ‖θ − V Vᵀ θ‖ / ‖θ‖.
"""
function _in_span(V::Matrix, θ::Vector)
    isempty(V) && return norm(θ) / (norm(θ) + eps())
    θ_proj = V * (V' * θ)
    return norm(θ_proj - θ) / (norm(θ) + eps())
end

# ============================================================================
@testset "InverseLQGames solver" begin

    # ── Test 1: solver_capabilities ──────────────────────────────────────────
    @testset "solver_capabilities" begin
        caps = solver_capabilities(InverseLQGames)
        @test :InverseLQGame in caps
        @test :LQGame         in caps
        @test :FeedbackPolicies in caps
        @test :InfiniteHorizon  in caps
    end

    # ── Test 2: struct defaults ───────────────────────────────────────────────
    @testset "InverseLQGames defaults" begin
        s = InverseLQGames()
        @test s.tol     ≈ 1e-6
        @test s.K_ridge ≈ 0.0
        @test s.svd_tol ≈ 1e-8

        s2 = InverseLQGames(tol=1e-4, K_ridge=0.1, svd_tol=1e-6)
        @test s2.tol     ≈ 1e-4
        @test s2.K_ridge ≈ 0.1
        @test s2.svd_tol ≈ 1e-6
    end

    # ── Test 3: InverseLQGame problem type ───────────────────────────────────
    @testset "InverseLQGame construction" begin
        prob = InverseLQGame(DI_A, [DI_B1, DI_B2], [DI_K1_STAR, DI_K2_STAR])
        @test prob.n_players  == 2
        @test prob.state_dim  == 2
        @test prob.control_dims == [1, 1]
        @test has_exact_K(prob)
        @test !has_trajectory_data(prob)

        # Trait queries
        @test is_inverse(prob)
    end

    # ── Test 4: solve returns InverseLQGameSolution ──────────────────────────
    @testset "solve returns InverseLQGameSolution" begin
        prob   = InverseLQGame(DI_A, [DI_B1, DI_B2], [DI_K1_STAR, DI_K2_STAR])
        solver = InverseLQGames()
        sol    = solve(prob, solver)

        @test sol isa InverseLQGameSolution
        @test sol.problem === prob
        @test length(sol.K_estimated) == 2
        @test length(sol.M)           == 2
        @test length(sol.kernels)     == 2
        @test length(sol.kernel_dims) == 2
        @test length(sol.residuals)   == 2
        @test length(sol.theta_min_norm) == 2
        @test haskey(sol.solver_info, :singular_values)
        @test haskey(sol.solver_info, :F_eigenvalues)
        @test haskey(sol.solver_info, :is_hurwitz)
    end

    # ── Test 5: closed-loop matrix F is Hurwitz ───────────────────────────────
    @testset "closed-loop F is Hurwitz (paper §IV.A)" begin
        prob = InverseLQGame(DI_A, [DI_B1, DI_B2], [DI_K1_STAR, DI_K2_STAR])
        sol  = solve(prob, InverseLQGames())

        @test sol.solver_info[:is_hurwitz]
        ev = sol.solver_info[:F_eigenvalues]
        @test all(real(λ) < 0 for λ in ev)

        # F = A - B1 K1 - B2 K2  (verify manually)
        F_manual = DI_A - DI_B1 * DI_K1_STAR - DI_B2 * DI_K2_STAR
        @test sol.F ≈ F_manual  atol=1e-10
    end

    # ── Test 6: M matrices have correct dimensions ───────────────────────────
    @testset "M matrix dimensions" begin
        # n=2, p_1=p_2=1, L = n² + p_1² + p_2² = 4 + 1 + 1 = 6
        prob = InverseLQGame(DI_A, [DI_B1, DI_B2], [DI_K1_STAR, DI_K2_STAR])
        sol  = solve(prob, InverseLQGames())

        n  = prob.state_dim
        ps = prob.control_dims
        L  = n^2 + sum(p^2 for p in ps)   # = 6

        for i in 1:2
            @test size(sol.M[i], 1) == n * ps[i]   # n·pᵢ = 2·1 = 2
            @test size(sol.M[i], 2) == L            # 6
        end
    end

    # ── Test 7: ground-truth parameters are in ker(Mᵢ) ──────────────────────
    #
    # Theorem 2 guarantees Mᵢ θᵢ^true = 0 for the correct Nash equilibrium.
    # This is the core correctness test.
    @testset "ground-truth θᵢ ∈ ker(Mᵢ)" begin
        prob = InverseLQGame(DI_A, [DI_B1, DI_B2], [DI_K1_STAR, DI_K2_STAR])
        sol  = solve(prob, InverseLQGames())

        θ1_true = DI_THETA1_TRUE   # [vec(Q₁); R₁₁; R₁₂]
        θ2_true = DI_THETA2_TRUE   # [vec(Q₂); R₂₁; R₂₂]

        # Residuals must be numerically zero (exact K* → machine-precision residuals)
        @test norm(sol.M[1] * θ1_true) < 1e-10
        @test norm(sol.M[2] * θ2_true) < 1e-10

        # θ_true must lie in the span of the computed kernel basis
        err1 = _in_span(sol.kernels[1], θ1_true)
        err2 = _in_span(sol.kernels[2], θ2_true)
        @test err1 < 1e-8
        @test err2 < 1e-8
    end

    # ── Test 8: kernel is non-trivial (solution set exists) ──────────────────
    @testset "kernel is non-trivial" begin
        prob = InverseLQGame(DI_A, [DI_B1, DI_B2], [DI_K1_STAR, DI_K2_STAR])
        sol  = solve(prob, InverseLQGames())

        # Mᵢ is 2×6 with rank ≤ 2, so kernel has dimension ≥ 4
        @test sol.kernel_dims[1] >= 1
        @test sol.kernel_dims[2] >= 1

        # Kernel basis vectors are orthonormal
        for i in 1:2
            V = sol.kernels[i]
            if size(V, 2) > 0
                @test V' * V ≈ I  atol=1e-10
            end
        end
    end

    # ── Test 9: min-norm kernel element has near-zero residual ───────────────
    @testset "min-norm representative residual ≈ 0" begin
        prob = InverseLQGame(DI_A, [DI_B1, DI_B2], [DI_K1_STAR, DI_K2_STAR])
        sol  = solve(prob, InverseLQGames())

        @test sol.converged
        for i in 1:2
            @test sol.residuals[i] < 1e-6
            θ_mn = sol.theta_min_norm[i]
            @test norm(sol.M[i] * θ_mn) < 1e-6
        end
    end

    # ── Test 10: extract_cost_matrices returns correct shapes ─────────────────
    @testset "extract_cost_matrices dimensions" begin
        prob = InverseLQGame(DI_A, [DI_B1, DI_B2], [DI_K1_STAR, DI_K2_STAR])
        sol  = solve(prob, InverseLQGames())

        n = 2; ps = [1, 1]
        for i in 1:2
            Q, R_vec = extract_cost_matrices(sol, i)
            @test size(Q)       == (n, n)
            @test length(R_vec) == 2
            @test size(R_vec[1]) == (ps[1], ps[1])
            @test size(R_vec[2]) == (ps[2], ps[2])
        end
    end

    # ── Test 11: get_kernel and get_M accessor functions ─────────────────────
    @testset "accessors" begin
        prob = InverseLQGame(DI_A, [DI_B1, DI_B2], [DI_K1_STAR, DI_K2_STAR])
        sol  = solve(prob, InverseLQGames())

        @test get_kernel(sol, 1) === sol.kernels[1]
        @test get_kernel(sol, 2) === sol.kernels[2]
        @test get_M(sol, 1) === sol.M[1]
        @test get_M(sol, 2) === sol.M[2]
    end

    # ── Test 12: Corollary 1 — scaling invariance ─────────────────────────────
    #
    # Corollary 1: if θᵢ is a valid parameter set, so is κᵢ θᵢ for any κᵢ > 0.
    @testset "scaling invariance (Corollary 1)" begin
        prob = InverseLQGame(DI_A, [DI_B1, DI_B2], [DI_K1_STAR, DI_K2_STAR])
        sol  = solve(prob, InverseLQGames())

        κ = 3.7
        for i in 1:2
            scaled = κ * DI_THETA1_TRUE   # κ · θ₁^true
            res = norm(sol.M[1] * scaled)
            @test res < 1e-9    # κ θᵢ is also in ker(Mᵢ) (scales with ‖M₁θ₁^true‖ ≈ 1e-15)
        end
    end

    # ── Test 13: 3-player example — uncoupled quadratic ───────────────────────
    #
    # 3-player, 3-state game where players are dynamically decoupled.
    # Each player i controls a 1D integrator: ẋᵢ = uᵢ.
    # Joint state x = [x₁; x₂; x₃], A = 0₃ₓ₃, Bᵢ = eᵢ (i-th basis vector).
    # Cost: Jᵢ = ½ ∫ (qᵢ xᵢ² + rᵢ uᵢ²) dt
    # Nash K*ᵢ = √(qᵢ/rᵢ) eᵢ^T (scalar gain on own state).
    @testset "3-player decoupled integrators" begin
        T  = Float64
        n  = 3   # one state per player
        np = 3

        A = zeros(T, n, n)                           # no coupling in dynamics
        B = [Matrix{T}(I, n, n)[:, [i]] for i in 1:np]  # Bᵢ = eᵢ ∈ ℝ³

        # True costs: qᵢ, rᵢ
        q = [2.0, 3.0, 1.0]
        r = [1.0, 2.0, 0.5]

        # Nash gain: Kᵢ* = √(qᵢ/rᵢ) eᵢ^T ∈ ℝ^{1×n}
        K_star = [
            sqrt(q[i]/r[i]) * Matrix{T}(I, n, n)[[i], :]
            for i in 1:np
        ]

        # True parameter vectors: θᵢ = [vec(Qᵢ); R_{i1}; R_{i2}; R_{i3}]
        # Only qᵢ eᵢeᵢ^T contributes to Qᵢ; Rᵢⱼ=0 for j≠i, Rᵢᵢ=rᵢ
        theta_true = map(1:np) do i
            Qi = zeros(T, n, n); Qi[i, i] = q[i]
            Rij = [zeros(T, 1, 1) for _ in 1:np]
            Rij[i] = fill(r[i], 1, 1)
            [vec(Qi); vcat([vec(Rij[j]) for j in 1:np]...)]
        end

        prob   = InverseLQGame(A, B, K_star)
        solver = InverseLQGames()
        sol    = solve(prob, solver)

        @test sol isa InverseLQGameSolution
        @test sol.solver_info[:is_hurwitz]

        # L = n² + Σ pⱼ² = 9 + 3 = 12
        L = n^2 + np   # pⱼ = 1 for all j
        for i in 1:np
            @test size(sol.M[i], 1) == n * 1   # npᵢ = 3
            @test size(sol.M[i], 2) == L
            # Ground truth in kernel (exact K* → machine-precision residuals)
            @test norm(sol.M[i] * theta_true[i]) < 1e-10
            @test _in_span(sol.kernels[i], theta_true[i]) < 1e-8
        end
    end

    # ── Test 14: trajectory mode — K estimation from noiseless data ──────────
    @testset "trajectory mode (noiseless)" begin
        # Noiseless trajectories: xₖ is arbitrary, uᵢ = -Kᵢ xₖ exactly
        T      = Float64
        K_samp = 50    # number of samples
        rng_x  = [randn(T, 2) for _ in 1:K_samp]

        X_traj = hcat(rng_x...)                           # 2×50
        U1_traj = -DI_K1_STAR * X_traj                   # 1×50
        U2_traj = -DI_K2_STAR * X_traj                   # 1×50

        prob_traj = InverseLQGame(
            DI_A,
            [DI_B1, DI_B2],
            X_traj, [U1_traj, U2_traj]
        )
        @test !has_exact_K(prob_traj)
        @test has_trajectory_data(prob_traj)

        sol_traj = solve(prob_traj, InverseLQGames(K_ridge=1e-12))

        # Estimated K* should closely match the true K* (noiseless data)
        @test sol_traj.K_estimated[1] ≈ DI_K1_STAR  atol=1e-8
        @test sol_traj.K_estimated[2] ≈ DI_K2_STAR  atol=1e-8

        # Ground-truth parameters still in kernel (K* estimated precisely from noiseless data)
        @test norm(sol_traj.M[1] * DI_THETA1_TRUE) < 1e-8
        @test norm(sol_traj.M[2] * DI_THETA2_TRUE) < 1e-8
    end

    # ── Test 15: InverseLQGame from GameProblem convenience constructor ───────
    @testset "InverseLQGame(game, K_star)" begin
        # Build a simple LQ GameProblem
        A    = DI_A
        B    = [DI_B1, DI_B2]
        Q    = [DI_Q1, DI_Q2]
        R    = [DI_R11, DI_R22]
        Qf   = [zeros(2, 2), zeros(2, 2)]
        x0   = [1.0, 0.0]
        game = LQGameProblem(A, B, Q, R, Qf, x0, 1.0; dt=0.01)

        prob = InverseLQGame(game, [DI_K1_STAR, DI_K2_STAR])
        @test prob isa InverseLQGame{Float64}
        @test prob.n_players == 2
        @test prob.A         ≈ A
        @test prob.B[1]      ≈ B[1]
        @test prob.B[2]      ≈ B[2]
    end

end # @testset "InverseLQGames solver"
