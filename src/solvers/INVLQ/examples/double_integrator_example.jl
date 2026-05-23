# ============================================================================
# double_integrator_example.jl
#
# Reproduction of the 2-player scenario from:
#   Inga, J., Bischoff, E., Molloy, T.L., Flad, M., Hohmann, S. (2019).
#   "Solution sets for inverse non-cooperative linear-quadratic differential
#   games."  IEEE Control Systems Letters, 3(4), 871–876.
#   Section IV.A (pp. 874–875)
#
# ── Problem ──────────────────────────────────────────────────────────────────
#   System (eq. 23):
#     ẋ(t) = [0 1] x(t) + [0] u₁(t) + [0] u₂(t)
#             [0 0]        [1]          [1]
#
#   Ground truth cost functions:
#     Q₁ = diag(1, 2),   Q₂ = diag(1, 0.7)
#     R₁₁ = 1, R₁₂ = 0, R₂₁ = 0, R₂₂ = 1
#
#   Nash equilibrium (infinite-horizon, from solving coupled Riccati equations):
#     K₁* = [0.5773502691896252  1.2828983775598228]
#     K₂* = [0.5773502691896261  0.5880716343657777]
#   (Paper reports 4 d.p.: K₁* ≈ [0.5773  1.2827], K₂* ≈ [0.5774  0.5882])
#
# ── Verification ─────────────────────────────────────────────────────────────
#   The paper (Theorem 2) proves that for every player i,
#     Mᵢ θᵢ^true = 0
#   where θᵢ^true = [vec(Qᵢ); vec(Rᵢ₁); vec(Rᵢ₂)].
#
#   This example verifies that the InverseLQGames solver correctly recovers
#   ker(Mᵢ) containing the ground truth parameters.
#
# ── Usage ────────────────────────────────────────────────────────────────────
#   using DifferentialGamesBase, DifferentialGamesBaseSolvers
#   include("double_integrator_example.jl")
#   sol = solve_double_integrator_inverse(verbose=true)
#   analyse_inverse_solution(sol)
# ============================================================================

using LinearAlgebra
using DifferentialGamesBase
using DifferentialGamesBaseSolvers

# ── System matrices ───────────────────────────────────────────────────────────

const DI_A  = Float64[0 1; 0 0]           # 2×2 double integrator
const DI_B1 = Float64[0; 1][:, :]         # 2×1
const DI_B2 = Float64[0; 1][:, :]         # 2×1

# ── Ground truth cost parameters ─────────────────────────────────────────────

const DI_Q1  = diagm([1.0, 2.0])    # Q₁ = diag(1, 2)
const DI_Q2  = diagm([1.0, 0.7])    # Q₂ = diag(1, 0.7)
const DI_R11 = fill(1.0, 1, 1)      # R₁₁ = 1 (1×1 matrix)
const DI_R12 = fill(0.0, 1, 1)      # R₁₂ = 0
const DI_R21 = fill(0.0, 1, 1)      # R₂₁ = 0
const DI_R22 = fill(1.0, 1, 1)      # R₂₂ = 1

# Ground-truth parameter vectors θᵢ = [vec(Qᵢ); vec(Rᵢ₁); vec(Rᵢ₂)]
# Julia uses column-major order for vec, so vec([a b; c d]) = [a, c, b, d].
const DI_THETA1_TRUE = [vec(DI_Q1); vec(DI_R11); vec(DI_R12)]  # 6×1
const DI_THETA2_TRUE = [vec(DI_Q2); vec(DI_R21); vec(DI_R22)]  # 6×1

# ── Nash equilibrium feedback matrices (from coupled algebraic Riccati equations) ─
# Computed via Lyapunov iteration to machine precision (|Riccati residual| ≈ 2e-15).
# The paper (Section IV.A) reports the rounded 4 d.p. values
#   K₁* ≈ [0.5773  1.2827],  K₂* ≈ [0.5774  0.5882].
# We store the full-precision values so that ‖Mᵢ θᵢ^true‖ ≈ 2e-15.

const DI_K1_STAR = Float64[0.5773502691896252  1.2828983775598228]   # 1×2
const DI_K2_STAR = Float64[0.5773502691896261  0.5880716343657777]   # 1×2

# ── Kernel basis vectors from the paper (eqs. 24–25) ──────────────────────────
# With the REDUCED 3-parameter diagonal structure (not full vec), the paper gives:
#   v₁^(1) = [0.4083, 0.8165, 0.4083]
#   v₂^(1) = [0.6337, 0.4437, 0.6337]
# These correspond to the diagonal-Q, self-R-only parameterisation.
# Our full-vec implementation uses a 6-dimensional parameterisation.
const DI_V1_PAPER = [0.4083, 0.8165, 0.4083]
const DI_V2_PAPER = [0.6337, 0.4437, 0.6337]

# ── Solver wrapper ────────────────────────────────────────────────────────────

"""
    solve_double_integrator_inverse(; verbose, svd_tol, tol)
        -> InverseLQGameSolution{Float64}

Run the InverseLQGames solver on the 2-player double-integrator example
from Section IV.A of Inga et al. (2019).

Returns the full solution including per-player constraint matrices and
null-space bases.
"""
function solve_double_integrator_inverse(;
    verbose ::Bool    = false,
    svd_tol ::Float64 = 1e-8,
    tol     ::Float64 = 1e-6
)
    prob = InverseLQGame(
        DI_A,
        [DI_B1, DI_B2],
        [DI_K1_STAR, DI_K2_STAR]
    )

    solver = InverseLQGames(; tol=tol, svd_tol=svd_tol)
    return solve(prob, solver; verbose=verbose)
end

# ── Analysis helper ───────────────────────────────────────────────────────────

"""
    analyse_inverse_solution(sol; θ_true=[DI_THETA1_TRUE, DI_THETA2_TRUE])

Print a detailed analysis of the inverse solution:
- Closed-loop eigenvalues (should be negative real part)
- Residuals ‖Mᵢ θᵢ^true‖ for the ground-truth parameters
- Projection error: how well θᵢ^true is captured by ker(Mᵢ)
- Representative cost matrices from the minimum-norm kernel element
"""
function analyse_inverse_solution(
    sol    ::InverseLQGameSolution{Float64};
    θ_true ::Vector{Vector{Float64}} = [DI_THETA1_TRUE, DI_THETA2_TRUE]
)
    np = sol.problem.n_players
    println("=" ^ 60)
    println("Inverse LQ Game — Double Integrator (Inga et al. 2019)")
    println("=" ^ 60)

    # Closed-loop stability
    ev = sol.solver_info[:F_eigenvalues]
    println("\nClosed-loop eigenvalues of F = A − Σᵢ BᵢKᵢ*:")
    for λ in ev
        tag = real(λ) < 0 ? "✓ stable" : "✗ UNSTABLE"
        println("  λ = $(round(real(λ), digits=4)) + $(round(imag(λ), digits=4))i  [$tag]")
    end

    for i in 1:np
        Mi = sol.M[i]
        θt = θ_true[i]
        kern = sol.kernels[i]

        println("\n── Player $i ─────────────────────────────────────────")
        println("  Mᵢ size            : $(size(Mi,1)) × $(size(Mi,2))")
        println("  ker(Mᵢ) dimension  : $(sol.kernel_dims[i])")
        println("  Singular values    : ",
                join(["$(round(s, sigdigits=3))" for s in sol.solver_info[:singular_values][i]], ", "))

        # Residual for ground-truth parameters
        res_true = norm(Mi * θt)
        println("\n  ‖Mᵢ θᵢ^true‖         : $(round(res_true, sigdigits=4))",
                res_true < 1e-4 ? "  ✓ (in kernel)" : "  ✗ (NOT in kernel)")

        # Projection of θ_true onto ker(Mᵢ)
        if sol.kernel_dims[i] > 0
            θ_proj = kern * (kern' * θt)
            proj_err = norm(θ_proj - θt) / (norm(θt) + eps())
            println("  Projection error   : $(round(proj_err, sigdigits=4))",
                    proj_err < 1e-3 ? "  ✓ (θᵢ^true ∈ span of kernel)" : "  (partial)")
        end

        # Show representative solution (min-norm kernel element)
        println("\n  Min-norm kernel representative θᵢ^(mn):")
        θ_mn = sol.theta_min_norm[i]
        Q, R_vec = extract_cost_matrices(sol, i)
        println("    Q_$i =\n", "    ", replace(string(round.(Q, digits=4)), "\n" => "\n    "))
        for j in 1:np
            println("    R_{$i,$j} = ", round.(R_vec[j], digits=4))
        end
        println("  ‖Mᵢ θᵢ^(mn)‖ = $(round(sol.residuals[i], sigdigits=4))")
    end
    println()
end
