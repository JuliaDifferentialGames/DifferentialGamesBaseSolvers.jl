# ============================================================================
# invlq.jl — Inverse Linear-Quadratic Differential Game Solver
#
# Implements the algebraic solution-set algorithm of:
#   Inga, J., Bischoff, E., Molloy, T.L., Flad, M., Hohmann, S. (2019).
#   "Solution sets for inverse non-cooperative linear-quadratic differential
#   games."  IEEE Control Systems Letters, 3(4), 871–876.
#   DOI: 10.1109/LCSYS.2019.2919271
#
# Problem (Section II):
#   System:  ẋ(t) = Ax(t) + Σᵢ Bᵢuᵢ(t)
#   Control: uᵢ(t) = −Kᵢx(t)   (linear feedback Nash equilibrium)
#   Cost:    Jᵢ = ½ ∫₀^∞ (xᵀQᵢx + Σⱼ uⱼᵀ Rᵢⱼ uⱼ) dt
#
# Given K* = (K₁*,…,Kₙ*), find ALL (Qᵢ, Rᵢ₁,…,RᵢN) consistent with Nash.
#
# Key equations:
#   F      = A − Σᵢ BᵢKᵢ                  (closed-loop matrix)          (3)
#   F⊕     = F^T ⊗ Iₙ + Iₙ ⊗ F^T         (Kronecker sum)               (9)
#   Sᵢ     = (Iₙ ⊗ Bᵢ^T) F⊕^{−1}         (n·pᵢ × n²)                  (10)
#   Kⱼ^⊗   = Kⱼ^T ⊗ Kⱼ^T                 (n² × pⱼ²)                   (11)
#   θᵢ     = [vec(Qᵢ); vec(Rᵢ₁);…;vec(RᵢN)]                            (12)
#   Mᵢ θᵢ = 0  ⟹ θᵢ ∈ ker(Mᵢ)          (canonical parameter set)      (13)
#   Mᵢ = [Sᵢ | SᵢK₁^⊗ | ⋯ | (SᵢKᵢ^⊗ + Kᵢ^T ⊗ Ipᵢ) | ⋯ | SᵢKₙ^⊗]  (14)
#   Θ  = ∪ᵢ ker(Mᵢ)   with Rᵢᵢ ≻ 0                                     (19)
#
# If K* is not known, it is estimated from observed trajectories via LS (22):
#   K̂ᵢ = argmin Σₖ ‖Kᵢ x^[k] + uᵢ^[k]‖₂²
# ============================================================================

using LinearAlgebra
import DifferentialGamesBase: InverseLQGame, InverseLQGameSolution,
                               has_exact_K, has_trajectory_data,
                               extract_cost_matrices

# ============================================================================
# Solver struct
# ============================================================================

"""
    InverseLQGames(; tol, K_ridge, svd_tol)

Algebraic solver for inverse infinite-horizon LQ differential games
(Inga et al. 2019, IEEE Control Systems Letters).

Computes the canonical parameter set
    Θ = ∩ᵢ ker(Mᵢ)  with Rᵢᵢ ≻ 0
characterising **all** cost-function matrices consistent with an observed
Nash equilibrium (Theorem 2).

# Algorithm steps (exact K* mode)
1. Compute closed-loop matrix `F = A − Σᵢ BᵢKᵢ*`  (eq. 3)
2. Form Kronecker sum `F⊕ = F^T⊗Iₙ + Iₙ⊗F^T` and invert  (eq. 9)
3. For each player i compute `Sᵢ = (Iₙ⊗Bᵢ^T) F⊕⁻¹`  (eq. 10)
4. Form `Mᵢ` via Kronecker products of feedback gains  (eq. 14)
5. Compute `ker(Mᵢ)` via SVD — columns of `V` for near-zero singular values

# Algorithm steps (trajectory mode)
0. Estimate K̂ᵢ from observed (X, Uᵢ) via regularised least squares (eq. 22)
Then proceed as above with K̂ᵢ in place of Kᵢ*.

# Options
- `tol::Float64 = 1e-6`    : residual threshold; `converged` is `true` when
                              all `‖Mᵢθᵢ‖ ≤ tol`
- `K_ridge::Float64 = 0.0` : Tikhonov regularisation λ for K estimation;
                              set > 0 for noisy trajectory data
- `svd_tol::Float64 = 1e-8`: relative singular-value threshold for null-space
                              computation; a singular value σ is treated as zero
                              when `σ < svd_tol * σ_max`
"""
@kwdef struct InverseLQGames <: GameSolver
    tol     ::Float64 = 1e-6
    K_ridge ::Float64 = 0.0
    svd_tol ::Float64 = 1e-8
end

solver_capabilities(::Type{InverseLQGames}) = [
    :InverseLQGame,
    :LQGame,
    :FeedbackPolicies,
    :InfiniteHorizon,
]

# ============================================================================
# Internal helpers
# ============================================================================

"""
    _ilq_estimate_K(X, U_all, control_dims, ridge) -> Vector{Matrix{T}}

Estimate feedback matrices K̂ = [K̂₁,…,K̂ₙ] from trajectory samples.

Solves (per player i, eq. 22 of Inga et al. 2019):
    K̂ᵢ = argmin Σₖ ‖Kᵢ x^[k] + uᵢ^[k]‖₂²

which has the closed-form solution
    K̂ᵢ = −Uᵢ Xᵀ (X Xᵀ + λI)⁻¹,   λ = ridge

# Arguments
- `X`    : n × K_samples state matrix
- `U_all`: Vector of pᵢ × K_samples control matrices
- `ridge`: Tikhonov regularisation λ ≥ 0
"""
function _ilq_estimate_K(
    X           ::Matrix{T},
    U_all       ::Vector{Matrix{T}},
    ridge       ::Float64
) where {T}
    n, K_samples = size(X)
    np = length(U_all)

    # Gram matrix G = X Xᵀ + λIₙ  (n × n)
    G = X * X'
    if ridge > 0
        G += T(ridge) * Matrix{T}(I, n, n)
    end
    G_inv = inv(Symmetric(G))      # symmetric: use Symmetric wrapper for stability

    # For each player: K̂ᵢ = −Uᵢ Xᵀ G⁻¹
    map(U_all) do Ui
        -Ui * X' * G_inv
    end
end

"""
    _ilq_check_hurwitz(F) -> (is_hurwitz, eigenvalues)

Check that the closed-loop matrix F is Hurwitz (all eigenvalues have
negative real parts), which is required for the Kronecker sum F⊕ to be
invertible (Lemma 1, Inga et al. 2019).
"""
function _ilq_check_hurwitz(F::Matrix{T}) where {T}
    ev = eigvals(F)
    return all(real(λ) < 0 for λ in ev), ev
end

"""
    _ilq_build_Mi(Si, K, i, control_dims) -> Matrix{T}

Build constraint matrix Mᵢ ∈ ℝ^{n·pᵢ × L} (eq. 14).

Columns are ordered as:
  [Sᵢ | SᵢK₁^⊗ | ⋯ | (SᵢKᵢ^⊗ + Kᵢ^T⊗Ipᵢ) | ⋯ | SᵢKₙ^⊗]

mapping θᵢ = [vec(Qᵢ); vec(Rᵢ₁); …; vec(RᵢN)] into the residual space.
"""
function _ilq_build_Mi(
    Si          ::Matrix{T},
    K           ::Vector{Matrix{T}},
    i           ::Int,
    control_dims::Vector{Int}
) where {T}
    np = length(K)
    pᵢ = control_dims[i]

    # Pre-compute Kⱼ^⊗ = Kⱼ^T ⊗ Kⱼ^T  for all j  (n² × pⱼ²)
    K_kron = [kron(K[j]', K[j]') for j in 1:np]

    # Column blocks: first block is Sᵢ (for vec(Qᵢ))
    blocks = Matrix{T}[Si]

    # Remaining blocks: one per j = 1:N (for vec(Rᵢⱼ))
    Ipᵢ = Matrix{T}(I, pᵢ, pᵢ)
    for j in 1:np
        if j == i
            # (SᵢKᵢ^⊗ + Kᵢ^T ⊗ Ipᵢ)  ∈ ℝ^{n·pᵢ × pᵢ²}
            push!(blocks, Si * K_kron[i] + kron(K[i]', Ipᵢ))
        else
            # SᵢKⱼ^⊗  ∈ ℝ^{n·pᵢ × pⱼ²}
            push!(blocks, Si * K_kron[j])
        end
    end

    return hcat(blocks...)
end

"""
    _ilq_nullspace(Mi, svd_tol) -> (kernel_basis, kernel_dim, singular_values)

Compute an orthonormal basis for ker(Mᵢ) via full SVD.

Uses `svd(Mi, full=true)` so that V is (nc × nc) and all right singular
vectors are available — including those beyond the rank of Mᵢ that span
the null space.  A singular value σ is treated as zero when σ < svd_tol * σ_max.

For a matrix Mᵢ with nr rows and nc > nr columns, the last (nc − rank) columns
of V are always in the null space regardless of the singular-value threshold.
"""
function _ilq_nullspace(Mi::Matrix{T}, svd_tol::Float64) where {T}
    nr, nc = size(Mi)

    # Full SVD: V is nc × nc; S has length min(nr, nc) = nr
    sv = svd(Mi; full=true)          # U:(nr×nr), S:(nr,), Vt:(nr×nc), V:(nc×nc)
    σ_max = isempty(sv.S) ? T(1) : sv.S[1]

    # Indices of singular values below threshold (numerically zero)
    small_sv_idx = findall(σ -> σ < T(svd_tol) * σ_max, sv.S)

    # The null space basis consists of:
    #   (a) columns of V corresponding to small singular values
    #   (b) all columns of V beyond index nr (structurally null — no singular value)
    structural_null_start = nr + 1
    null_col_idx = unique(vcat(small_sv_idx, structural_null_start:nc))

    if isempty(null_col_idx)
        return zeros(T, nc, 0), 0, sv.S
    end

    # Columns of the full V give the null-space basis
    kern = sv.V[:, null_col_idx]
    return kern, length(null_col_idx), sv.S
end

# ============================================================================
# Main solve dispatch
# ============================================================================

"""
    solve(prob::InverseLQGame{T}, solver::InverseLQGames; verbose=false)
        -> InverseLQGameSolution{T}

Compute the canonical parameter set Θ = ∩ᵢ ker(Mᵢ) (Theorem 2).

# Keyword Arguments
- `verbose::Bool = false` : print progress and diagnostic information

# Returns
[`InverseLQGameSolution{T}`](@ref) containing per-player constraint matrices,
null-space bases, minimum-norm representatives, and diagnostics.

# Errors / Warnings
- Warns if F is not Hurwitz (Kronecker sum may be singular)
- Warns if any Mᵢ has an empty null space (no consistent parameters found)
"""
function solve(
    prob   ::InverseLQGame{T},
    solver ::InverseLQGames;
    verbose::Bool = false
) where {T}
    np = prob.n_players
    n  = prob.state_dim
    ps = prob.control_dims
    A  = prob.A
    B  = prob.B

    t_start = time()

    # ── Step 1: obtain Nash feedback matrices ─────────────────────────────────
    K_est_residuals = nothing
    if has_exact_K(prob)
        K = prob.K_star
        verbose && println("[InverseLQGames] Using exact K* ($(np) players)")
    elseif has_trajectory_data(prob)
        X_traj = prob.state_trajectories
        U_traj = prob.control_trajectories
        verbose && print("[InverseLQGames] Estimating K from trajectories (ridge=$(solver.K_ridge))...")
        K = _ilq_estimate_K(X_traj, U_traj, solver.K_ridge)
        # LS residuals
        K_est_residuals = map(1:np) do i
            norm(K[i] * X_traj + U_traj[i]) / sqrt(size(X_traj, 2))
        end
        verbose && println(" done. LS residuals: $(round.(K_est_residuals; sigdigits=3))")
    else
        error("InverseLQGame has neither K_star nor trajectory data. " *
              "Provide K_star or (state_trajectories, control_trajectories).")
    end

    # ── Step 2: closed-loop matrix F = A − Σᵢ BᵢKᵢ ──────────────────────────
    F = copy(A)
    for i in 1:np
        F -= B[i] * K[i]
    end
    verbose && println("[InverseLQGames] F computed. Checking Hurwitz stability...")

    is_hurwitz, F_ev = _ilq_check_hurwitz(F)
    if !is_hurwitz
        @warn "Closed-loop matrix F is not Hurwitz — the Kronecker sum F⊕ may be " *
              "singular and the inverse may be ill-conditioned. " *
              "Check that K* corresponds to a stabilizing Nash equilibrium. " *
              "Eigenvalues: $(round.(F_ev; sigdigits=4))"
    end
    verbose && println("[InverseLQGames] F eigenvalues: $(round.(F_ev; sigdigits=4))")

    # ── Step 3: Kronecker sum F⊕ = F^T ⊗ Iₙ + Iₙ ⊗ F^T  (n² × n²) ──────────
    Iₙ    = Matrix{T}(I, n, n)
    F_T   = F'
    F_pls = kron(F_T, Iₙ) + kron(Iₙ, F_T)     # n² × n²

    # Invert via LU (F_pls is guaranteed invertible when F is Hurwitz)
    F_pls_inv = inv(F_pls)

    verbose && println("[InverseLQGames] F⊕ inverted (condition number: ",
                       round(cond(F_pls); sigdigits=3), ")")

    # ── Step 4 & 5: for each player build Sᵢ then Mᵢ, then compute ker(Mᵢ) ──
    M_vec          = Vector{Matrix{T}}(undef, np)
    kernels        = Vector{Matrix{T}}(undef, np)
    kernel_dims    = Vector{Int}(undef, np)
    theta_min_norm = Vector{Vector{T}}(undef, np)
    residuals      = Vector{T}(undef, np)
    sv_all         = Vector{Vector{T}}(undef, np)

    for i in 1:np
        pᵢ = ps[i]

        # Sᵢ = (Iₙ ⊗ Bᵢ^T) F⊕⁻¹  ∈ ℝ^{n·pᵢ × n²}
        # Iₙ ⊗ Bᵢ^T: (n × n) ⊗ (pᵢ × n) → (n·pᵢ × n²)
        IkBT = kron(Iₙ, B[i]')          # n·pᵢ × n²
        Si   = IkBT * F_pls_inv          # n·pᵢ × n²

        verbose && println("[InverseLQGames] Player $i: Sᵢ computed ($(n*pᵢ)×$(n^2))")

        # Mᵢ  ∈ ℝ^{n·pᵢ × L}  where L = n² + Σⱼ pⱼ²
        Mi = _ilq_build_Mi(Si, K, i, ps)
        M_vec[i] = Mi

        verbose && println("[InverseLQGames] Player $i: Mᵢ built ($(size(Mi,1))×$(size(Mi,2)))")

        # ker(Mᵢ) via SVD
        kern, kdim, σ = _ilq_nullspace(Mi, solver.svd_tol)
        kernels[i]     = kern
        kernel_dims[i] = kdim
        sv_all[i]      = σ

        verbose && println("[InverseLQGames] Player $i: ker(Mᵢ) dim=$(kdim), " *
                           "σ_min=$(round(minimum(σ); sigdigits=4))")

        # Minimum-norm representative: take first basis vector (or zero)
        if kdim == 0
            @warn "Player $i: ker(Mᵢ) is empty — no consistent cost parameters found. " *
                  "Check that K* is an exact Nash equilibrium."
            θ_mn = zeros(T, size(Mi, 2))
        else
            θ_mn = kern[:, 1]    # first orthonormal basis vector
        end
        theta_min_norm[i] = θ_mn
        residuals[i]      = norm(Mi * θ_mn)
    end

    converged = all(r < T(solver.tol) for r in residuals)

    if verbose
        println("[InverseLQGames] Done in $(round(time()-t_start; digits=3))s.")
        println("[InverseLQGames] Converged: $converged")
        for i in 1:np
            println("  Player $i: ker dim=$(kernel_dims[i]), residual=$(round(residuals[i]; sigdigits=4))")
        end
    end

    info = Dict{Symbol, Any}(
        :singular_values        => sv_all,
        :F_eigenvalues          => F_ev,
        :is_hurwitz             => is_hurwitz,
        :K_estimated_residuals  => K_est_residuals,
        :solve_time             => time() - t_start,
    )

    return InverseLQGameSolution{T}(
        prob,
        K,
        F,
        M_vec,
        kernels,
        kernel_dims,
        theta_min_norm,
        residuals,
        converged,
        info
    )
end

# ============================================================================
# Convenience: solve from a GameProblem + observed K* (LQ game round-trip)
# ============================================================================

"""
    InverseLQGame(game::GameProblem{T}, K_star) -> InverseLQGame{T}

Construct an `InverseLQGame` directly from a finite-horizon `GameProblem`.
Extracts `A` and the `B` matrices from `game.dynamics` (must be
`LinearDynamics`). Useful for testing the round-trip: solve forward → invert.

# Example
```julia
fwd  = LQGameProblem(A, B, Q, R, Qf, x0, 10.0; dt=0.01)
sol  = solve(fwd, FNELQ())
K_fwd = [sol.strategy.gains[i][:, :, 1] for i in 1:2]   # feedback gains at t=0
iprob = InverseLQGame(fwd, K_fwd)
isol  = solve(iprob, InverseLQGames())
```
"""
function InverseLQGame(
    game  ::GameProblem{T},
    K_star::Vector{Matrix{T}}
) where {T}
    dyn = game.dynamics
    @assert dyn isa LinearDynamics "InverseLQGame(game, K_star) requires LinearDynamics"
    A = get_A(dyn, 1)
    B = [get_B(dyn, i, 1) for i in 1:game.n_players]
    return InverseLQGame(A, B, K_star)
end
