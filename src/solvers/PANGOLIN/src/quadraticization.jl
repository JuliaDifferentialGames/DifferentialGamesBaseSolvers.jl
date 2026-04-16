# # ============================================================================
# # quadraticization.jl
# #
# # Phase 2 — Quadraticization Layer for iLQGames++
# #
# # Provides:
# #   - quadraticize(cost, x_ref, u_ref, k)  → LQ coefficient matrices
# #   - linearize(dyn, x_ref, u_ref, k)      → (A_k, B_vecs_k, d_k)
# #   - LTVLQApproximation                    → pre-allocated solver scratch
# #   - OperatingPoint                        → current trajectory iterate
# #   - lq_approximation!(lqg, game, op)      → orchestrator
# #
# # Design principles:
# #   - quadraticize/linearize are pure functions (no mutation); the orchestrator
# #     lq_approximation! owns all mutation via LTVLQApproximation.
# #   - LTVLQApproximation holds raw Matrix{T} arrays, not LQStageCost wrappers,
# #     to keep the iLQGames inner loop allocation-free.
# #   - ForwardDiff AD is gated behind GeneralStageCost and CoupledNonlinearDynamics;
# #     LinearDynamics and LQStageCost paths are purely algebraic (zero AD overhead).
# # ============================================================================

# # ============================================================================
# # OperatingPoint — current trajectory iterate
# # ============================================================================

# """
#     OperatingPoint{T}

# Current trajectory iterate for the iLQGames outer loop.

# Holds the reference trajectory around which the LQ approximation is built
# at each iteration. Written by `rollout!`; read by `lq_approximation!`.

# # Fields
# - `x` : State trajectory; `x[k]` is the state at step k, length N+1
# - `u` : Per-player control sequences; `u[i][k]` is player i's control at step k, length N

# # Indexing convention
# - `x[1]` = x(0) = initial state
# - `x[N+1]` = x(N) = terminal state
# - `u[i][k]` = uᵢ(k-1), control applied at step k-1 producing x[k]
# """
# struct OperatingPoint{T}
#     x::Vector{Vector{T}}          # length N+1
#     u::Vector{Vector{Vector{T}}}  # n_players × N

#     function OperatingPoint{T}(
#         x::Vector{Vector{T}},
#         u::Vector{Vector{Vector{T}}}
#     ) where {T}
#         N = length(x) - 1
#         @assert N > 0 "Operating point must have at least one timestep"
#         n_players = length(u)
#         @assert n_players > 0 "Must have at least one player"
#         for i in 1:n_players
#             @assert length(u[i]) == N "u[$i] must have length N=$N"
#         end
#         new{T}(x, u)
#     end
# end

# """
#     OperatingPoint(game::GameProblem{T}) -> OperatingPoint{T}

# Allocate a zero-initialized operating point consistent with `game` dimensions.
# Used for cold-start initialization.
# """
# function OperatingPoint(game::GameProblem{T}) where {T}
#     N         = n_steps(game)
#     n         = total_state_dim(game.dynamics)
#     n_p       = game.n_players
#     ctrl_dims = game.metadata.control_dims   # works for all dynamics types

#     x = [zeros(T, n) for _ in 1:N+1]
#     x[1] .= game.initial_state

#     u = [[zeros(T, ctrl_dims[i]) for _ in 1:N] for i in 1:n_p]
#     OperatingPoint{T}(x, u)
# end

# """
#     OperatingPoint(trajectories::Vector{<:Trajectory{T}}, n::Int) -> OperatingPoint{T}

# Convert a `Vector{Trajectory}` (from a prior solve) into an `OperatingPoint`.
# `n` is the total state dimension (shared state space).
# """
# function OperatingPoint(trajectories::Vector{<:Trajectory{T}}, n::Int) where {T}
#     N        = size(trajectories[1].states, 2) - 1
#     n_p      = length(trajectories)

#     # Shared state — all players reference the same x trajectory
#     x = [trajectories[1].states[:, k] for k in 1:N+1]
#     u = [[trajectories[i].controls[:, k] for k in 1:N] for i in 1:n_p]
#     OperatingPoint{T}(x, u)
# end

# # ============================================================================
# # LTVLQApproximation — pre-allocated solver scratch
# # ============================================================================

# """
#     LTVLQApproximation{T}

# Pre-allocated scratch storage for the LTV LQ approximation of a game.

# Built once at solver construction time (matching iLQGames.jl `_lq_mem` pattern)
# and written into in-place by `lq_approximation!` at each outer iteration.
# All fields are raw `Matrix{T}` / `Vector{T}` arrays — no `LQStageCost` wrapping —
# to keep the iLQGames backward pass allocation-free.

# # Fields (all length-N sequences unless noted)
# - `A`    : Linearized system matrices A(k), each (n × n)
# - `B`    : Per-player linearized control matrices B_i(k), each (n × mᵢ)
# - `d`    : Affine dynamics defect d(k), each (n,); zero for linear dynamics
# - `Q`    : Quadraticized state cost Q_i(k), each (n × n) per player
# - `R`    : Quadraticized control cost R_i(k), each (mᵢ × mᵢ) per player
# - `M`    : Cross-term M_i(k), each (n × mᵢ) per player
# - `q`    : Linear state cost q_i(k), each (n,) per player
# - `r`    : Linear control cost r_i(k), each (mᵢ,) per player
# - `Qf`   : Terminal cost matrices Qf_i, each (n × n) per player (length n_players)
# - `qf`   : Linear terminal cost qf_i, each (n,) per player

# # Indexing
# `A[k]`, `B[i][k]`, `Q[i][k]`, etc. — 1-based, k ∈ 1:N.
# """
# struct LTVLQApproximation{T}
#     A::Vector{Matrix{T}}                  # [N]   (n × n)
#     B::Vector{Vector{Matrix{T}}}          # [n_players][N]  (n × mᵢ)
#     d::Vector{Vector{T}}                  # [N]   (n,)
#     Q::Vector{Vector{Matrix{T}}}          # [n_players][N]  (n × n)
#     R::Vector{Vector{Matrix{T}}}          # [n_players][N]  (mᵢ × mᵢ)
#     M::Vector{Vector{Matrix{T}}}          # [n_players][N]  (n × mᵢ)
#     q::Vector{Vector{Vector{T}}}          # [n_players][N]  (n,)
#     r::Vector{Vector{Vector{T}}}          # [n_players][N]  (mᵢ,)
#     Qf::Vector{Matrix{T}}                 # [n_players]     (n × n)
#     qf::Vector{Vector{T}}                 # [n_players]     (n,)

#     n::Int          # state dimension
#     N::Int          # number of timesteps
#     n_players::Int
#     control_dims::Vector{Int}
# end

# """
#     LTVLQApproximation(game::GameProblem{T}) -> LTVLQApproximation{T}

# Pre-allocate scratch storage sized to `game`. All matrices initialised to zero.
# Call once at solver construction; reuse across iterations via `lq_approximation!`.
# """
# function LTVLQApproximation(game::GameProblem{T}) where {T}
#     N   = n_steps(game)
#     n   = total_state_dim(game.dynamics)
#     n_p = game.n_players

#     # Use metadata control_dims — always populated regardless of dynamics type
#     # (game.dynamics.control_dims only exists on LinearDynamics/SeparableDynamics)
#     ctrl_dims = game.metadata.control_dims

#     A  = [zeros(T, n, n)                         for _ in 1:N]
#     B  = [[zeros(T, n, ctrl_dims[i]) for _ in 1:N] for i in 1:n_p]
#     d  = [zeros(T, n)                             for _ in 1:N]
#     Q  = [[zeros(T, n, n)              for _ in 1:N] for _ in 1:n_p]
#     R  = [[zeros(T, ctrl_dims[i], ctrl_dims[i]) for _ in 1:N] for i in 1:n_p]
#     M  = [[zeros(T, n, ctrl_dims[i]) for _ in 1:N] for i in 1:n_p]
#     q  = [[zeros(T, n)               for _ in 1:N] for _ in 1:n_p]
#     r  = [[zeros(T, ctrl_dims[i])    for _ in 1:N] for i in 1:n_p]
#     Qf = [zeros(T, n, n)             for _ in 1:n_p]
#     qf = [zeros(T, n)                for _ in 1:n_p]

#     LTVLQApproximation{T}(A, B, d, Q, R, M, q, r, Qf, qf, n, N, n_p, ctrl_dims)
# end

# # ============================================================================
# # GeneralStageCost — AD-based wrapper for iLQGames quadraticization
# # ============================================================================

# """
#     GeneralStageCost{F, T} <: AbstractStageCost

# Wraps any callable `func(x, u) -> scalar` for use in the iLQGames quadraticization
# pipeline. Derivatives are computed via `ForwardDiff` unless analytical overrides
# are provided.

# This type is the entry point for nonlinear costs in iLQGames++. It is distinct
# from `NonlinearStageCost` (which uses the `(x, u, p, t)` signature with explicit
# parameters and time); `GeneralStageCost` uses a closure over any parameters,
# which is more compatible with ForwardDiff's type requirements.

# # Fields
# - `func`      : Cost function `(x::Vector, u::Vector) -> scalar`; must be
#                 ForwardDiff-compatible (no mutation, no type-unstable branches)
# - `grad`      : Optional `(x, u) -> (∇ₓ, ∇ᵤ)`; if `nothing`, AD is used
# - `hess`      : Optional `(x, u) -> (Hₓₓ, Hᵤᵤ, Hₓᵤ)`; if `nothing`, AD is used
# - `n_x`       : State dimension (required for pre-allocated AD config)
# - `n_u`       : Control dimension

# # ForwardDiff compatibility requirements
# - `func` must accept `AbstractVector` inputs (not concrete `Vector{Float64}`)
# - No mutation inside `func`
# - Branches on `length(x)` are fine; branches on element values are not
# - Use `sum(x .^ 2)` not `dot(x, x)` for quadratic terms if AD dual numbers
#   cause issues (both work; the former is marginally more AD-friendly)

# # Example
# ```julia
# # Quadratic cost with obstacle avoidance
# obs = [5.0, 3.0]
# cost = GeneralStageCost(
#     (x, u) -> 0.5 * sum(x[1:2] .^ 2) + 0.5 * sum(u .^ 2) + 1.0 / max(norm(x[1:2] - obs), 0.1),
#     4, 2
# )
# ```
# """
# struct GeneralStageCost{F, GF, HF, T} <: AbstractStageCost
#     func::F
#     grad::GF    # Union{Nothing, Function}
#     hess::HF    # Union{Nothing, Function}
#     n_x::Int
#     n_u::Int
#     _T::Type{T} # numeric type, for AD config

#     function GeneralStageCost(
#         func::F,
#         n_x::Int,
#         n_u::Int;
#         grad::GF = nothing,
#         hess::HF = nothing,
#         T::Type{TT} = Float64
#     ) where {F, GF, HF, TT}
#         @assert n_x > 0 && n_u > 0
#         new{F, GF, HF, TT}(func, grad, hess, n_x, n_u, TT)
#     end
# end

# """
#     GeneralTerminalCost{F, GF, HF, T} <: AbstractTerminalCost

# Wraps any callable `func(x) -> scalar` for ForwardDiff-based quadraticization
# at the terminal timestep.
# """
# struct GeneralTerminalCost{F, GF, HF, T} <: AbstractTerminalCost
#     func::F
#     grad::GF
#     hess::HF
#     n_x::Int
#     _T::Type{T}

#     function GeneralTerminalCost(
#         func::F,
#         n_x::Int;
#         grad::GF = nothing,
#         hess::HF = nothing,
#         T::Type{TT} = Float64
#     ) where {F, GF, HF, TT}
#         @assert n_x > 0
#         new{F, GF, HF, TT}(func, grad, hess, n_x, TT)
#     end
# end

# # ============================================================================
# # 2.1 — quadraticize: cost → LQ coefficients at operating point
# # ============================================================================

# """
#     quadraticize(cost::AbstractStageCost, x_ref, u_ref, k) -> (Q, R, M, q, r, c)

# Compute the second-order Taylor expansion of `cost` at `(x_ref, u_ref)` at
# timestep `k`. Returns LQ coefficient matrices suitable for the iLQGames inner loop.

# # Returns
# - `Q  :: Matrix{T}` : (n × n) state Hessian, symmetric PSD
# - `R  :: Matrix{T}` : (m × m) control Hessian, symmetric PD
# - `M  :: Matrix{T}` : (n × m) cross term
# - `q  :: Vector{T}` : (n,) linear state term
# - `r  :: Vector{T}` : (m,) linear control term
# - `c  :: T`         : scalar constant

# # Mathematical form
# ℓ̂(x, u) ≈ ½ xᵀQx + ½ uᵀRu + xᵀMu + qᵀx + rᵀu + c

# # Dispatch
# - `LQStageCost`     : identity — returns stored matrices directly (no AD)
# - `DiagonalLQStageCost` : constructs full matrices from diagonal vectors
# - `GeneralStageCost` : ForwardDiff Hessian of `func` w.r.t. [x; u]
# - `NonlinearStageCost` : delegates to stored `hessian`/`gradient` or AD
# """
# function quadraticize end

# # ── LQStageCost (LTI and LTV) — identity, zero overhead ─────────────────────

# function quadraticize(
#     cost::LQStageCost,
#     x_ref::AbstractVector,
#     u_ref::AbstractVector,
#     k::Int
# )
#     # For LTI costs the accessor ignores k; for LTV it indexes the sequence.
#     # Returned matrices are views/aliases into the stored arrays — do not mutate.
#     Q = get_Q(cost, k)
#     R = get_R(cost, k)
#     M = get_M(cost, k)
#     q = get_q(cost, k)
#     r = get_r(cost, k)
#     c = cost.c
#     return Q, R, M, q, r, c
# end

# # ── DiagonalLQStageCost — promote to full matrices ───────────────────────────

# function quadraticize(
#     cost::DiagonalLQStageCost{T},
#     x_ref::AbstractVector,
#     u_ref::AbstractVector,
#     k::Int
# ) where {T}
#     n_x = length(cost.qx)
#     n_u = length(cost.ru)
#     Q   = Diagonal(cost.qx)     # LinearAlgebra.Diagonal — acts as Matrix in BLAS
#     R   = Diagonal(cost.ru)
#     M   = zeros(T, n_x, n_u)
#     q   = zeros(T, n_x)
#     r   = zeros(T, n_u)
#     c   = cost.c
#     return Q, R, M, q, r, c
# end

# # ── GeneralStageCost — ForwardDiff ───────────────────────────────────────────

# function quadraticize(
#     cost::GeneralStageCost{F, GF, HF, T},
#     x_ref::AbstractVector,
#     u_ref::AbstractVector,
#     k::Int
# ) where {F, GF, HF, T}
#     n_x = cost.n_x
#     n_u = cost.n_u

#     z_ref     = vcat(x_ref, u_ref)
#     f_stacked = z -> cost.func(z[1:n_x], z[n_x+1:end])

#     # ── Single-pass gradient + Hessian via DiffResults ────────────────────────
#     # Matches iLQGames.jl _quadraticize_ad pattern: one ForwardDiff.hessian! call
#     # populates both the gradient and Hessian, avoiding a redundant gradient pass.
#     if cost.hess !== nothing && cost.grad !== nothing
#         Hxx, Huu, Hxu = cost.hess(x_ref, u_ref)
#         ∇x, ∇u        = cost.grad(x_ref, u_ref)
#     elseif cost.hess !== nothing
#         # Hessian supplied analytically, gradient still needs AD
#         Hxx, Huu, Hxu = cost.hess(x_ref, u_ref)
#         diff = DiffResults.GradientResult(z_ref)
#         diff = ForwardDiff.gradient!(diff, f_stacked, z_ref)
#         ∇z   = DiffResults.gradient(diff)
#         ∇x   = ∇z[1:n_x]; ∇u = ∇z[n_x+1:end]
#     else
#         # Full AD: single hessian! call returns both gradient and Hessian
#         diff   = DiffResults.HessianResult(z_ref)
#         diff   = ForwardDiff.hessian!(diff, f_stacked, z_ref)
#         ∇z     = DiffResults.gradient(diff)
#         H_full = DiffResults.hessian(diff)
#         ∇x  = ∇z[1:n_x];            ∇u  = ∇z[n_x+1:end]
#         Hxx = H_full[1:n_x, 1:n_x]; Huu = H_full[n_x+1:end, n_x+1:end]
#         Hxu = H_full[1:n_x, n_x+1:end]
#         if cost.hess !== nothing   # override with analytical if partially supplied
#             Hxx, Huu, Hxu = cost.hess(x_ref, u_ref)
#         end
#         if cost.grad !== nothing
#             ∇x, ∇u = cost.grad(x_ref, u_ref)
#         end
#     end

#     # ── Symmetry enforcement ──────────────────────────────────────────────────
#     Q = (Hxx + Hxx') / 2
#     R = (Huu + Huu') / 2
#     M = Hxu

#     # ── R positive definiteness guard ─────────────────────────────────────────
#     λ_min_R = minimum(real.(eigvals(Symmetric(R))))
#     if λ_min_R < 0
#         ε = -λ_min_R + sqrt(eps(T))
#         @warn "quadraticize(GeneralStageCost): R not PD at k=$k (λ_min=$λ_min_R); regularizing by ε=$ε"
#         R = R + ε * I
#     end

#     # ── Linear terms in absolute coordinates ─────────────────────────────────
#     q = ∇x - Q * x_ref - M  * u_ref
#     r = ∇u - R * u_ref - M' * x_ref
#     c = T(cost.func(x_ref, u_ref))

#     return Q, R, M, q, r, c
# end

# # ── NonlinearStageCost — delegate to stored derivatives or AD ─────────────────

# function quadraticize(
#     cost::NonlinearStageCost,
#     x_ref::AbstractVector{T},
#     u_ref::AbstractVector{T},
#     k::Int
# ) where {T}
#     n_x = length(x_ref)
#     n_u = length(u_ref)
#     z_ref     = vcat(x_ref, u_ref)
#     f_stacked = z -> cost.func(z[1:n_x], z[n_x+1:end], nothing, k)

#     if cost.hessian !== nothing && cost.gradient !== nothing
#         Hxx, Huu, Hxu = cost.hessian(x_ref, u_ref, nothing, k)
#         ∇x, ∇u        = cost.gradient(x_ref, u_ref, nothing, k)
#     elseif cost.hessian !== nothing
#         Hxx, Huu, Hxu = cost.hessian(x_ref, u_ref, nothing, k)
#         diff = DiffResults.GradientResult(z_ref)
#         diff = ForwardDiff.gradient!(diff, f_stacked, z_ref)
#         ∇z   = DiffResults.gradient(diff)
#         ∇x   = ∇z[1:n_x]; ∇u = ∇z[n_x+1:end]
#     else
#         # Single DiffResults pass: gradient + Hessian together
#         diff   = DiffResults.HessianResult(z_ref)
#         diff   = ForwardDiff.hessian!(diff, f_stacked, z_ref)
#         ∇z     = DiffResults.gradient(diff)
#         H_full = DiffResults.hessian(diff)
#         ∇x  = ∇z[1:n_x];            ∇u  = ∇z[n_x+1:end]
#         Hxx = H_full[1:n_x, 1:n_x]; Huu = H_full[n_x+1:end, n_x+1:end]
#         Hxu = H_full[1:n_x, n_x+1:end]
#         if cost.gradient !== nothing
#             ∇x, ∇u = cost.gradient(x_ref, u_ref, nothing, k)
#         end
#     end

#     Q = (Hxx + Hxx') / 2
#     R = (Huu + Huu') / 2
#     M = Hxu
#     q = ∇x - Q * x_ref - M  * u_ref
#     r = ∇u - R * u_ref - M' * x_ref
#     c = T(cost.func(x_ref, u_ref, nothing, k))
#     return Q, R, M, q, r, c
# end

# # ── Terminal cost quadraticization ───────────────────────────────────────────

# """
#     quadraticize(cost::AbstractTerminalCost, x_ref) -> (Qf, qf)

# Compute second-order expansion of the terminal cost at `x_ref`.

# # Returns
# - `Qf :: Matrix{T}` : (n × n) terminal Hessian, symmetric PSD
# - `qf :: Vector{T}` : (n,) linear terminal term
# """
# function quadraticize(cost::LQTerminalCost, x_ref::AbstractVector)
#     # Identity: return stored matrices directly
#     return cost.Qf, cost.qf
# end

# function quadraticize(cost::GeneralTerminalCost{F, GF, HF, T}, x_ref::AbstractVector) where {F, GF, HF, T}
#     n_x = cost.n_x

#     if cost.hess !== nothing
#         Hxx = cost.hess(x_ref)
#     else
#         Hxx = ForwardDiff.hessian(cost.func, x_ref)
#     end
#     Qf = (Hxx + Hxx') / 2

#     if cost.grad !== nothing
#         ∇x = cost.grad(x_ref)
#     else
#         ∇x = ForwardDiff.gradient(cost.func, x_ref)
#     end
#     qf = ∇x - Qf * x_ref

#     return Qf, qf
# end

# function quadraticize(cost::NonlinearTerminalCost, x_ref::AbstractVector{T}) where {T}
#     if cost.hessian !== nothing
#         Hxx = cost.hessian(x_ref, nothing)
#     else
#         Hxx = ForwardDiff.hessian(x -> cost.func(x, nothing), x_ref)
#     end
#     Qf = (Hxx + Hxx') / 2

#     if cost.gradient !== nothing
#         ∇x = cost.gradient(x_ref, nothing)
#     else
#         ∇x = ForwardDiff.gradient(x -> cost.func(x, nothing), x_ref)
#     end
#     qf = ∇x - Qf * x_ref
#     return Qf, qf
# end

# # ============================================================================
# # 2.2 — linearize: dynamics → (A_k, B_vecs_k, d_k) at operating point
# # ============================================================================

# """
#     linearize(dyn::DynamicsSpec, x_ref, u_refs, k) -> (A_k, B_vecs_k, d_k)

# First-order linearization of dynamics at `(x_ref, u_refs)` at timestep `k`.

# # Arguments
# - `dyn`     : Dynamics specification
# - `x_ref`   : Reference state, length n
# - `u_refs`  : Per-player reference controls; `u_refs[i]` is length mᵢ
# - `k`       : Timestep index (1-based)

# # Returns
# - `A_k`      : (n × n) Jacobian ∂f/∂x at (x_ref, u_refs)
# - `B_vecs_k` : Vector of (n × mᵢ) Jacobians ∂f/∂uᵢ; one per player
# - `d_k`      : (n,) affine defect: f(x_ref, u_refs) - A_k x_ref - Σᵢ Bᵢ uᵢ_ref
#                Zero for linear dynamics; nonzero for nonlinear dynamics.

# # Defect interpretation
# The linearized dynamics are: x(k+1) ≈ A_k x(k) + Σᵢ Bᵢ(k) uᵢ(k) + d_k
# `d_k` accounts for the first-order mismatch between the nonlinear dynamics
# and their linear approximation at the reference point.
# """
# function linearize end

# # ── LinearDynamics (LTI and LTV) — identity, zero defect ────────────────────

# function linearize(
#     dyn::LinearDynamics,
#     x_ref::AbstractVector,
#     u_refs::Vector{<:AbstractVector},
#     k::Int
# )
#     A_k      = get_A(dyn, k)
#     B_vecs_k = [get_B(dyn, i, k) for i in 1:length(dyn.control_dims)]
#     d_k      = zeros(eltype(x_ref), dyn.state_dim)  # Linear dynamics: defect is zero
#     return A_k, B_vecs_k, d_k
# end

# # ── CoupledNonlinearDynamics — ForwardDiff Jacobian ──────────────────────────

# """
#     linearize(dyn::CoupledNonlinearDynamics, x_ref, u_refs, k) -> (A_k, B_vecs_k, d_k)

# Linearize nonlinear dynamics via ForwardDiff Jacobians.

# # Continuous-time warning
# If `dyn` encodes continuous-time dynamics ẋ = f(x, u), this linearization
# produces the Jacobian of ẋ, not of x(k+1). Zero-order hold (ZOH) discretization
# is required before this can be used in the discrete-time iLQGames backward pass:
#     A_d = exp(A_c · dt),  B_d = A_c⁻¹(A_d - I) B_c
# This is NOT performed here — see `linearize_zoh` for the discretized version.
# Budgeted as a separate implementation step.

# # ForwardDiff compatibility
# `dyn.func` must accept `AbstractVector` arguments (dual numbers).
# The full control vector u = [u₁; u₂; ...; uₙ] is concatenated before
# Jacobian evaluation; per-player B matrices are then extracted by column slicing.
# """
# function linearize(
#     dyn::CoupledNonlinearDynamics,
#     x_ref::AbstractVector{T},
#     u_refs::Vector{<:AbstractVector},
#     k::Int
# ) where {T}
#     n   = dyn.state_dim
#     n_p = length(u_refs)

#     # Concatenate controls into a single vector for joint Jacobian
#     u_cat  = vcat(u_refs...)
#     ctrl_dims = [length(u_refs[i]) for i in 1:n_p]
#     ctrl_ranges = [sum(ctrl_dims[1:i-1])+1 : sum(ctrl_dims[1:i]) for i in 1:n_p]

#     # Joint state-control vector z = [x; u] for one Jacobian call
#     z_ref = vcat(x_ref, u_cat)
#     n_z   = length(z_ref)

#     # Wrap dynamics to accept stacked z
#     # Note: dyn.func signature is f(x, u, p, t); pass nothing/k for p/t
#     f_stacked = z -> dyn.func(z[1:n], z[n+1:end], nothing, k)

#     # Full Jacobian [∂f/∂x | ∂f/∂u] — (n × n_z)
#     if dyn.jacobian !== nothing
#         # User supplied analytical Jacobian — preferred for performance
#         J_x, J_u = dyn.jacobian(x_ref, u_cat, nothing, k)
#         J = hcat(J_x, J_u)
#     else
#         J = ForwardDiff.jacobian(f_stacked, z_ref)   # (n × n_z)
#     end

#     A_k      = J[:, 1:n]                          # (n × n)
#     B_cat    = J[:, n+1:end]                      # (n × total_controls)
#     B_vecs_k = [B_cat[:, ctrl_ranges[i]] for i in 1:n_p]

#     # Affine defect: f(x_ref, u_ref) - A x_ref - Σ Bᵢ uᵢ_ref
#     f_ref = dyn.func(x_ref, u_cat, nothing, k)
#     d_k   = f_ref - A_k * x_ref - sum(B_vecs_k[i] * u_refs[i] for i in 1:n_p)

#     return A_k, B_vecs_k, d_k
# end

# # ── SeparableDynamics — per-player linearization ─────────────────────────────

# function linearize(
#     dyn::SeparableDynamics,
#     x_ref::AbstractVector{T},
#     u_refs::Vector{<:AbstractVector},
#     k::Int
# ) where {T}
#     n_p       = length(dyn.player_dynamics)
#     state_dims = dyn.state_dims
#     state_offs = [0; cumsum(state_dims)[1:end-1]]
#     n_total    = sum(state_dims)

#     # Block-diagonal A, per-player Bᵢ
#     A_k = zeros(T, n_total, n_total)
#     B_vecs_k = Vector{Matrix{T}}(undef, n_p)
#     d_k = zeros(T, n_total)

#     for i in 1:n_p
#         xi_ref  = x_ref[state_offs[i]+1 : state_offs[i]+state_dims[i]]
#         ui_ref  = u_refs[i]
#         ni      = state_dims[i]
#         mi      = dyn.control_dims[i]
#         fi      = dyn.player_dynamics[i]

#         # Jacobians of fᵢ(xᵢ, uᵢ) w.r.t. xᵢ and uᵢ
#         zi_ref    = vcat(xi_ref, ui_ref)
#         fi_stacked = z -> fi(z[1:ni], z[ni+1:end], nothing, k)
#         Ji         = ForwardDiff.jacobian(fi_stacked, zi_ref)

#         Ai = Ji[:, 1:ni]
#         Bi = Ji[:, ni+1:end]

#         rows = state_offs[i]+1 : state_offs[i]+ni
#         A_k[rows, rows] = Ai
#         B_vecs_k[i]     = zeros(T, n_total, mi)
#         B_vecs_k[i][rows, :] = Bi

#         fi_ref = fi(xi_ref, ui_ref, nothing, k)
#         d_k[rows] = fi_ref - Ai * xi_ref - Bi * ui_ref
#     end

#     return A_k, B_vecs_k, d_k
# end

# # ============================================================================
# # 2.3 — lq_approximation!: orchestrator
# # ============================================================================

# """
#     lq_approximation!(lqg::LTVLQApproximation, game::GameProblem, op::OperatingPoint)

# Build the LTV LQ approximation of `game` around the operating point `op`.
# Results are written in-place into the pre-allocated `lqg` scratch struct.

# Called at the start of each iLQGames outer iteration. After this call,
# `lqg` holds a complete LTV LQ game that can be solved by the FNELQ backward pass.

# # Loop structure
# For k = 1, …, N:
#   1. `linearize(dyn, x_ref[k], u_refs[:,k], k)`  → A(k), B(k), d(k)
#   2. `quadraticize(cost_i, x_ref[k], u_ref_i[k], k)` → Q_i(k), R_i(k), ...
#   3. Write results into `lqg.A[k]`, `lqg.B[i][k]`, etc.

# Terminal costs are quadraticized at `op.x[N+1]`.

# # In-place write semantics
# All matrix entries in `lqg` are overwritten; no allocation occurs if `lqg`
# was constructed from the same `game`. The calling solver owns `lqg` and is
# responsible for not aliasing it across concurrent calls.
# """
# function lq_approximation!(
#     lqg::LTVLQApproximation{T},
#     game::GameProblem{T},
#     op::OperatingPoint{T}
# ) where {T}
#     N   = lqg.N
#     n_p = lqg.n_players
#     dyn = game.dynamics

#     is_sep        = dyn isa SeparableDynamics
#     state_offsets = game.metadata.state_offsets
#     state_dims    = game.metadata.state_dims

#     for k in 1:N
#         x_k = op.x[k]
#         u_k = [op.u[i][k] for i in 1:n_p]

#         A_k, B_vecs_k, d_k = linearize(dyn, x_k, u_k, k)
#         copyto!(lqg.A[k], A_k)
#         copyto!(lqg.d[k], d_k)
#         for i in 1:n_p
#             copyto!(lqg.B[i][k], B_vecs_k[i])
#         end

#         for i in 1:n_p
#             sc   = game.objectives[i].stage_cost
#             xi_k = is_sep ? x_k[state_offsets[i]+1 : state_offsets[i]+state_dims[i]] : x_k

#             Q_k, R_k, M_k, q_k, r_k, _ = quadraticize(sc, xi_k, op.u[i][k], k)

#             if is_sep
#                 ni   = state_dims[i]
#                 rows = state_offsets[i]+1 : state_offsets[i]+ni
#                 fill!(lqg.Q[i][k], zero(T))
#                 lqg.Q[i][k][rows, rows] .= Q_k
#                 fill!(lqg.M[i][k], zero(T))
#                 lqg.M[i][k][rows, :] .= M_k
#                 fill!(lqg.q[i][k], zero(T))
#                 lqg.q[i][k][rows] .= q_k
#             else
#                 copyto!(lqg.Q[i][k], Q_k)
#                 copyto!(lqg.M[i][k], M_k)
#                 copyto!(lqg.q[i][k], q_k)
#             end
#             copyto!(lqg.R[i][k], R_k)
#             copyto!(lqg.r[i][k], r_k)
#         end
#     end

#     x_N = op.x[N+1]
#     for i in 1:n_p
#         tc   = game.objectives[i].terminal_cost
#         xi_N = is_sep ? x_N[state_offsets[i]+1 : state_offsets[i]+state_dims[i]] : x_N
#         Qf_i, qf_i = quadraticize(tc, xi_N)
#         if is_sep
#             ni   = state_dims[i]
#             rows = state_offsets[i]+1 : state_offsets[i]+ni
#             fill!(lqg.Qf[i], zero(T))
#             lqg.Qf[i][rows, rows] .= Qf_i
#             fill!(lqg.qf[i], zero(T))
#             lqg.qf[i][rows] .= qf_i
#         else
#             copyto!(lqg.Qf[i], Qf_i)
#             copyto!(lqg.qf[i], qf_i)
#         end
#     end

#     return lqg
# end
