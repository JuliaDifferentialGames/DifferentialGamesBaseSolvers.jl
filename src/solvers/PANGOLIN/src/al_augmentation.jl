# ============================================================================
# al_augmentation.jl
#
# Phase 3 — Augmented Lagrangian Constraint Layer for iLQGames++
#
# Provides:
#   - ALOptions                          → penalty schedule parameters
#   - ALAugmentedObjective               → wraps PlayerObjective + multipliers
#   - augmented_stage_cost               → evaluates AL-penalized cost
#   - quadraticize(ALAugmentedObjective) → LQ coefficients with active-set updates
#   - update_multipliers!                → dual ascent + penalty schedule
#   - constraint_violation               → per-player violation norms for logging
#
# Design principles:
#   - λ_private[j] corresponds to private constraints for this player only.
#     Length = sum of dims of all private constraints for this player.
#   - λ_shared[j] is a single centralized vector owned by ALAugmentedObjective,
#     consistent with FALCON's AL structure. Players share the same λ_shared
#     reference; the solver is responsible for synchronization post-iteration.
#   - Active-set detection at (x, u): constraint j is active iff
#     g_j(x,u) + λ_j/ρ > 0 (inequality). Equality constraints are always active.
#   - Penalty schedule mirrors FALCON: geometric ρ ← min(ρ_max, φ·ρ), triggered
#     only when violation has not decreased by at least violation_tol.
#   - OperatingPoint is the trajectory representation throughout — no Trajectory
#     allocation in the hot loop.
# ============================================================================

# ============================================================================
# ALOptions — penalty schedule parameters
# ============================================================================

"""
    ALOptions{T}

Parameters controlling the augmented Lagrangian penalty schedule.

# Fields
- `ρ_init`        : Initial penalty weight (default 1.0)
- `ρ_max`         : Maximum penalty weight; caps geometric growth (default 1e4)
- `φ`             : Geometric growth factor φ ∈ (1, 10] (default 10.0)
- `violation_tol` : Fractional violation reduction below which ρ is increased.
                    If `‖g_new‖ > (1 - violation_tol) · ‖g_old‖`, penalty fires.
                    (default 0.25 — matches FALCON outer AL logic)
- `dual_tol`      : Convergence threshold on ‖Δλ‖_∞ (default 1e-4)

# Notes
φ = 10 is aggressive and appropriate for well-scaled problems. For poorly
scaled constraints or near-singular S matrices, reduce to φ ∈ (1, 3].
ρ_max = 1e4 is a hard cap; values above ~1e4 cause floating-point cancellation
in the ΔL computation inside the backward pass (matches FALCON ω cap lesson).
"""
struct ALOptions{T}
    ρ_init::T
    ρ_max::T
    φ::T
    violation_tol::T
    dual_tol::T

    function ALOptions(;
        ρ_init::T        = 1.0,
        ρ_max::T         = 1e4,
        φ::T             = 10.0,
        violation_tol::T = 0.25,
        dual_tol::T      = 1e-4
    ) where {T <: AbstractFloat}
        @assert ρ_init > 0     "ρ_init must be positive"
        @assert ρ_max >= ρ_init "ρ_max must be ≥ ρ_init"
        @assert 1 < φ <= 10    "φ must be in (1, 10]"
        @assert 0 < violation_tol < 1 "violation_tol must be in (0, 1)"
        @assert dual_tol > 0   "dual_tol must be positive"
        new{T}(ρ_init, ρ_max, φ, violation_tol, dual_tol)
    end
end

# Default constructor — explicit Float64 values bypass the parametric inner constructor
ALOptions() = ALOptions(;
    ρ_init        = 1.0,
    ρ_max         = 1e4,
    φ             = 10.0,
    violation_tol = 0.25,
    dual_tol      = 1e-4
)

# ============================================================================
# Constraint dimension helpers
# ============================================================================

"""
    constraint_dim(c::AbstractConstraint) -> Int

Return the output dimension of a constraint function (number of scalar inequalities).
"""
constraint_dim(c::LinearConstraint)     = size(c.A, 1)
constraint_dim(c::BoundConstraint)      = 2 * length(c.lower)
constraint_dim(c::NormConstraint)       = 1
constraint_dim(c::NonlinearConstraint)  = c.dim

constraint_dim(pc::PrivateConstraint)   = constraint_dim(pc.constraint)
constraint_dim(sc::SharedConstraint)    = constraint_dim(sc.constraint)

"""
    is_equality_constraint(c::AbstractConstraint) -> Bool

Return true if c encodes an equality constraint.
"""
is_equality_constraint(::LinearConstraint)    = false
is_equality_constraint(::BoundConstraint)     = false
is_equality_constraint(::NormConstraint)      = false
is_equality_constraint(c::NonlinearConstraint) = c.constraint_type == :equality
is_equality_constraint(pc::PrivateConstraint)  = is_equality_constraint(pc.constraint)
is_equality_constraint(sc::SharedConstraint)   = is_equality_constraint(sc.constraint)

# ============================================================================
# ALAugmentedObjective
# ============================================================================

"""
    ALAugmentedObjective{T, OBJ}

Wraps a `PlayerObjective` with Augmented Lagrangian multipliers and penalty weight.

# Mathematical Form
The AL-augmented stage cost for player i at step k is:

    ℓᵢᴬᴸ(x, uᵢ, k) = ℓᵢ(x, uᵢ, k)
        + Σⱼ∈private  [λⱼᵖ gⱼ(x,uᵢ,k) + (ρ/2) max(0, gⱼ + λⱼᵖ/ρ)²]
        + Σⱼ∈shared   [λⱼˢ gⱼ(x,u,k)  + (ρ/2) max(0, gⱼ + λⱼˢ/ρ)²]

where the outer max(0, ·) is the active-set gate for inequality constraints;
equality constraints always contribute (no gating).

# Fields
- `base`             : Underlying `PlayerObjective`
- `private_constraints` : Private constraints for this player (from GameProblem)
- `shared_constraints`  : Shared constraints involving this player
- `λ_private`        : Multipliers for private constraints; length = total private dim
- `λ_shared`         : Multipliers for shared constraints; length = total shared dim.
                        Centralized: all players reference the same underlying vector
                        via the solver's ALSolverState.
- `ρ`                : Current penalty weight (mutable — updated by schedule)
- `_priv_offsets`    : Precomputed offsets into λ_private per private constraint
- `_shared_offsets`  : Precomputed offsets into λ_shared per shared constraint
- `_priv_eq`         : Bool vector: true if private constraint j is equality
- `_shared_eq`       : Bool vector: true if shared constraint j is equality
"""
mutable struct ALAugmentedObjective{T, OBJ <: PlayerObjective}
    base::OBJ
    private_constraints::Vector{<:PrivateConstraint}
    shared_constraints::Vector{<:SharedConstraint}

    λ_private::Vector{T}    # length = Σ dim(private_j)
    λ_shared::Vector{T}     # length = Σ dim(shared_j); centralized, shared by reference

    ρ::T                    # current penalty weight

    # Precomputed metadata — set at construction, never mutated
    _priv_offsets::Vector{Int}   # start index of constraint j in λ_private
    _shared_offsets::Vector{Int} # start index of constraint j in λ_shared
    _priv_eq::Vector{Bool}
    _shared_eq::Vector{Bool}
    _priv_dims::Vector{Int}
    _shared_dims::Vector{Int}
end

"""
    ALAugmentedObjective(base, private_constraints, shared_constraints, λ_shared, opts)

Construct an `ALAugmentedObjective` for one player.

`λ_shared` must be a reference to the solver-level centralized multiplier vector —
all players for the same shared constraints should reference the same `Vector{T}`.

# Arguments
- `base`                : `PlayerObjective` for this player
- `private_constraints` : This player's private constraints (filtered from GameProblem)
- `shared_constraints`  : Shared constraints involving this player
- `λ_shared`            : Centralized shared multiplier vector (solver-owned)
- `opts`                : `ALOptions` providing ρ_init
"""
function ALAugmentedObjective(
    base::OBJ,
    private_constraints::Vector{<:PrivateConstraint},
    shared_constraints::Vector{<:SharedConstraint},
    λ_shared::Vector{T},
    opts::ALOptions{T}
) where {T, OBJ <: PlayerObjective}

    # ── Private constraint metadata ──────────────────────────────────────────
    priv_dims    = [constraint_dim(pc) for pc in private_constraints]
    priv_offsets = [1; cumsum(priv_dims)[1:end-1] .+ 1]  # 1-based start indices
    λ_private    = zeros(T, sum(priv_dims; init=0))
    priv_eq      = [is_equality_constraint(pc) for pc in private_constraints]

    # ── Shared constraint metadata ────────────────────────────────────────────
    shared_dims    = [constraint_dim(sc) for sc in shared_constraints]
    shared_offsets = [1; cumsum(shared_dims)[1:end-1] .+ 1]
    shared_eq      = [is_equality_constraint(sc) for sc in shared_constraints]

    # Validate λ_shared length
    total_shared_dim = sum(shared_dims; init=0)
    @assert length(λ_shared) == total_shared_dim "λ_shared length $(length(λ_shared)) must equal total shared constraint dim $total_shared_dim"

    ALAugmentedObjective{T, OBJ}(
        base,
        private_constraints, shared_constraints,
        λ_private, λ_shared, opts.ρ_init,
        priv_offsets, shared_offsets,
        priv_eq, shared_eq,
        priv_dims, shared_dims
    )
end

# ── Helpers: extract λ slice for constraint j ─────────────────────────────────

@inline function _priv_λ_range(obj::ALAugmentedObjective, j::Int)
    start = obj._priv_offsets[j]
    stop  = start + obj._priv_dims[j] - 1
    return start:stop
end

@inline function _shared_λ_range(obj::ALAugmentedObjective, j::Int)
    start = obj._shared_offsets[j]
    stop  = start + obj._shared_dims[j] - 1
    return start:stop
end

# ============================================================================
# augmented_stage_cost
# ============================================================================

"""
    augmented_stage_cost(obj::ALAugmentedObjective{T}, x, u, k) -> T

Evaluate the AL-augmented stage cost at (x, u, k).

    ℓᴬᴸ = ℓ_base + Σⱼ AL_term(gⱼ, λⱼ, ρ)

where for inequality constraints:
    AL_term = (ρ/2) max(0, gⱼ + λⱼ/ρ)² · dim   (the "shifted clamp" form)

and for equality constraints:
    AL_term = λⱼᵀgⱼ + (ρ/2) ‖gⱼ‖²   (always active)

# Note on shifted-clamp equivalence
The standard AL penalty (ρ/2) max(0, g + λ/ρ)² is equivalent to
λᵀg + (ρ/2) max(0, g)² when the active set is correct, but the shifted
form avoids a discontinuity in the gradient at g=0 when λ>0.
We use the shifted form throughout for gradient continuity.
"""
function augmented_stage_cost(
    obj::ALAugmentedObjective{T},
    x::AbstractVector,
    u::AbstractVector,
    k::Int
) where {T}
    # Base cost
    cost = T(evaluate_stage_cost(obj.base.stage_cost, x, u, nothing, k)) * obj.base.scaling

    ρ = obj.ρ

    # ── Private constraints ───────────────────────────────────────────────────
    for (j, pc) in enumerate(obj.private_constraints)
        g = evaluate_constraint(pc.constraint, x, u, nothing, k)  # Vector{T}
        λ = obj.λ_private[_priv_λ_range(obj, j)]

        if obj._priv_eq[j]
            # Equality: always penalized — λᵀg + (ρ/2)‖g‖²
            cost += dot(λ, g) + (ρ/2) * dot(g, g)
        else
            # Inequality: shifted clamp — (ρ/2) ‖max(0, g + λ/ρ)‖²
            shifted = g + λ / ρ
            active  = max.(shifted, zero(T))
            cost   += (ρ/2) * dot(active, active)
        end
    end

    # ── Shared constraints ────────────────────────────────────────────────────
    for (j, sc) in enumerate(obj.shared_constraints)
        g = evaluate_constraint(sc.constraint, x, u, nothing, k)
        λ = obj.λ_shared[_shared_λ_range(obj, j)]

        if obj._shared_eq[j]
            cost += dot(λ, g) + (ρ/2) * dot(g, g)
        else
            shifted = g + λ / ρ
            active  = max.(shifted, zero(T))
            cost   += (ρ/2) * dot(active, active)
        end
    end

    return cost
end

# ============================================================================
# 3.2 — AL Quadraticization
# ============================================================================

"""
    quadraticize(obj::ALAugmentedObjective, x_ref, u_ref, k) -> (Q, R, M, q, r, c)

Compute the second-order Taylor expansion of the AL-augmented cost at (x_ref, u_ref, k).

# Algorithm
1. Quadraticize the base cost: (Q₀, R₀, M₀, q₀, r₀, c₀)
2. For each active constraint j (active iff shifted_j > 0 for inequality,
   always for equality):
   a. Compute constraint Jacobian (Jˣⱼ, Jᵘⱼ) at (x_ref, u_ref)
   b. Add rank-1 Hessian update: Q += ρ · Jˣⱼᵀ Jˣⱼ, R += ρ · Jᵘⱼᵀ Jᵘⱼ, M += ρ · Jˣⱼᵀ Jᵘⱼ
   c. Add gradient update: q += μⱼ · Jˣⱼᵀ, r += μⱼ · Jᵘⱼᵀ
      where μⱼ = (λⱼ + ρ·gⱼ) for the active rows (or λⱼ + ρ·gⱼ for equality)

# Active-set gating
For inequality constraint j, row d is active iff:
    gⱼ[d](x_ref, u_ref) + λⱼ[d] / ρ > 0
Active rows contribute rank-1 updates; inactive rows are skipped entirely.
This gives the correct second-order behaviour of the shifted-clamp penalty.

# Symmetry
Q and R are symmetrized after all updates. The rank-1 Hessian updates preserve
symmetry by construction (Jᵀ J is always symmetric PSD).
"""
function quadraticize(
    obj::ALAugmentedObjective{T},
    x_ref::AbstractVector,
    u_ref::AbstractVector,
    k::Int
) where {T}
    # ── Step 1: base quadraticization ─────────────────────────────────────────
    Q, R, M, q, r, c = quadraticize(obj.base.stage_cost, x_ref, u_ref, k)

    # Copy into mutable arrays — base may return aliased views (e.g. LQStageCost identity)
    Q = Matrix{T}(Q)
    R = Matrix{T}(R)
    M = Matrix{T}(M)
    q = Vector{T}(q)
    r = Vector{T}(r)
    c = T(c)

    ρ = obj.ρ

    # ── Step 2: private constraint contributions ──────────────────────────────
    for (j, pc) in enumerate(obj.private_constraints)
        g  = evaluate_constraint(pc.constraint, x_ref, u_ref, nothing, k)
        λj = obj.λ_private[_priv_λ_range(obj, j)]
        Jx, Ju = constraint_jacobian(pc.constraint, x_ref, u_ref, nothing, k)

        _add_al_quadratic_update!(Q, R, M, q, r, g, λj, Jx, Ju, ρ, obj._priv_eq[j])
    end

    # ── Step 3: shared constraint contributions ───────────────────────────────
    for (j, sc) in enumerate(obj.shared_constraints)
        g  = evaluate_constraint(sc.constraint, x_ref, u_ref, nothing, k)
        λj = obj.λ_shared[_shared_λ_range(obj, j)]
        Jx, Ju = constraint_jacobian(sc.constraint, x_ref, u_ref, nothing, k)

        _add_al_quadratic_update!(Q, R, M, q, r, g, λj, Jx, Ju, ρ, obj._shared_eq[j])
    end

    # ── Symmetrize ────────────────────────────────────────────────────────────
    Q = (Q + Q') / 2
    R = (R + R') / 2

    # Regularize R if rank-1 updates have damaged positive definiteness
    λ_min_R = minimum(real.(eigvals(Symmetric(R))))
    if λ_min_R < 0
        ε = -λ_min_R + sqrt(eps(T))
        R = R + ε * I
    end

    return Q, R, M, q, r, c
end

"""
    _add_al_quadratic_update!(Q, R, M, q, r, g, λ, Jx, Ju, ρ, is_eq)

In-place rank-1 Hessian and gradient updates for a single constraint block.

For each row d of the constraint (one scalar inequality/equality per row):
- Active (or equality): add ρ · jˣ[d] jˣ[d]ᵀ to Q, ρ · jᵘ[d] jᵘ[d]ᵀ to R,
                        ρ · jˣ[d] jᵘ[d]ᵀ to M, μ[d] · jˣ[d] to q, μ[d] · jᵘ[d] to r
  where μ[d] = λ[d] + ρ · g[d]
- Inactive (inequality only): skip row entirely

Operates row-by-row to correctly handle partial active sets (e.g., BoundConstraint
where some bound directions are active and others are not).
"""
function _add_al_quadratic_update!(
    Q::Matrix{T}, R::Matrix{T}, M::Matrix{T},
    q::Vector{T}, r::Vector{T},
    g::Vector{T}, λ::Vector{T},
    Jx::AbstractMatrix, Ju::AbstractMatrix,
    ρ::T, is_eq::Bool
) where {T}
    n_rows = length(g)
    @assert length(λ) == n_rows
    @assert size(Jx, 1) == n_rows
    @assert size(Ju, 1) == n_rows

    for d in 1:n_rows
        # Active-set check for inequality constraints
        is_active = is_eq || (g[d] + λ[d] / ρ > zero(T))
        is_active || continue

        jx = Jx[d, :]   # (n_x,) — row d of state Jacobian
        ju = Ju[d, :]   # (n_u,) — row d of control Jacobian
        μd = λ[d] + ρ * g[d]

        # Rank-1 Hessian updates: ρ · jᵀj
        # jx * jx' is a matrix (outer product) — cannot use @. here
        Q .+= ρ .* (jx * jx')   # (n_x × n_x)
        R .+= ρ .* (ju * ju')   # (n_u × n_u)
        M .+= ρ .* (jx * ju')   # (n_x × n_u)

        # Gradient updates: μ · j — scalar × vector, broadcast is correct
        q .+= μd .* jx
        r .+= μd .* ju
    end
    return nothing
end

# ============================================================================
# 3.3 — Dual Update and Penalty Schedule
# ============================================================================

"""
    ALSolverState{T}

Solver-level state for the AL outer loop. Owns the centralized shared multiplier
vector and tracks violation history for the penalty schedule.

# Fields
- `λ_shared`          : Centralized shared multiplier vector (all players reference this)
- `ρ`                 : Current penalty weight (synchronized across all player AL objectives)
- `prev_violation`    : ‖g_shared‖_∞ from previous outer iteration (for schedule trigger)
- `opts`              : ALOptions

# Usage
One `ALSolverState` per solve. At construction, `λ_shared` is passed by reference
into each player's `ALAugmentedObjective`. After each iLQGames inner solve,
call `update_multipliers!` once on each player's AL objective (private multipliers),
then call `update_shared_multipliers!` once on the `ALSolverState` (shared multipliers),
then call `maybe_update_penalty!` to trigger ρ growth if warranted.
"""
mutable struct ALSolverState{T}
    λ_shared::Vector{T}
    ρ::T
    prev_violation::T
    opts::ALOptions{T}
    _shared_eq::Vector{Bool}       # equality flags per shared constraint block
    _shared_offsets::Vector{Int}   # offsets into λ_shared per shared constraint
    _shared_dims::Vector{Int}
end

"""
    ALSolverState(shared_constraints, opts) -> ALSolverState{T}

Construct solver-level AL state for the given shared constraints.
"""
function ALSolverState(
    shared_constraints::Vector{<:SharedConstraint},
    opts::ALOptions{T}
) where {T}
    shared_dims    = [constraint_dim(sc) for sc in shared_constraints]
    shared_offsets = [1; cumsum(shared_dims)[1:end-1] .+ 1]
    shared_eq      = [is_equality_constraint(sc) for sc in shared_constraints]
    λ_shared       = zeros(T, sum(shared_dims; init=0))

    ALSolverState{T}(
        λ_shared, opts.ρ_init, T(Inf), opts,
        shared_eq, shared_offsets, shared_dims
    )
end

# ── Private multiplier update ─────────────────────────────────────────────────

"""
    update_multipliers!(obj::ALAugmentedObjective, op::OperatingPoint)

Update private Lagrange multipliers using the AL dual ascent rule over the
trajectory encoded in `op`.

For each timestep k and each private constraint j:
  - Inequality: λⱼ[d] ← max(0, λⱼ[d] + ρ · ḡⱼ[d])
                where ḡⱼ[d] = (1/N) Σₖ gⱼ[d](x(k), uᵢ(k))
  - Equality:   λⱼ[d] ← λⱼ[d] + ρ · ḡⱼ[d]   (unclamped)

Returns the ‖Δλ‖_∞ for convergence checking.
"""
function update_multipliers!(
    obj::ALAugmentedObjective{T},
    op::OperatingPoint{T}
) where {T}
    N   = length(op.x) - 1
    ρ   = obj.ρ
    i   = obj.base.player_id
    Δλ_max = zero(T)

    for (j, pc) in enumerate(obj.private_constraints)
        rng = _priv_λ_range(obj, j)
        dim = obj._priv_dims[j]
        ḡ   = zeros(T, dim)

        # Time-average constraint violation over trajectory
        for k in 1:N
            x_k = op.x[k]
            u_k = op.u[i][k]
            ḡ .+= evaluate_constraint(pc.constraint, x_k, u_k, nothing, k)
        end
        ḡ ./= N

        λ_old = copy(obj.λ_private[rng])

        if obj._priv_eq[j]
            obj.λ_private[rng] .+= ρ .* ḡ
        else
            obj.λ_private[rng] .= max.(zero(T), obj.λ_private[rng] .+ ρ .* ḡ)
        end

        Δλ_max = max(Δλ_max, maximum(abs.(obj.λ_private[rng] .- λ_old)))
    end

    return Δλ_max
end

# ── Shared multiplier update (centralized) ────────────────────────────────────

"""
    update_shared_multipliers!(state::ALSolverState, objectives, op::OperatingPoint)

Update the centralized shared multiplier vector using violation averaged over
all players involved in each shared constraint.

For each shared constraint j with players Pⱼ:
  ḡⱼ = (1/|Pⱼ|) Σᵢ∈Pⱼ (1/N) Σₖ gⱼ(x(k), uᵢ(k))

The player average accounts for the fact that shared constraints may depend
on multiple players' controls (e.g., collision avoidance between agents).

`objectives` is the `Vector{ALAugmentedObjective}` — one per player, indexed 1:n_players.
"""
function update_shared_multipliers!(
    state::ALSolverState{T},
    shared_constraints::Vector{<:SharedConstraint},
    objectives::Vector{<:ALAugmentedObjective},
    op::OperatingPoint{T}
) where {T}
    N = length(op.x) - 1

    Δλ_max         = zero(T)
    total_violation = zero(T)

    for (j, sc) in enumerate(shared_constraints)
        rng     = state._shared_offsets[j] : state._shared_offsets[j] + state._shared_dims[j] - 1
        dim     = state._shared_dims[j]
        players = sc.players
        ḡ       = zeros(T, dim)

        # Average violation over trajectory and over involved players
        for i in players
            for k in 1:N
                x_k = op.x[k]
                u_k = op.u[i][k]
                ḡ .+= evaluate_constraint(sc.constraint, x_k, u_k, nothing, k)
            end
        end
        ḡ ./= (length(players) * N)

        λ_old = copy(state.λ_shared[rng])

        if state._shared_eq[j]
            state.λ_shared[rng] .+= state.ρ .* ḡ
        else
            state.λ_shared[rng] .= max.(zero(T), state.λ_shared[rng] .+ state.ρ .* ḡ)
        end

        Δλ_max          = max(Δλ_max, maximum(abs.(state.λ_shared[rng] .- λ_old)))
        total_violation = max(total_violation, maximum(abs.(ḡ)))
    end

    return Δλ_max, total_violation
end

# ── Penalty schedule ──────────────────────────────────────────────────────────

"""
    maybe_update_penalty!(state::ALSolverState, objectives, current_violation) -> Bool

Apply the geometric penalty update ρ ← min(ρ_max, φ·ρ) if the constraint
violation has not decreased sufficiently since the last call.

Mirrors FALCON's outer AL penalty schedule exactly:
- Fires if `current_violation > (1 - violation_tol) · prev_violation`
- Caps at `ρ_max` (prevents floating-point cancellation in backward pass)

Updates `ρ` in-place on `state` and synchronizes ρ to all player AL objectives.
Returns `true` if ρ was updated.
"""
function maybe_update_penalty!(
    state::ALSolverState{T},
    objectives::Vector{<:ALAugmentedObjective},
    current_violation::T
) where {T}
    threshold = (1 - state.opts.violation_tol) * state.prev_violation
    updated   = false

    if current_violation > threshold && state.ρ < state.opts.ρ_max
        new_ρ   = min(state.opts.ρ_max, state.opts.φ * state.ρ)
        state.ρ = new_ρ
        # Synchronize ρ to all player AL objectives
        for obj in objectives
            obj.ρ = new_ρ
        end
        updated = true
    end

    state.prev_violation = current_violation
    return updated
end

# ── Constraint violation utilities ───────────────────────────────────────────

"""
    constraint_violation(obj::ALAugmentedObjective, op::OperatingPoint) -> (priv_viol, shared_viol)

Compute ‖g‖_∞ for private and shared constraints over the trajectory.
Used for convergence logging and penalty schedule triggering.

Returns:
- `priv_viol`   : Maximum violation across private constraints and timesteps
- `shared_viol` : Maximum violation across shared constraints and timesteps
"""
function constraint_violation(
    obj::ALAugmentedObjective{T},
    op::OperatingPoint{T}
) where {T}
    N = length(op.x) - 1
    i = obj.base.player_id

    priv_viol   = zero(T)
    shared_viol = zero(T)

    for (j, pc) in enumerate(obj.private_constraints)
        for k in 1:N
            g = evaluate_constraint(pc.constraint, op.x[k], op.u[i][k], nothing, k)
            # Violation for inequality: max(0, g); for equality: |g|
            if obj._priv_eq[j]
                priv_viol = max(priv_viol, maximum(abs.(g)))
            else
                priv_viol = max(priv_viol, maximum(max.(g, zero(T))))
            end
        end
    end

    for (j, sc) in enumerate(obj.shared_constraints)
        for k in 1:N
            g = evaluate_constraint(sc.constraint, op.x[k], op.u[i][k], nothing, k)
            if obj._shared_eq[j]
                shared_viol = max(shared_viol, maximum(abs.(g)))
            else
                shared_viol = max(shared_viol, maximum(max.(g, zero(T))))
            end
        end
    end

    return priv_viol, shared_viol
end

"""
    reset_multipliers!(obj::ALAugmentedObjective{T})

Zero out private multipliers and reset penalty. Used to restart AL from scratch.
Does not touch shared multipliers (those are owned by ALSolverState).
"""
function reset_multipliers!(obj::ALAugmentedObjective{T}, opts::ALOptions{T}) where {T}
    fill!(obj.λ_private, zero(T))
    obj.ρ = opts.ρ_init
end

"""
    reset_al_state!(state::ALSolverState{T})

Reset all shared multipliers and penalty to initial values.
"""
function reset_al_state!(state::ALSolverState{T}) where {T}
    fill!(state.λ_shared, zero(T))
    state.ρ             = state.opts.ρ_init
    state.prev_violation = T(Inf)
end