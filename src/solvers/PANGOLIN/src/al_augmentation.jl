# ============================================================================
# al_augmentation.jl
#
# Phase 3 — Augmented Lagrangian Constraint Layer for iLQGames++
#
# Per-timestep shared multipliers (ALGAMES design):
#   ALSolverState.λ_shared_traj[j] is a (dim_j × N) Matrix{T}.
#   Column k holds the multiplier for shared constraint j at timestep k.
#   This avoids trajectory-averaging dilution when constraints are violated
#   at only a few timesteps (e.g., collision at the crossing point).
# ============================================================================

# ============================================================================
# ALOptions
# ============================================================================

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
        @assert ρ_init > 0
        @assert ρ_max >= ρ_init
        @assert 1 < φ <= 10
        @assert 0 < violation_tol < 1
        @assert dual_tol > 0
        new{T}(ρ_init, ρ_max, φ, violation_tol, dual_tol)
    end
end

ALOptions() = ALOptions(;
    ρ_init=1.0, ρ_max=1e4, φ=10.0, violation_tol=0.25, dual_tol=1e-4
)

# ============================================================================
# Constraint dimension / type helpers
# ============================================================================

constraint_dim(c::LinearConstraint)     = size(c.A, 1)
constraint_dim(c::BoundConstraint)      = 2 * length(c.lower)
constraint_dim(c::NormConstraint)       = 1
constraint_dim(c::NonlinearConstraint)  = c.dim
constraint_dim(pc::PrivateConstraint)   = constraint_dim(pc.constraint)
constraint_dim(sc::SharedConstraint)    = constraint_dim(sc.constraint)

is_equality_constraint(::LinearConstraint)     = false
is_equality_constraint(::BoundConstraint)      = false
is_equality_constraint(::NormConstraint)       = false
is_equality_constraint(c::NonlinearConstraint) = c.constraint_type == :equality
is_equality_constraint(pc::PrivateConstraint)  = is_equality_constraint(pc.constraint)
is_equality_constraint(sc::SharedConstraint)   = is_equality_constraint(sc.constraint)

# ============================================================================
# ALAugmentedObjective
# ============================================================================

"""
    ALAugmentedObjective{T, OBJ}

Wraps a `PlayerObjective` with per-timestep Augmented Lagrangian multipliers.

# Multiplier storage
- Private:  `λ_private[j]`      — `(dim_j × N)` matrix; column k = multiplier at step k
- Shared:   `λ_shared_traj[j]`  — `(dim_j × N)` matrix; centralized, shared by reference
              via `ALSolverState`. All player objectives for the same constraint reference
              the same underlying matrix.

Per-timestep multipliers eliminate the trajectory-averaging dilution that defeats
time-averaged AL when constraints are violated at only a few timesteps.
"""
mutable struct ALAugmentedObjective{T, OBJ <: PlayerObjective}
    base::OBJ
    private_constraints::Vector{<:PrivateConstraint}
    shared_constraints::Vector{<:SharedConstraint}

    λ_private::Vector{Matrix{T}}    # [j] → (dim_j × N)
    λ_shared_traj::Vector{Matrix{T}} # [j] → (dim_j × N); centralized

    ρ::T

    _priv_eq::Vector{Bool}
    _shared_eq::Vector{Bool}
    _priv_dims::Vector{Int}
    _shared_dims::Vector{Int}
end

function ALAugmentedObjective(
    base::OBJ,
    private_constraints::Vector{<:PrivateConstraint},
    shared_constraints::Vector{<:SharedConstraint},
    λ_shared_traj::Vector{Matrix{T}},
    N::Int,
    opts::ALOptions{T}
) where {T, OBJ <: PlayerObjective}

    priv_dims  = [constraint_dim(pc) for pc in private_constraints]
    priv_eq    = [is_equality_constraint(pc) for pc in private_constraints]
    λ_private  = [zeros(T, d, N) for d in priv_dims]

    shared_dims = [constraint_dim(sc) for sc in shared_constraints]
    shared_eq   = [is_equality_constraint(sc) for sc in shared_constraints]

    @assert length(λ_shared_traj) == length(shared_constraints)
    for (j, sc) in enumerate(shared_constraints)
        d = constraint_dim(sc)
        @assert size(λ_shared_traj[j]) == (d, N) "λ_shared_traj[$j] must be ($d × $N), got $(size(λ_shared_traj[j]))"
    end

    ALAugmentedObjective{T, OBJ}(
        base,
        private_constraints, shared_constraints,
        λ_private, λ_shared_traj,
        opts.ρ_init,
        priv_eq, shared_eq,
        priv_dims, shared_dims
    )
end

# ============================================================================
# augmented_stage_cost
# ============================================================================

"""
    augmented_stage_cost(obj, x, u, k; x_joint) -> T

AL-augmented stage cost at (x, u, k) using per-timestep multipliers.
"""
function augmented_stage_cost(
    obj::ALAugmentedObjective{T},
    x::AbstractVector,
    u::AbstractVector,
    k::Int;
    x_joint::Union{Nothing, AbstractVector} = nothing
) where {T}
    cost = T(evaluate_stage_cost(obj.base.stage_cost, x, u, nothing, k)) * obj.base.scaling
    ρ    = obj.ρ

    for (j, pc) in enumerate(obj.private_constraints)
        g = evaluate_constraint(pc.constraint, x, u, nothing, k)
        λ = @view obj.λ_private[j][:, k]
        if obj._priv_eq[j]
            cost += dot(λ, g) + (ρ/2) * dot(g, g)
        else
            shifted = g .+ λ ./ ρ
            active  = max.(shifted, zero(T))
            cost   += (ρ/2) * dot(active, active)
        end
    end

    x_s = isnothing(x_joint) ? x : x_joint
    for (j, sc) in enumerate(obj.shared_constraints)
        g = evaluate_constraint(sc.constraint, x_s, u, nothing, k)
        λ = @view obj.λ_shared_traj[j][:, k]
        if obj._shared_eq[j]
            cost += dot(λ, g) + (ρ/2) * dot(g, g)
        else
            shifted = g .+ λ ./ ρ
            active  = max.(shifted, zero(T))
            cost   += (ρ/2) * dot(active, active)
        end
    end

    return cost
end

# ============================================================================
# AL Quadraticization
# ============================================================================

"""
    quadraticize(obj::ALAugmentedObjective, x_ref, u_ref, k; x_joint, player_state_cols)

Second-order Taylor expansion of the AL cost at (x_ref, u_ref, k) using
per-timestep multipliers `λ[j][:, k]`.
"""
function quadraticize(
    obj::ALAugmentedObjective{T},
    x_ref::AbstractVector,
    u_ref::AbstractVector,
    k::Int;
    x_joint::Union{Nothing, AbstractVector} = nothing,
    player_state_cols::Union{Nothing, UnitRange{Int}} = nothing
) where {T}

    Q, R, M, q, r, c = quadraticize(obj.base.stage_cost, x_ref, u_ref, k)
    Q = Matrix{T}(Q); R = Matrix{T}(R); M = Matrix{T}(M)
    q = Vector{T}(q); r = Vector{T}(r); c = T(c)

    ρ = obj.ρ

    for (j, pc) in enumerate(obj.private_constraints)
        g  = evaluate_constraint(pc.constraint, x_ref, u_ref, nothing, k)
        λj = @view obj.λ_private[j][:, k]
        Jx, Ju = constraint_jacobian(pc.constraint, x_ref, u_ref, nothing, k)
        _add_al_quadratic_update!(Q, R, M, q, r, g, λj, Jx, Ju, ρ, obj._priv_eq[j])
    end

    x_s = isnothing(x_joint) ? x_ref : x_joint
    for (j, sc) in enumerate(obj.shared_constraints)
        g  = evaluate_constraint(sc.constraint, x_s, u_ref, nothing, k)
        λj = @view obj.λ_shared_traj[j][:, k]
        Jx_full, Ju = constraint_jacobian(sc.constraint, x_s, u_ref, nothing, k)
        Jx_local = isnothing(player_state_cols) ? Jx_full :
                   Jx_full[:, player_state_cols]
        _add_al_quadratic_update!(Q, R, M, q, r, g, λj, Jx_local, Ju, ρ, obj._shared_eq[j])
    end

    Q = (Q + Q') / 2
    R = (R + R') / 2
    λ_min_R = minimum(real.(eigvals(Symmetric(R))))
    if λ_min_R < 0
        R = R + (-λ_min_R + sqrt(eps(T))) * I
    end

    return Q, R, M, q, r, c
end

"""
    _add_al_quadratic_update!(Q, R, M, q, r, g, λ, Jx, Ju, ρ, is_eq)

In-place rank-1 Hessian and gradient update for one constraint block.
Row-by-row active-set gating for correct partial active-set handling.
"""
function _add_al_quadratic_update!(
    Q::Matrix{T}, R::Matrix{T}, M::Matrix{T},
    q::Vector{T}, r::Vector{T},
    g::Vector{T}, λ::AbstractVector,
    Jx::AbstractMatrix, Ju::AbstractMatrix,
    ρ::T, is_eq::Bool
) where {T}
    n_rows = length(g)
    @assert length(λ) == n_rows
    @assert size(Jx, 1) == n_rows
    @assert size(Ju, 1) == n_rows

    for d in 1:n_rows
        active = is_eq || (g[d] + λ[d] / ρ > 0)
        active || continue

        μd = λ[d] + ρ * g[d]
        jx = @view Jx[d, :]
        ju = @view Ju[d, :]

        Q .+= ρ .* (jx * jx')
        R .+= ρ .* (ju * ju')
        M .+= ρ .* (jx * ju')
        q .+= μd .* jx
        r .+= μd .* ju
    end
    return nothing
end

# ============================================================================
# ALSolverState — per-timestep shared multipliers
# ============================================================================

"""
    ALSolverState{T}

Solver-level AL state owning centralized per-timestep shared multipliers.

# Fields
- `λ_shared_traj`  : `Vector{Matrix{T}}`; `λ_shared_traj[j]` is `(dim_j × N)`.
                     All player `ALAugmentedObjective` structs reference these matrices.
- `ρ`              : Current penalty weight (synchronized across all objectives)
- `prev_violation` : Max pointwise active violation from previous outer iteration
- `opts`           : ALOptions
- `_shared_eq`     : Equality flags per shared constraint
- `_shared_dims`   : Output dimension per shared constraint
"""
mutable struct ALSolverState{T}
    λ_shared_traj::Vector{Matrix{T}}
    ρ::T
    prev_violation::T
    opts::ALOptions{T}
    _shared_eq::Vector{Bool}
    _shared_dims::Vector{Int}
end

"""
    ALSolverState(shared_constraints, opts, N) -> ALSolverState{T}

Construct solver-level AL state for N timesteps.
"""
function ALSolverState(
    shared_constraints::Vector{<:SharedConstraint},
    opts::ALOptions{T},
    N::Int
) where {T}
    shared_dims   = [constraint_dim(sc) for sc in shared_constraints]
    shared_eq     = [is_equality_constraint(sc) for sc in shared_constraints]
    λ_shared_traj = [zeros(T, d, N) for d in shared_dims]
    ALSolverState{T}(λ_shared_traj, opts.ρ_init, T(Inf), opts, shared_eq, shared_dims)
end

# ============================================================================
# Dual updates
# ============================================================================

"""
    update_multipliers!(obj, op, game) -> Δλ_max

Per-timestep private dual ascent.
`λ_private[j][:, k] ← max(0, λ + ρ·g(x_k, u_k))` for each k.
"""
function update_multipliers!(
    obj::ALAugmentedObjective{T},
    op::OperatingPoint{T},
    game::Union{Nothing, GameProblem} = nothing
) where {T}
    N      = length(op.x) - 1
    ρ      = obj.ρ
    i      = obj.base.player_id
    Δλ_max = zero(T)

    is_sep       = !isnothing(game) && game.dynamics isa SeparableDynamics
    state_offset = is_sep ? game.metadata.state_offsets[i] : 0
    state_dim_i  = is_sep ? game.metadata.state_dims[i]    : 0

    for (j, pc) in enumerate(obj.private_constraints)
        for k in 1:N
            x_k = is_sep ? op.x[k][state_offset+1:state_offset+state_dim_i] : op.x[k]
            u_k = op.u[i][k]
            g_k = evaluate_constraint(pc.constraint, x_k, u_k, nothing, k)

            λ_old = copy(obj.λ_private[j][:, k])
            if obj._priv_eq[j]
                obj.λ_private[j][:, k] .+= ρ .* g_k
            else
                obj.λ_private[j][:, k] .= max.(zero(T), obj.λ_private[j][:, k] .+ ρ .* g_k)
            end
            Δλ_max = max(Δλ_max, maximum(abs.(obj.λ_private[j][:, k] .- λ_old)))
        end
    end

    return Δλ_max
end

"""
    update_shared_multipliers!(state, shared_constraints, objectives, op) -> (Δλ_max, total_violation)

Per-timestep shared dual ascent.

For each constraint j and timestep k:
    `λ_shared_traj[j][:, k] ← max(0, λ + ρ · g(x_k, u_k))`

Violation = `max_{j,k} max(0, g_j(x_k))` — pointwise, no trajectory averaging.
This ensures instantaneous constraint violations are not masked by satisfied
timesteps earlier or later in the trajectory.
"""
function update_shared_multipliers!(
    state::ALSolverState{T},
    shared_constraints::Vector{<:SharedConstraint},
    objectives::Vector{<:ALAugmentedObjective},
    op::OperatingPoint{T}
) where {T}
    N = length(op.x) - 1

    Δλ_max          = zero(T)
    total_violation = zero(T)

    for (j, sc) in enumerate(shared_constraints)
        players = sc.players
        is_eq   = state._shared_eq[j]
        Λ       = state.λ_shared_traj[j]   # (dim_j × N); mutated in-place

        for k in 1:N
            x_k = op.x[k]
            # Use first player's control as representative (collision/hallway
            # constraints don't depend on u; control-coupling constraints should
            # use a proper aggregation, but that's problem-specific)
            u_k = op.u[players[1]][k]

            g_k = evaluate_constraint(sc.constraint, x_k, u_k, nothing, k)

            λ_old = copy(Λ[:, k])
            if is_eq
                # Equality: unclamped dual ascent
                Λ[:, k] .+= state.ρ .* g_k
                active_viol = maximum(abs.(g_k))
            else
                # Inequality: standard AL dual ascent (ALGAMES Eq. 10)
                #   λ ← max(0, λ + ρ·g)
                # This is the correct update from AL theory. When g < 0 (constraint
                # satisfied), λ decreases toward zero — this is correct complementary
                # slackness behaviour, not a bug.
                Λ[:, k] .= max.(zero(T), Λ[:, k] .+ state.ρ .* g_k)
                active_viol = maximum(max.(g_k, zero(T)))
            end

            Δλ_max          = max(Δλ_max, maximum(abs.(Λ[:, k] .- λ_old)))
            total_violation = max(total_violation, active_viol)
        end
    end

    return Δλ_max, total_violation
end

# ============================================================================
# Penalty schedule
# ============================================================================

"""
    update_penalty!(state, objectives) -> Bool

Unconditional geometric penalty update ρ ← min(ρ_max, φ·ρ), applied every
outer iteration per ALGAMES Algorithm 3 / Equation 11. Returns true if ρ
was increased (i.e., not yet capped at ρ_max).

Unlike a violation-triggered schedule, the unconditional increase ensures
ρ → ρ_max in finite outer iterations, which is required for AL convergence:
the dual update λ ← max(0, λ + ρg) only approximates the true multiplier
to within O(1/ρ), so ρ must grow to drive this error to zero.
"""
function update_penalty!(
    state::ALSolverState{T},
    objectives::Vector{<:ALAugmentedObjective}
) where {T}
    updated = false
    if state.ρ < state.opts.ρ_max
        state.ρ = min(state.opts.ρ_max, state.opts.φ * state.ρ)
        for obj in objectives
            obj.ρ = state.ρ
        end
        updated = true
    end
    return updated
end

# Keep maybe_update_penalty! as an alias for backward compatibility with tests
function maybe_update_penalty!(
    state::ALSolverState{T},
    objectives::Vector{<:ALAugmentedObjective},
    current_violation::T   # retained in signature for API compatibility; unused
) where {T}
    return update_penalty!(state, objectives)
end

# ============================================================================
# Utilities
# ============================================================================

"""
    constraint_violation(obj, op) -> (priv_viol, shared_viol)

Pointwise max active violation over trajectory (no averaging).
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
            v = obj._priv_eq[j] ? maximum(abs.(g)) : maximum(max.(g, zero(T)))
            priv_viol = max(priv_viol, v)
        end
    end

    for (j, sc) in enumerate(obj.shared_constraints)
        for k in 1:N
            g = evaluate_constraint(sc.constraint, op.x[k], op.u[i][k], nothing, k)
            v = obj._shared_eq[j] ? maximum(abs.(g)) : maximum(max.(g, zero(T)))
            shared_viol = max(shared_viol, v)
        end
    end

    return priv_viol, shared_viol
end

"""
    reset_multipliers!(obj, opts)

Zero private per-timestep multipliers and reset ρ.
"""
function reset_multipliers!(obj::ALAugmentedObjective{T}, opts::ALOptions{T}) where {T}
    for λ in obj.λ_private; fill!(λ, zero(T)); end
    obj.ρ = opts.ρ_init
end

"""
    reset_al_state!(state)

Zero all shared per-timestep multipliers and reset ρ.
"""
function reset_al_state!(state::ALSolverState{T}) where {T}
    for Λ in state.λ_shared_traj; fill!(Λ, zero(T)); end
    state.ρ              = state.opts.ρ_init
    state.prev_violation = T(Inf)
end