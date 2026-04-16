# # ============================================================================
# # al_augmentation.jl
# #
# # Phase 3 — Augmented Lagrangian Constraint Layer for iLQGames++
# #
# # Per-timestep shared multipliers (ALGAMES design):
# #   ALSolverState.λ_shared_traj[j] is a (dim_j × N) Matrix{T}.
# #   Column k holds the multiplier for shared constraint j at timestep k.
# #   This avoids trajectory-averaging dilution when constraints are violated
# #   at only a few timesteps (e.g., collision at the crossing point).
# # ============================================================================

# # ============================================================================
# # ALOptions
# # ============================================================================

# struct ALOptions{T}
#     ρ_init::T
#     ρ_max::T
#     φ::T
#     violation_tol::T
#     dual_tol::T
#     φ_threshold::T              # violation must improve by this factor to skip ρ increase
#     max_consecutive_increases::Int  # hard cap: pause ρ growth for this many iters after increase

#     function ALOptions(;
#         ρ_init::T        = 1.0,
#         ρ_max::T         = 1e4,
#         φ::T             = 10.0,
#         violation_tol::T = 1e-3,
#         dual_tol::T      = 1e-4,
#         φ_threshold::T   = 0.75,
#         max_consecutive_increases::Int = 3
#     ) where {T <: AbstractFloat}
#         @assert ρ_init > 0
#         @assert ρ_max >= ρ_init
#         @assert 1 < φ <= 20
#         @assert 0 < violation_tol < 1
#         @assert dual_tol > 0
#         @assert 0 < φ_threshold < 1
#         @assert max_consecutive_increases >= 1
#         new{T}(ρ_init, ρ_max, φ, violation_tol, dual_tol,
#                φ_threshold, max_consecutive_increases)
#     end
# end

# ALOptions() = ALOptions(;
#     ρ_init=1.0, ρ_max=1e4, φ=10.0, violation_tol=1e-3, dual_tol=1e-4,
#     φ_threshold=0.75, max_consecutive_increases=3
# )

# # ============================================================================
# # Constraint dimension / type helpers
# # ============================================================================

# constraint_dim(c::LinearConstraint)     = size(c.A, 1)
# constraint_dim(c::BoundConstraint)      = 2 * length(c.lower)
# constraint_dim(c::NormConstraint)       = 1
# constraint_dim(c::NonlinearConstraint)  = c.dim
# constraint_dim(pc::PrivateConstraint)   = constraint_dim(pc.constraint)
# constraint_dim(sc::SharedConstraint)    = constraint_dim(sc.constraint)

# is_equality_constraint(::LinearConstraint)     = false
# is_equality_constraint(::BoundConstraint)      = false
# is_equality_constraint(::NormConstraint)       = false
# is_equality_constraint(c::NonlinearConstraint) = c.constraint_type == :equality
# is_equality_constraint(pc::PrivateConstraint)  = is_equality_constraint(pc.constraint)
# is_equality_constraint(sc::SharedConstraint)   = is_equality_constraint(sc.constraint)

# # ============================================================================
# # ALAugmentedObjective
# # ============================================================================

# """
#     ALAugmentedObjective{T, OBJ}

# Wraps a `PlayerObjective` with per-timestep Augmented Lagrangian multipliers.

# # Multiplier storage
# - Private:  `λ_private[j]`      — `(dim_j × N)` matrix; column k = multiplier at step k
# - Shared:   `λ_shared_traj[j]`  — `(dim_j × N)` matrix; centralized, shared by reference
#               via `ALSolverState`. All player objectives for the same constraint reference
#               the same underlying matrix.

# Per-timestep multipliers eliminate the trajectory-averaging dilution that defeats
# time-averaged AL when constraints are violated at only a few timesteps.
# """
# mutable struct ALAugmentedObjective{T, OBJ <: PlayerObjective}
#     base::OBJ
#     private_constraints::Vector{<:PrivateConstraint}
#     shared_constraints::Vector{<:SharedConstraint}

#     λ_private::Vector{Matrix{T}}    # [j] → (dim_j × N)
#     λ_shared_traj::Vector{Matrix{T}} # [j] → (dim_j × N); centralized

#     ρ::T

#     _priv_eq::Vector{Bool}
#     _shared_eq::Vector{Bool}
#     _priv_dims::Vector{Int}
#     _shared_dims::Vector{Int}
# end

# function ALAugmentedObjective(
#     base::OBJ,
#     private_constraints::Vector{<:PrivateConstraint},
#     shared_constraints::Vector{<:SharedConstraint},
#     λ_shared_traj::Vector{Matrix{T}},
#     N::Int,
#     opts::ALOptions{T}
# ) where {T, OBJ <: PlayerObjective}

#     priv_dims  = [constraint_dim(pc) for pc in private_constraints]
#     priv_eq    = [is_equality_constraint(pc) for pc in private_constraints]
#     # N+1 columns: k=1:N for stage, k=N+1 for terminal state
#     λ_private  = [zeros(T, d, N+1) for d in priv_dims]

#     shared_dims = [constraint_dim(sc) for sc in shared_constraints]
#     shared_eq   = [is_equality_constraint(sc) for sc in shared_constraints]

#     @assert length(λ_shared_traj) == length(shared_constraints)
#     for (j, sc) in enumerate(shared_constraints)
#         d = constraint_dim(sc)
#         @assert size(λ_shared_traj[j]) == (d, N) "λ_shared_traj[$j] must be ($d × $N), got $(size(λ_shared_traj[j]))"
#     end

#     ALAugmentedObjective{T, OBJ}(
#         base,
#         private_constraints, shared_constraints,
#         λ_private, λ_shared_traj,
#         opts.ρ_init,
#         priv_eq, shared_eq,
#         priv_dims, shared_dims
#     )
# end

# # ============================================================================
# # augmented_stage_cost
# # ============================================================================

# """
#     augmented_stage_cost(obj, x, u, k; x_joint) -> T

# AL-augmented stage cost at (x, u, k) using per-timestep multipliers.
# """
# function augmented_stage_cost(
#     obj::ALAugmentedObjective{T},
#     x::AbstractVector,
#     u::AbstractVector,
#     k::Int;
#     x_joint::Union{Nothing, AbstractVector} = nothing
# ) where {T}
#     cost = T(evaluate_stage_cost(obj.base.stage_cost, x, u, nothing, k)) * obj.base.scaling
#     ρ    = obj.ρ

#     for (j, pc) in enumerate(obj.private_constraints)
#         g = evaluate_constraint(pc.constraint, x, u, nothing, k)
#         λ = @view obj.λ_private[j][:, k]
#         if obj._priv_eq[j]
#             cost += dot(λ, g) + (ρ/2) * dot(g, g)
#         else
#             shifted = g .+ λ ./ ρ
#             active  = max.(shifted, zero(T))
#             cost   += (ρ/2) * dot(active, active)
#         end
#     end

#     x_s = isnothing(x_joint) ? x : x_joint
#     for (j, sc) in enumerate(obj.shared_constraints)
#         g = evaluate_constraint(sc.constraint, x_s, u, nothing, k)
#         λ = @view obj.λ_shared_traj[j][:, k]
#         if obj._shared_eq[j]
#             cost += dot(λ, g) + (ρ/2) * dot(g, g)
#         else
#             shifted = g .+ λ ./ ρ
#             active  = max.(shifted, zero(T))
#             cost   += (ρ/2) * dot(active, active)
#         end
#     end

#     return cost
# end

# # ============================================================================
# # AL Quadraticization
# # ============================================================================

# """
#     quadraticize(obj::ALAugmentedObjective, x_ref, u_ref, k; x_joint, player_state_cols)

# Second-order Taylor expansion of the AL cost at (x_ref, u_ref, k) using
# per-timestep multipliers `λ[j][:, k]`.
# """
# function quadraticize(
#     obj::ALAugmentedObjective{T},
#     x_ref::AbstractVector,
#     u_ref::AbstractVector,
#     k::Int;
#     x_joint::Union{Nothing, AbstractVector} = nothing,
#     player_state_cols::Union{Nothing, UnitRange{Int}} = nothing
# ) where {T}

#     Q, R, M, q, r, c = quadraticize(obj.base.stage_cost, x_ref, u_ref, k)
#     Q = Matrix{T}(Q); R = Matrix{T}(R); M = Matrix{T}(M)
#     q = Vector{T}(q); r = Vector{T}(r); c = T(c)

#     ρ = obj.ρ

#     for (j, pc) in enumerate(obj.private_constraints)
#         g  = evaluate_constraint(pc.constraint, x_ref, u_ref, nothing, k)
#         λj = @view obj.λ_private[j][:, k]
#         Jx, Ju = constraint_jacobian(pc.constraint, x_ref, u_ref, nothing, k)
#         _add_al_quadratic_update!(Q, R, M, q, r, g, λj, Jx, Ju, ρ, obj._priv_eq[j])
#     end

#     x_s = isnothing(x_joint) ? x_ref : x_joint
#     for (j, sc) in enumerate(obj.shared_constraints)
#         g  = evaluate_constraint(sc.constraint, x_s, u_ref, nothing, k)
#         λj = @view obj.λ_shared_traj[j][:, k]
#         Jx_full, Ju = constraint_jacobian(sc.constraint, x_s, u_ref, nothing, k)
#         Jx_local = isnothing(player_state_cols) ? Jx_full :
#                    Jx_full[:, player_state_cols]
#         _add_al_quadratic_update!(Q, R, M, q, r, g, λj, Jx_local, Ju, ρ, obj._shared_eq[j])
#     end

#     Q = (Q + Q') / 2
#     R = (R + R') / 2
#     λ_min_R = minimum(real.(eigvals(Symmetric(R))))
#     if λ_min_R < 0
#         R = R + (-λ_min_R + sqrt(eps(T))) * I
#     end

#     return Q, R, M, q, r, c
# end

# """
#     _add_al_quadratic_update!(Q, R, M, q, r, g, λ, Jx, Ju, ρ, is_eq)

# In-place rank-1 Hessian and gradient update for one constraint block.
# Row-by-row active-set gating for correct partial active-set handling.
# """
# function _add_al_quadratic_update!(
#     Q::Matrix{T}, R::Matrix{T}, M::Matrix{T},
#     q::Vector{T}, r::Vector{T},
#     g::Vector{T}, λ::AbstractVector,
#     Jx::AbstractMatrix, Ju::AbstractMatrix,
#     ρ::T, is_eq::Bool
# ) where {T}
#     n_rows = length(g)
#     @assert length(λ) == n_rows
#     @assert size(Jx, 1) == n_rows
#     @assert size(Ju, 1) == n_rows

#     for d in 1:n_rows
#         active = is_eq || (g[d] + λ[d] / ρ > 0)
#         active || continue

#         μd = λ[d] + ρ * g[d]
#         jx = @view Jx[d, :]
#         ju = @view Ju[d, :]

#         Q .+= ρ .* (jx * jx')
#         R .+= ρ .* (ju * ju')
#         M .+= ρ .* (jx * ju')
#         q .+= μd .* jx
#         r .+= μd .* ju
#     end
#     return nothing
# end

# # ============================================================================
# # ALSolverState — per-timestep shared multipliers
# # ============================================================================

# """
#     ALSolverState{T}

# Solver-level AL state owning centralized per-timestep shared multipliers.

# # Fields
# - `λ_shared_traj`         : `Vector{Matrix{T}}`; `λ_shared_traj[j]` is `(dim_j × N)`.
# - `ρ`                     : Current penalty weight (synchronized across all objectives)
# - `prev_violation`        : Max pointwise active violation from previous outer iteration
# - `opts`                  : ALOptions
# - `_shared_eq`            : Equality flags per shared constraint
# - `_shared_dims`          : Output dimension per shared constraint
# - `consecutive_increases` : Count of consecutive ρ increases; reset to 0 when held.
#                             When this reaches `opts.max_consecutive_increases`, ρ is
#                             held for one iteration to let the LQ approximation stabilize
#                             before resuming increases.
# """
# mutable struct ALSolverState{T}
#     λ_shared_traj::Vector{Matrix{T}}
#     ρ::T
#     prev_violation::T
#     opts::ALOptions{T}
#     _shared_eq::Vector{Bool}
#     _shared_dims::Vector{Int}
#     consecutive_increases::Int
# end

# """
#     ALSolverState(shared_constraints, opts, N) -> ALSolverState{T}

# Construct solver-level AL state for N timesteps.
# """
# function ALSolverState(
#     shared_constraints::Vector{<:SharedConstraint},
#     opts::ALOptions{T},
#     N::Int
# ) where {T}
#     shared_dims   = [constraint_dim(sc) for sc in shared_constraints]
#     shared_eq     = [is_equality_constraint(sc) for sc in shared_constraints]
#     λ_shared_traj = [zeros(T, d, N) for d in shared_dims]
#     ALSolverState{T}(λ_shared_traj, opts.ρ_init, T(Inf), opts, shared_eq, shared_dims, 0)
# end

# # ============================================================================
# # Dual updates
# # ============================================================================

# """
#     update_multipliers!(obj, op, game) -> Δλ_max

# Per-timestep private dual ascent.
# `λ_private[j][:, k] ← max(0, λ + ρ·g(x_k, u_k))` for each k.
# """
# function update_multipliers!(
#     obj::ALAugmentedObjective{T},
#     op::OperatingPoint{T},
#     game::Union{Nothing, GameProblem} = nothing
# ) where {T}
#     N      = length(op.x) - 1
#     ρ      = obj.ρ
#     i      = obj.base.player_id
#     Δλ_max = zero(T)

#     is_sep       = !isnothing(game) && game.dynamics isa SeparableDynamics
#     state_offset = is_sep ? game.metadata.state_offsets[i] : 0
#     state_dim_i  = is_sep ? game.metadata.state_dims[i]    : 0

#     for (j, pc) in enumerate(obj.private_constraints)
#         for k in 1:N+1
#             x_k = is_sep ? op.x[k][state_offset+1:state_offset+state_dim_i] : op.x[k]
#             u_k = k <= N ? op.u[i][k] : zeros(T, length(op.u[i][1]))
#             g_k = evaluate_constraint(pc.constraint, x_k, u_k, nothing, k)

#             λ_old = copy(obj.λ_private[j][:, k])
#             if obj._priv_eq[j]
#                 obj.λ_private[j][:, k] .+= ρ .* g_k
#             else
#                 obj.λ_private[j][:, k] .= max.(zero(T), obj.λ_private[j][:, k] .+ ρ .* g_k)
#             end
#             Δλ_max = max(Δλ_max, maximum(abs.(obj.λ_private[j][:, k] .- λ_old)))
#         end
#     end

#     return Δλ_max
# end

# """
#     update_shared_multipliers!(state, shared_constraints, objectives, op) -> (Δλ_max, total_violation)

# Per-timestep shared dual ascent.

# For each constraint j and timestep k:
#     `λ_shared_traj[j][:, k] ← max(0, λ + ρ · g(x_k, u_k))`

# Violation = `max_{j,k} max(0, g_j(x_k))` — pointwise, no trajectory averaging.
# This ensures instantaneous constraint violations are not masked by satisfied
# timesteps earlier or later in the trajectory.
# """
# function update_shared_multipliers!(
#     state::ALSolverState{T},
#     shared_constraints::Vector{<:SharedConstraint},
#     objectives::Vector{<:ALAugmentedObjective},
#     op::OperatingPoint{T}
# ) where {T}
#     N = length(op.x) - 1

#     Δλ_max          = zero(T)
#     total_violation = zero(T)

#     for (j, sc) in enumerate(shared_constraints)
#         players = sc.players
#         is_eq   = state._shared_eq[j]
#         Λ       = state.λ_shared_traj[j]   # (dim_j × N); mutated in-place

#         for k in 1:N
#             x_k = op.x[k]
#             u_k = op.u[players[1]][k]

#             g_k = evaluate_constraint(sc.constraint, x_k, u_k, nothing, k)

#             λ_old = copy(Λ[:, k])
#             if is_eq
#                 Λ[:, k] .+= state.ρ .* g_k
#                 active_viol = maximum(abs.(g_k))
#             else
#                 # Standard AL dual ascent (ALGAMES Eq. 10):
#                 #   λ ← max(0, λ + ρ·g)
#                 # When g < 0 (satisfied), λ decreases toward zero — correct
#                 # complementary slackness behaviour.
#                 Λ[:, k] .= max.(zero(T), Λ[:, k] .+ state.ρ .* g_k)
#                 active_viol = maximum(max.(g_k, zero(T)))
#             end

#             Δλ_max          = max(Δλ_max, maximum(abs.(Λ[:, k] .- λ_old)))
#             total_violation = max(total_violation, active_viol)
#         end
#     end

#     return Δλ_max, total_violation
# end

# # ============================================================================
# # Penalty schedule
# # ============================================================================

# """
#     update_penalty!(state, objectives, current_violation) -> Bool

# Violation-gated geometric penalty update with consecutive-increase cap.

# Two conditions must both be true to increase ρ:
#   1. Violation has **not** improved by `φ_threshold` fraction vs. previous iter.
#   2. The number of consecutive ρ increases has not reached `max_consecutive_increases`.

# Condition 1 (violation-gating): only escalate ρ when the dual update alone is
# stalling. Lets multipliers do cheap work at fixed ρ before escalating.

# Condition 2 (consecutive cap): after `max_consecutive_increases` consecutive
# increases, hold ρ for one iteration. This gives the LQ approximation time to
# re-linearize around the new, higher-ρ cost landscape before ρ grows again.
# Without this cap, rapid consecutive increases (e.g., violation oscillates near
# the threshold) produce ill-conditioned LQ subproblems whose solutions involve
# impulsive controls — large α[k] at a single timestep to avoid the constraint —
# which the deviation-based line search cannot detect because the deviation is
# measured from the previous (also-impulsive) operating point.

# Returns `true` if ρ was increased.
# """
# function update_penalty!(
#     state::ALSolverState{T},
#     objectives::Vector{<:ALAugmentedObjective},
#     current_violation::T
# ) where {T}
#     improved = current_violation < state.opts.φ_threshold * state.prev_violation
#     state.prev_violation = current_violation

#     # Determine whether to increase ρ this iteration
#     want_increase   = !improved && state.ρ < state.opts.ρ_max
#     budget_remains  = state.consecutive_increases < state.opts.max_consecutive_increases

#     updated = false
#     if want_increase && budget_remains
#         state.ρ = min(state.opts.ρ_max, state.opts.φ * state.ρ)
#         for obj in objectives
#             obj.ρ = state.ρ
#         end
#         state.consecutive_increases += 1
#         updated = true
#     else
#         # Either improved (good), hit the consecutive cap (pause), or already at ρ_max.
#         # Reset the counter so the cap enforces a minimum 1-iter hold after each burst.
#         state.consecutive_increases = 0
#     end
#     return updated
# end

# """
#     maybe_update_penalty!(state, objectives, current_violation) -> Bool

# Alias for `update_penalty!` — retained for backward compatibility.
# """
# function maybe_update_penalty!(
#     state::ALSolverState{T},
#     objectives::Vector{<:ALAugmentedObjective},
#     current_violation::T
# ) where {T}
#     return update_penalty!(state, objectives, current_violation)
# end

# # ============================================================================
# # Utilities
# # ============================================================================

# """
#     constraint_violation(obj, op) -> (priv_viol, shared_viol)

# Pointwise max active violation over trajectory (no averaging).
# """
# function constraint_violation(
#     obj::ALAugmentedObjective{T},
#     op::OperatingPoint{T}
# ) where {T}
#     N = length(op.x) - 1
#     i = obj.base.player_id

#     priv_viol   = zero(T)
#     shared_viol = zero(T)

#     for (j, pc) in enumerate(obj.private_constraints)
#         for k in 1:N
#             g = evaluate_constraint(pc.constraint, op.x[k], op.u[i][k], nothing, k)
#             v = obj._priv_eq[j] ? maximum(abs.(g)) : maximum(max.(g, zero(T)))
#             priv_viol = max(priv_viol, v)
#         end
#     end

#     for (j, sc) in enumerate(obj.shared_constraints)
#         for k in 1:N
#             g = evaluate_constraint(sc.constraint, op.x[k], op.u[i][k], nothing, k)
#             v = obj._shared_eq[j] ? maximum(abs.(g)) : maximum(max.(g, zero(T)))
#             shared_viol = max(shared_viol, v)
#         end
#     end

#     return priv_viol, shared_viol
# end

# """
#     reset_multipliers!(obj, opts)

# Zero private per-timestep multipliers and reset ρ.
# """
# function reset_multipliers!(obj::ALAugmentedObjective{T}, opts::ALOptions{T}) where {T}
#     for λ in obj.λ_private; fill!(λ, zero(T)); end
#     obj.ρ = opts.ρ_init
# end

# """
#     reset_al_state!(state)

# Zero all shared per-timestep multipliers and reset ρ.
# """
# function reset_al_state!(state::ALSolverState{T}) where {T}
#     for Λ in state.λ_shared_traj; fill!(Λ, zero(T)); end
#     state.ρ                     = state.opts.ρ_init
#     state.prev_violation        = T(Inf)
#     state.consecutive_increases = 0
# end