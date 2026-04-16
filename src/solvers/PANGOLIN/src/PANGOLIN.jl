# # ============================================================================
# # PANGOLIN.jl
# #
# # Penalty-Augmented Nash Game via Operating-point Linearization and
# # Iterative Nash — iLQGames++ solver with AL constraint handling.
# #
# # Solver struct follows the DifferentialGamesBase GameSolver protocol:
# #   struct PANGOLIN <: GameSolver  (immutable parameters)
# #   PANGOLINCache{T}               (mutable pre-allocated scratch, per-game)
# #   _solve(game, solver, warmstart, verbose) -> GameSolution
# #
# # Dependencies:
# #   quadraticization.jl   → OperatingPoint, LTVLQApproximation, lq_approximation!
# #   al_augmentation.jl    → ALOptions, ALAugmentedObjective, ALSolverState
# #   FeedbackStrategy, solve_lq_game!, rollout!, backtrack_scale!, has_converged
# #   (defined in this file — they are PANGOLIN-internal types/functions)
# # ============================================================================

# # ============================================================================
# # FeedbackStrategy
# # ============================================================================

# """
#     FeedbackStrategy{T}

# Affine feedback strategy for N players over a finite horizon.

#     uᵢ(k) = -Pᵢ(k) x(k) - αᵢ(k)

# # Fields
# - `P` : `P[i][k]` is (mᵢ × n) feedback gain for player i at step k
# - `α` : `α[i][k]` is (mᵢ,) feed-forward offset for player i at step k
# """
# struct FeedbackStrategy{T}
#     P::Vector{Vector{Matrix{T}}}    # [n_players][N]  (mᵢ × n)
#     α::Vector{Vector{Vector{T}}}    # [n_players][N]  (mᵢ,)
# end

# function FeedbackStrategy(game::GameProblem{T}) where {T}
#     N         = n_steps(game)
#     n         = total_state_dim(game.dynamics)
#     ctrl_dims = game.metadata.control_dims
#     n_p       = game.n_players
#     P = [[zeros(T, ctrl_dims[i], n) for _ in 1:N] for i in 1:n_p]
#     α = [[zeros(T, ctrl_dims[i])    for _ in 1:N] for i in 1:n_p]
#     FeedbackStrategy{T}(P, α)
# end

# """
#     scale_feedforward!(strat::FeedbackStrategy, sf)

# Multiply all feed-forward terms α[i][k] by scalar `sf` in-place.
# Feedback gains P are unchanged — this is the iLQGames line-search operation.
# """
# function scale_feedforward!(strat::FeedbackStrategy{T}, sf::Real) where {T}
#     for αi in strat.α, αk in αi
#         αk .*= T(sf)
#     end
#     return strat
# end

# function Base.copyto!(dst::FeedbackStrategy{T}, src::FeedbackStrategy{T}) where {T}
#     for i in eachindex(dst.P), k in eachindex(dst.P[i])
#         copyto!(dst.P[i][k], src.P[i][k])
#         copyto!(dst.α[i][k], src.α[i][k])
#     end
#     return dst
# end

# # ============================================================================
# # PANGOLINCache — per-game pre-allocated scratch
# # ============================================================================

# """
#     PANGOLINCache{T}

# Pre-allocated mutable scratch storage for a single PANGOLIN solve.
# Constructed once from a `GameProblem` and reused across iterations.
# """
# mutable struct PANGOLINCache{T}
#     lq_mem::LTVLQApproximation{T}
#     op_mem::OperatingPoint{T}
#     strat_mem::FeedbackStrategy{T}
# end

# function PANGOLINCache(game::GameProblem{T}) where {T}
#     PANGOLINCache{T}(
#         LTVLQApproximation(game),
#         OperatingPoint(game),
#         FeedbackStrategy(game)
#     )
# end

# # ============================================================================
# # PANGOLIN solver struct
# # ============================================================================

# """
#     PANGOLIN <: GameSolver

# **P**enalty-**A**ugmented **N**ash **G**ame via **O**perating-point
# **L**inearization and **I**terative **N**ash.

# Computes feedback Nash equilibria for general nonlinear N-player differential
# games with inequality/equality constraints via iterative LQ approximation
# (iLQGames) and an Augmented Lagrangian outer loop.

# # Algorithm (per outer iteration)
# 1. `lq_approximation!`   — linearize dynamics, quadraticize AL-penalized costs
# 2. `solve_lq_game!`      — backward DP on LTV LQ approximation → FeedbackStrategy
# 3. `backtrack_scale!`    — trajectory-deviation line search on feed-forward term α
# 4. Dual update           — `update_multipliers!` + `update_shared_multipliers!`
# 5. Penalty schedule      — `update_penalty!` (violation-gated; see al_augmentation.jl)
# 6. Convergence check     — primal + feasibility + cost stability

# # Fields
# - `max_iter`              : Maximum outer iterations (default 200)
# - `max_line_search_iter`  : Maximum backtracking steps per iteration (default 30)
# - `α_init`                : Initial feed-forward scale factor (default 0.5)
# - `α_step`                : Geometric reduction per backtrack step (default 0.5)
# - `convergence_tol`       : ‖x_new - x_old‖_∞ threshold (default 1e-3)
# - `state_regularization`  : Tikhonov regularization added to Q (default 0.0)
# - `control_regularization`: Tikhonov regularization added to R (default 0.0)
# - `constraint_opts`       : `ALOptions` controlling penalty schedule
# - `init_mode`             : Initialization ∈ {:zero_control, :straight_line, :warmstart}
# - `use_strategy_warmstart`: Initialize feedback gains from warmstart (default false)
# - `dyn_consistency_tol`   : Warmstart dynamic consistency tolerance (default 1e-3)
# - `check_singularity`     : Warn on ill-conditioned S matrix (default true)
# - `rcond_threshold`       : rcond threshold for singularity warning (default 1e-10)
# """
# struct PANGOLIN <: GameSolver
#     max_iter::Int
#     max_line_search_iter::Int
#     α_init::Float64
#     α_step::Float64
#     armijo_param::Float64
#     convergence_tol::Float64
#     max_elwise_diff_step::Float64
#     state_regularization::Float64
#     control_regularization::Float64
#     constraint_opts::ALOptions{Float64}
#     init_mode::Symbol
#     use_strategy_warmstart::Bool
#     dyn_consistency_tol::Float64
#     check_singularity::Bool
#     rcond_threshold::Float64
#     control_limit::Float64   # reject rollouts whose max |u| exceeds this

#     function PANGOLIN(;
#         max_iter::Int               = 200,
#         max_line_search_iter::Int   = 30,
#         α_init::Float64             = 0.5,
#         α_step::Float64             = 0.5,
#         armijo_param::Float64       = 1e-4,
#         convergence_tol::Float64    = 1e-3,
#         max_elwise_diff_step::Float64 = 0.0,
#         state_regularization::Float64  = 0.0,
#         control_regularization::Float64 = 0.0,
#         constraint_opts::ALOptions{Float64} = ALOptions(),
#         init_mode::Symbol           = :straight_line,
#         use_strategy_warmstart::Bool = false,
#         dyn_consistency_tol::Float64 = 1e-3,
#         check_singularity::Bool     = true,
#         rcond_threshold::Float64    = 1e-10,
#         control_limit::Float64      = Inf
#     )
#         @assert max_iter > 0
#         @assert max_line_search_iter > 0
#         @assert 0 < α_init <= 1
#         @assert 0 < α_step <= 1
#         @assert 0 < armijo_param < 0.5  "armijo_param β must be in (0, 0.5)"
#         @assert convergence_tol > 0
#         @assert init_mode in (:zero_control, :straight_line, :warmstart)
#         @assert control_limit > 0

#         step_limit = max_elwise_diff_step > 0 ? max_elwise_diff_step : 30.0 * convergence_tol

#         new(
#             max_iter, max_line_search_iter,
#             α_init, α_step, armijo_param,
#             convergence_tol, step_limit,
#             state_regularization, control_regularization,
#             constraint_opts, init_mode,
#             use_strategy_warmstart, dyn_consistency_tol,
#             check_singularity, rcond_threshold,
#             control_limit
#         )
#     end
# end

# solver_capabilities(::Type{PANGOLIN}) = [
#     :NonlinearGame,
#     :FeedbackPolicies,
#     :ConstrainedGame,
#     :DiscreteTime
# ]

# # ============================================================================
# # solve_lq_game! — backward DP on LTVLQApproximation
# # ============================================================================

# """
#     solve_lq_game!(strat::FeedbackStrategy{T}, lqg::LTVLQApproximation{T}, solver::PANGOLIN)

# Solve the LTV LQ game in `lqg` for feedback Nash equilibrium via dynamic
# programming. Results written in-place into `strat`.
# """
# function solve_lq_game!(
#     strat::FeedbackStrategy{T},
#     lqg::LTVLQApproximation{T},
#     solver::PANGOLIN
# ) where {T}
#     N         = lqg.N
#     n         = lqg.n
#     n_p       = lqg.n_players
#     ctrl_dims = lqg.control_dims

#     total_controls = sum(ctrl_dims)
#     ctrl_ranges    = [sum(ctrl_dims[1:i-1])+1 : sum(ctrl_dims[1:i]) for i in 1:n_p]

#     state_reg   = T(solver.state_regularization)
#     control_reg = T(solver.control_regularization)

#     Z = [copy(lqg.Qf[i]) for i in 1:n_p]
#     ζ = [copy(lqg.qf[i]) for i in 1:n_p]

#     for k in N:-1:1
#         A_k = lqg.A[k]
#         d_k = lqg.d[k]

#         Q_k = [lqg.Q[i][k] + state_reg   * I for i in 1:n_p]
#         R_k = [lqg.R[i][k] + control_reg * I for i in 1:n_p]
#         M_k = [lqg.M[i][k]                   for i in 1:n_p]
#         q_k = [lqg.q[i][k]                   for i in 1:n_p]
#         r_k = [lqg.r[i][k]                   for i in 1:n_p]

#         B_k = hcat([lqg.B[i][k] for i in 1:n_p]...)

#         S  = zeros(T, total_controls, total_controls)
#         YP = zeros(T, total_controls, n)
#         Yα = zeros(T, total_controls)

#         for i in 1:n_p
#             rng_i = ctrl_ranges[i]
#             Bi    = lqg.B[i][k]
#             Ri    = R_k[i]
#             Zi    = Z[i]
#             Mi    = M_k[i]
#             BiZi  = Bi' * Zi

#             for j in 1:n_p
#                 rng_j = ctrl_ranges[j]
#                 Bj    = lqg.B[j][k]
#                 S[rng_i, rng_j] = (i == j ? Ri : zeros(T, length(rng_i), length(rng_j))) +
#                                   BiZi * Bj
#             end

#             YP[rng_i, :]  = (BiZi + Mi') * A_k
#             Yα[rng_i]     = Bi' * (ζ[i] + Zi * d_k) + r_k[i] + Mi' * d_k
#         end

#         if solver.check_singularity
#             S_rcond = 1.0 / cond(S)
#             if S_rcond < solver.rcond_threshold
#                 @warn "PANGOLIN solve_lq_game!: S poorly conditioned at k=$k" rcond=S_rcond
#             end
#         end

#         if state_reg > 0
#             for d in 1:total_controls
#                 S[d, d] += state_reg
#             end
#         end

#         P_full = S \ YP
#         α_full = S \ Yα

#         for i in 1:n_p
#             copyto!(strat.P[i][k], P_full[ctrl_ranges[i], :])
#             copyto!(strat.α[i][k], α_full[ctrl_ranges[i]])
#         end

#         F = A_k - B_k * P_full
#         β = -B_k * α_full

#         for i in 1:n_p
#             Pi  = P_full[ctrl_ranges[i], :]
#             αi  = α_full[ctrl_ranges[i]]
#             Ri  = R_k[i]
#             PRi = Pi' * Ri

#             Z_new = F' * Z[i] * F + Q_k[i] + PRi * Pi
#             Z[i]  = (Z_new + Z_new') / 2

#             ζ[i] = F' * (ζ[i] + Z[i] * (β + d_k)) + q_k[i] + PRi * αi - Pi' * r_k[i]
#         end
#     end

#     return strat
# end

# # ============================================================================
# # rollout!
# # ============================================================================

# """
#     rollout!(op, game, strat, x0; max_deviation=Inf) -> Bool

# Simulate forward under the feedback strategy. Writes results into `op` in-place.
# Returns `false` and aborts if ‖x(k+1) - x_ref(k+1)‖_∞ > `max_deviation`.
# """
# function rollout!(
#     op::OperatingPoint{T},
#     game::GameProblem{T},
#     strat::FeedbackStrategy{T},
#     x0::Vector{T};
#     max_deviation::T = T(Inf),
#     max_control::T   = T(Inf)
# ) where {T}
#     N   = length(op.x) - 1
#     n_p = game.n_players
#     dyn = game.dynamics

#     x_ref = [copy(op.x[k]) for k in 1:N+1]
#     op.x[1] .= x0

#     for k in 1:N
#         x_k = op.x[k]
#         for i in 1:n_p
#             op.u[i][k] .= .-strat.P[i][k] * x_k .- strat.α[i][k]
#         end

#         # Reject if any control magnitude exceeds the limit.
#         # This catches impulsive solutions generated by an ill-conditioned LQ
#         # subproblem when ρ is large: the backward DP produces a huge α[k] to
#         # dodge a constraint, and the deviation check alone won't catch it because
#         # the deviation is measured from the previous (already-impulsive) op.
#         if max_control < T(Inf)
#             for i in 1:n_p
#                 maximum(abs.(op.u[i][k])) > max_control && return false
#             end
#         end

#         x_next = _step_dynamics(dyn, x_k, [op.u[i][k] for i in 1:n_p], k)

#         if max_deviation < T(Inf) &&
#            maximum(abs.(x_next .- x_ref[k+1])) > max_deviation
#             return false
#         end

#         op.x[k+1] .= x_next
#     end
#     return true
# end

# function _step_dynamics(dyn::LinearDynamics{T}, x::Vector{T}, us::Vector{<:Vector}, k::Int) where {T}
#     x_next = get_A(dyn, k) * x
#     for i in eachindex(us)
#         x_next = x_next + get_B(dyn, i, k) * us[i]
#     end
#     return x_next
# end

# function _step_dynamics(dyn::CoupledNonlinearDynamics, x::Vector, us::Vector{<:Vector}, k::Int)
#     return dyn.func(x, vcat(us...), nothing, k)
# end

# function _step_dynamics(dyn::SeparableDynamics, x::Vector, us::Vector{<:Vector}, k::Int)
#     offs   = [0; cumsum(dyn.state_dims)[1:end-1]]
#     x_next = similar(x)
#     for i in eachindex(dyn.player_dynamics)
#         xi = x[offs[i]+1 : offs[i]+dyn.state_dims[i]]
#         x_next[offs[i]+1 : offs[i]+dyn.state_dims[i]] .=
#             dyn.player_dynamics[i](xi, us[i], nothing, k)
#     end
#     return x_next
# end

# # ============================================================================
# # backtrack_scale! — deviation-based line search
# # ============================================================================

# """
#     backtrack_scale!(strat, op, last_op, game, al_objs, solver, x0;
#                      first_iter, sf_init) -> (Bool, Float64)

# Line search on the feed-forward term α using a trajectory-deviation criterion:
# accept the first scale factor sf for which ‖x_new(k) - x_ref(k)‖_∞ < max_elwise_diff_step.

# See PANGOLIN docstring for why cost-based Armijo is inappropriate in noncooperative games.
# """
# function backtrack_scale!(
#     strat::FeedbackStrategy{T},
#     op::OperatingPoint{T},
#     last_op::OperatingPoint{T},
#     game::GameProblem{T},
#     al_objs,
#     solver::PANGOLIN,
#     x0::Vector{T};
#     first_iter::Bool = false,
#     sf_init::T = T(solver.α_init)
# ) where {T}
#     max_dev  = T(solver.max_elwise_diff_step)
#     max_ctrl = T(solver.control_limit)
#     N        = length(op.x) - 1
#     n_p      = game.n_players

#     sf = sf_init
#     for i in 1:solver.max_line_search_iter
#         if i == 1
#             scale_feedforward!(strat, sf)
#         else
#             scale_feedforward!(strat, T(solver.α_step))
#             sf *= T(solver.α_step)
#         end

#         for k in 1:N+1; op.x[k] .= last_op.x[k]; end
#         for p in 1:n_p, k in 1:N; op.u[p][k] .= last_op.u[p][k]; end

#         dev = (first_iter && i == 1) ? T(Inf) : max_dev
#         rollout!(op, game, strat, x0;
#                  max_deviation = dev,
#                  max_control   = max_ctrl) && return (true, sf)
#     end
#     return (false, zero(T))
# end

# # ============================================================================
# # has_converged
# # ============================================================================

# """
#     has_converged(op_new, op_old, tol) -> Bool

# ‖x_new(k) - x_old(k)‖_∞ < tol for all k.
# """
# function has_converged(op_new::OperatingPoint{T}, op_old::OperatingPoint{T}, tol::T) where {T}
#     return all(maximum(abs.(op_new.x[k] .- op_old.x[k])) < tol
#                for k in eachindex(op_new.x))
# end

# # ============================================================================
# # Warmstart / Initialization
# # ============================================================================

# function _init_operating_point(
#     game::GameProblem{T},
#     solver::PANGOLIN,
#     warmstart::Union{Nothing, WarmstartData},
#     cache::PANGOLINCache{T}
# ) where {T}
#     N    = n_steps(game)
#     n    = total_state_dim(game.dynamics)
#     n_p  = game.n_players
#     x0   = game.initial_state

#     strat = cache.strat_mem

#     if solver.init_mode == :warmstart && warmstart !== nothing &&
#        hasproperty(warmstart, :trajectories) && !isnothing(warmstart.trajectories)

#         op = OperatingPoint(warmstart.trajectories, n)

#         max_defect = zero(T)
#         for k in 1:N
#             x_next_true = _step_dynamics(game.dynamics, op.x[k],
#                                          [op.u[i][k] for i in 1:n_p], k)
#             max_defect = max(max_defect, maximum(abs.(x_next_true .- op.x[k+1])))
#         end

#         if max_defect > solver.dyn_consistency_tol
#             @warn "PANGOLIN: warmstart not dynamically consistent " *
#                   "(max_defect=$max_defect > tol=$(solver.dyn_consistency_tol)); " *
#                   "re-rolling out from x0 using warmstart controls"
#             _load_strategy_warmstart!(strat, warmstart, game, solver)
#             op_out = OperatingPoint(game)
#             rollout!(op_out, game, strat, x0)
#             return op_out
#         end

#         _load_strategy_warmstart!(strat, warmstart, game, solver)
#         return op

#     elseif solver.init_mode == :straight_line
#         op = OperatingPoint(game)
#         rollout!(op, game, strat, x0)
#         return op

#     else   # :zero_control
#         op = OperatingPoint(game)
#         rollout!(op, game, strat, x0)
#         return op
#     end
# end

# function _load_strategy_warmstart!(
#     strat::FeedbackStrategy{T},
#     warmstart::WarmstartData,
#     game::GameProblem{T},
#     solver::PANGOLIN
# ) where {T}
#     for i in eachindex(strat.P), k in eachindex(strat.P[i])
#         fill!(strat.P[i][k], zero(T))
#         fill!(strat.α[i][k], zero(T))
#     end
#     if solver.use_strategy_warmstart &&
#        hasproperty(warmstart, :feedback_gains) && !isnothing(warmstart.feedback_gains)
#         P_ws = warmstart.feedback_gains
#         α_ws = hasproperty(warmstart, :affine_gains) ? warmstart.affine_gains : nothing
#         N    = n_steps(game)
#         n_p  = game.n_players
#         for i in 1:n_p, k in 1:N
#             copyto!(strat.P[i][k], P_ws[i][k])
#             isnothing(α_ws) || copyto!(strat.α[i][k], α_ws[i][k])
#         end
#     end
# end

# # ============================================================================
# # AL objective construction
# # ============================================================================

# function _build_al_objectives(
#     game::GameProblem{T},
#     al_state::ALSolverState{T},
#     opts::ALOptions{T}
# ) where {T}
#     N   = n_steps(game)
#     n_p = game.n_players
#     objs = Vector{ALAugmentedObjective}(undef, n_p)
#     for i in 1:n_p
#         priv   = Vector{PrivateConstraint}(filter(pc -> pc.player == i,
#                                                    game.private_constraints))
#         shared = Vector{SharedConstraint}(filter(sc -> i in sc.players,
#                                                   game.shared_constraints))
#         shared_idxs    = findall(sc -> i in sc.players, game.shared_constraints)
#         λ_shared_traj_i = al_state.λ_shared_traj[shared_idxs]

#         objs[i] = ALAugmentedObjective(game.objectives[i], priv, shared,
#                                         λ_shared_traj_i, N, opts)
#     end
#     return objs
# end

# # ============================================================================
# # _solve — main loop
# # ============================================================================

# """
#     _solve(game::GameProblem{T}, solver::PANGOLIN, warmstart, verbose) -> GameSolution

# Main PANGOLIN solve loop.
# """
# function _solve(
#     game::GameProblem{T},
#     solver::PANGOLIN,
#     warmstart::Union{Nothing, WarmstartData},
#     verbose::Bool
# ) where {T}
#     @assert game.time_horizon isa DiscreteTime "PANGOLIN requires discrete-time formulation"

#     N   = n_steps(game)
#     n_p = game.n_players
#     x0  = game.initial_state

#     cache    = PANGOLINCache(game)
#     lqg      = cache.lq_mem
#     last_op  = cache.op_mem
#     strat    = cache.strat_mem

#     opts     = solver.constraint_opts
#     al_state = ALSolverState(Vector{SharedConstraint}(game.shared_constraints), opts, N)
#     al_objs  = _build_al_objectives(game, al_state, opts)

#     current_op = _init_operating_point(game, solver, warmstart, cache)
#     rollout!(current_op, game, strat, x0)

#     for k in 1:N+1; last_op.x[k] .= current_op.x[k]; end
#     for i in 1:n_p, k in 1:N; last_op.u[i][k] .= current_op.u[i][k]; end

#     iter_costs     = Vector{Vector{T}}()
#     violation_hist = Vector{T}()
#     ρ_hist         = Vector{T}()
#     converged      = false
#     i_iter         = 0
#     sf_adapt       = T(solver.α_init)
#     t_start        = time()

#     verbose && @info "PANGOLIN: starting" n_players=n_p N=N max_iter=solver.max_iter constrained=!is_unconstrained(game)

#     while !(converged || i_iter >= solver.max_iter)

#         # 1. LQ approximation
#         if is_unconstrained(game)
#             lq_approximation!(lqg, game, current_op)
#         else
#             _lq_approx_al!(lqg, game, al_objs, current_op)
#         end

#         # 2. Backward DP
#         solve_lq_game!(strat, lqg, solver)

#         # 3. Save last_op, line search
#         for k in 1:N+1; last_op.x[k] .= current_op.x[k]; end
#         for i in 1:n_p, k in 1:N; last_op.u[i][k] .= current_op.u[i][k]; end

#         objs_for_ls = is_unconstrained(game) ? nothing : al_objs
#         success, sf_accepted = backtrack_scale!(
#             strat, current_op, last_op, game, objs_for_ls, solver, x0;
#             first_iter = (i_iter == 0),
#             sf_init    = sf_adapt
#         )

#         if !success
#             verbose && @warn "PANGOLIN: line search failed at iter $i_iter"
#             for k in 1:N+1; current_op.x[k] .= last_op.x[k]; end
#             for i in 1:n_p, k in 1:N; current_op.u[i][k] .= last_op.u[i][k]; end
#             break
#         end

#         if sf_accepted ≥ sf_adapt * T(0.99)
#             sf_adapt = min(T(solver.α_init), sf_adapt * T(1.5))
#         else
#             sf_adapt = sf_accepted
#         end

#         # 4. Dual update
#         current_viol = zero(T)
#         if !is_unconstrained(game)
#             for i in 1:n_p
#                 update_multipliers!(al_objs[i], current_op, game)
#             end
#             _, current_viol = update_shared_multipliers!(
#                 al_state,
#                 Vector{SharedConstraint}(game.shared_constraints),
#                 al_objs, current_op
#             )
#             push!(violation_hist, current_viol)
#             push!(ρ_hist, al_state.ρ)

#             # 5. Violation-gated penalty increase (see al_augmentation.jl)
#             update_penalty!(al_state, al_objs, current_viol)
#         else
#             push!(violation_hist, zero(T))
#             push!(ρ_hist, zero(T))
#         end

#         # 6. Convergence check
#         # Require primal convergence AND feasibility AND cost stability.
#         # Cost stability guards against the period-2 cycle where ‖Δx‖ is small
#         # on the feasible half of the oscillation while costs swing by O(10+).
#         primal_converged = has_converged(current_op, last_op, T(solver.convergence_tol))
#         feas_converged   = is_unconstrained(game) ||
#                            (isempty(violation_hist) ? true :
#                             violation_hist[end] ≤ solver.constraint_opts.violation_tol)

#         curr_costs = _trajectory_costs(game, current_op)
#         cost_stable = length(iter_costs) < 2 || begin
#             prev_costs = iter_costs[end]
#             # Use a tight absolute tolerance: cost changes < 10× primal tolerance
#             # indicate a genuine fixed point rather than a low-amplitude oscillation.
#             maximum(abs.(curr_costs .- prev_costs)) ≤ T(solver.convergence_tol) * 10
#         end

#         converged = primal_converged && feas_converged && cost_stable

#         push!(iter_costs, curr_costs)
#         i_iter += 1

#         if verbose && (i_iter % 10 == 0 || i_iter == 1)
#             viol_str = is_unconstrained(game) ? "n/a" :
#                        string(round(violation_hist[end]; sigdigits=3))
#             @info "  iter $i_iter" converged=converged violation=viol_str ρ=al_state.ρ
#         end
#     end

#     solve_time = time() - t_start
#     verbose && @info "PANGOLIN: done" iterations=i_iter converged=converged solve_time=solve_time

#     trajectories = _op_to_trajectories(current_op, game)

#     solver_info = Dict{Symbol, Any}(
#         :feedback_gains    => deepcopy(strat.P),
#         :affine_gains      => deepcopy(strat.α),
#         :iter_costs        => iter_costs,
#         :violation_history => violation_hist,
#         :ρ_history         => ρ_hist,
#         :al_state          => al_state,
#         :iterations        => i_iter
#     )

#     return GameSolution(
#         game, trajectories;
#         equilibrium_type = :FeedbackNash,
#         converged        = converged,
#         iterations       = i_iter,
#         solve_time       = solve_time,
#         solver_info      = solver_info
#     )
# end

# # ============================================================================
# # Helpers
# # ============================================================================

# function _lq_approx_al!(
#     lqg::LTVLQApproximation{T},
#     game::GameProblem{T},
#     al_objs::Vector{<:ALAugmentedObjective},
#     op::OperatingPoint{T}
# ) where {T}
#     N   = lqg.N
#     n_p = lqg.n_players
#     dyn = game.dynamics

#     is_sep = dyn isa SeparableDynamics
#     state_offsets = game.metadata.state_offsets
#     state_dims    = game.metadata.state_dims

#     for k in 1:N
#         x_k = op.x[k]
#         u_k = [op.u[i][k] for i in 1:n_p]

#         A_k, B_vecs_k, d_k = linearize(dyn, x_k, u_k, k)
#         copyto!(lqg.A[k], A_k)
#         copyto!(lqg.d[k], d_k)
#         for i in 1:n_p; copyto!(lqg.B[i][k], B_vecs_k[i]); end

#         for i in 1:n_p
#             xi_k = is_sep ? x_k[state_offsets[i]+1 : state_offsets[i]+state_dims[i]] : x_k

#             player_cols = is_sep ?
#                 (state_offsets[i]+1 : state_offsets[i]+state_dims[i]) : nothing

#             Q_k, R_k, M_k, q_k, r_k, _ = quadraticize(
#                 al_objs[i], xi_k, op.u[i][k], k;
#                 x_joint = is_sep ? x_k : nothing,
#                 player_state_cols = player_cols
#             )

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
#         xi_N = is_sep ? x_N[state_offsets[i]+1 : state_offsets[i]+state_dims[i]] : x_N
#         Qf_i, qf_i = quadraticize(game.objectives[i].terminal_cost, xi_N)

#         u_dummy = zeros(T, game.metadata.control_dims[i])
#         ρ_i  = al_objs[i].ρ
#         for (j, pc) in enumerate(al_objs[i].private_constraints)
#             g  = evaluate_constraint(pc.constraint, xi_N, u_dummy, nothing, N+1)
#             λj = @view al_objs[i].λ_private[j][:, N+1]
#             Jx, _ = constraint_jacobian(pc.constraint, xi_N, u_dummy, nothing, N+1)
#             is_eq  = al_objs[i]._priv_eq[j]
#             n_rows = length(g)
#             for d in 1:n_rows
#                 active = is_eq || (g[d] + λj[d] / ρ_i > 0)
#                 active || continue
#                 μd = λj[d] + ρ_i * g[d]
#                 jx = @view Jx[d, :]
#                 Qf_i .+= ρ_i .* (jx * jx')
#                 qf_i .+= μd .* jx
#             end
#         end
#         Qf_i = (Qf_i + Qf_i') / 2

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

# function _trajectory_costs(game::GameProblem{T}, op::OperatingPoint{T}) where {T}
#     N   = length(op.x) - 1
#     n_p = game.n_players
#     costs = zeros(T, n_p)

#     is_sep        = game.dynamics isa SeparableDynamics
#     state_offsets = game.metadata.state_offsets
#     state_dims    = game.metadata.state_dims

#     for i in 1:n_p
#         sc = game.objectives[i].stage_cost
#         tc = game.objectives[i].terminal_cost

#         for k in 1:N
#             xi_k = is_sep ? op.x[k][state_offsets[i]+1 : state_offsets[i]+state_dims[i]] :
#                             op.x[k]
#             costs[i] += T(evaluate_stage_cost(sc, xi_k, op.u[i][k], nothing, k))
#         end

#         xi_N = is_sep ? op.x[N+1][state_offsets[i]+1 : state_offsets[i]+state_dims[i]] :
#                         op.x[N+1]
#         costs[i] += T(evaluate_terminal_cost(tc, xi_N, nothing))
#         costs[i] *= game.objectives[i].scaling
#     end
#     return costs
# end

# function _op_to_trajectories(op::OperatingPoint{T}, game::GameProblem{T}) where {T}
#     N      = length(op.x) - 1
#     n_p    = game.n_players
#     t_vec  = collect(range(T(0), game.time_horizon.tf, length=N+1))
#     costs  = _trajectory_costs(game, op)

#     is_sep        = game.dynamics isa SeparableDynamics
#     state_offsets = game.metadata.state_offsets
#     state_dims    = game.metadata.state_dims

#     trajectories = Trajectory{T}[]
#     for i in 1:n_p
#         if is_sep
#             ni     = state_dims[i]
#             rows   = state_offsets[i]+1 : state_offsets[i]+ni
#             states = hcat([op.x[k][rows] for k in 1:N+1]...)
#         else
#             states = hcat(op.x...)
#         end
#         ctrl_mat = hcat(op.u[i]...)
#         push!(trajectories, Trajectory(i, states, ctrl_mat, t_vec, costs[i]))
#     end
#     return trajectories
# end