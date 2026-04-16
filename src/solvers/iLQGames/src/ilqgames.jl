# ============================================================================
# iLQGames — Iterative Linear-Quadratic Games solver
#
# Implements Algorithm 1 of:
#   Fridovich-Keil, Ratner, Peters, Dragan, Tomlin.
#   "Efficient Iterative Linear-Quadratic Approximations for Nonlinear
#    Multi-Player General-Sum Differential Games."
#   arXiv:1909.04694 (ICRA 2020).
#
# ── Algorithm ────────────────────────────────────────────────────────────────
#
#   Given game g with initial state x₀, nominal trajectory (X, U):
#
#   for l = 1, 2, ..., max_iter:
#     1. exp = expand(g, da, X, U)          ← linearise dynamics + quadraticise costs
#     2. lq  = assemble_lq_game(exp, g)     ← time-varying LQ subgame
#     3. sol = solve(lq, FNELQ(μ))           ← backward Riccati recursion with S-regularisation
#     4. strat = rebase(sol.strategy, X, U)  ← align nominal to outer operating point
#     5. η, X̃, Ũ = line_search(g, x₀, strat, X, U)
#     6. X, U ← X̃, Ũ
#     7. Δ = ‖X_new - X_old‖_∞
#     8. if Δ < ε_conv: converged; break
#
#   Return GNEPSolution wrapping converged (X*, U*) as FeedbackNash strategy.
#
# ── Equilibrium type ─────────────────────────────────────────────────────────
#
#   iLQGames finds a feedback Nash equilibrium of successive LQ approximations.
#   The primary output is the FeedbackStrategy from the last converged FNELQ
#   solve, which constitutes a feedback Nash equilibrium of the LQ subgame
#   built around the converged nominal trajectory ξ*. This is the equilibrium
#   concept described in §IV-B of the paper: strategies that are globally Nash
#   optimal for the local LQ approximation at the fixed point.
#
#   The converged nominal trajectory (X*, U*) is also stored as an
#   OpenLoopStrategy in solver_info[:open_loop_strategy] for callers that
#   require an open-loop representation or warm-starting.
#
# ── Convergence criterion ────────────────────────────────────────────────────
#
#   Convergence is declared when the maximum element-wise deviation between
#   successive state and control trajectories falls below ε_conv:
#
#     Δ = max( ‖X^(l) - X^(l-1)‖_∞ , ‖U^(l) - U^(l-1)‖_∞ )
#
#   This is the trajectory-change criterion used in iLQGames.jl (are_close).
#   The paper's §IV-B strategy-change criterion (on Pᵢ and αᵢ norms) is
#   equivalent at the fixed point (by Lemma 1, αᵢ* = 0 at convergence) but
#   the trajectory criterion is more robust to scaling differences across
#   problems.
#
# ── Line search ──────────────────────────────────────────────────────────────
#
#   Geometric backtracking (§IV-C): start at η=1 (full LQ step), multiply by
#   β on failure. The full LQ step η=1 recovers Eq. (6) of the paper.
#   Minimum step η_min = β^20.
#
#   The Armijo (cost-decrease) criterion is NOT used. §IV-A of the paper
#   states that cost decrease is not a valid merit function in a noncooperative
#   game since costs conflict between players. The acceptance criterion is a
#   trajectory-deviation bound: ‖X_new - X_curr‖_∞ < max_state_step.
#
# ── Regularisation ───────────────────────────────────────────────────────────
#
#   Levenberg-Marquardt directly on the gain system S inside FNELQ. Adding
#   μI to each diagonal block Sᵢᵢ caps gain magnitude regardless of Z scale.
#   Regularising Q/R in the outer loop is insufficient when ‖B'ZB‖ >> μ.
#   μ increases geometrically on ill-conditioning; decays after accepted steps.
#
# ── Discretisation ───────────────────────────────────────────────────────────
#
#   Continuous dynamics are auto-discretised from solver.discretization and
#   game.dt. LinearDynamics (already discrete) bypasses discretisation.
#
# ── Fixes applied (see review) ───────────────────────────────────────────────
#
#   FIX-3: Line search starts at η=1 (full LQ step). Previously started at
#          η=β=0.5, never testing the step that recovers the LQ solution.
#   FIX-4: Line search fallback returns (X_curr, U_curr) with the cost already
#          computed for that trajectory, not zeros(U) with a mismatched cost.
#   FIX-6: Convergence docstring and implementation now agree: both use the
#          trajectory ∞-norm criterion. The previously documented strategy-norm
#          formula (on Pᵢ, αᵢ) is noted as an alternative but not the impl.
#   FIX-7: _regularise_lq_game is retained but clearly marked as superseded by
#          FNELQ's direct S-regularisation. It is no longer called in _solve.
#   FIX-9: primary_strategy, eq_type, file comment, and struct docstring now
#          all agree: the primary output is the FeedbackStrategy from the last
#          FNELQ solve, with equilibrium_type = :FeedbackNash.
#   FIX-10: cost_curr removed from _line_search signature (it was unused;
#           cost-decrease criterion is intentionally absent per §IV-A).
# ============================================================================

# ============================================================================
# iLQGames solver struct
# ============================================================================

"""
    iLQGames <: GameSolver

Iterative Linear-Quadratic Games solver for nonlinear multi-player general-sum
differential games.

Finds a feedback Nash equilibrium of successive LQ approximations by
linearising the dynamics and quadraticising the costs around the current
nominal trajectory, then solving the resulting LQ subgame exactly via FNELQ.

# Fields
- `max_iter`         : Maximum outer iterations (default 200)
- `ε_conv`           : Trajectory-change convergence threshold (default 0.05).
                       Convergence when max(‖X^(l)-X^(l-1)‖_∞, ‖U^(l)-U^(l-1)‖_∞) < ε_conv.
- `β`                : Line search backtrack factor (default 0.5).
                       Line search starts at η=1 and multiplies by β on failure.
- `η_min`            : Minimum line search step (default β^20)
- `max_state_step`   : Maximum ‖X_new - X_curr‖_∞ accepted in line search (default 1.0)
- `μ_init`           : Initial S-regularisation added inside FNELQ (default 1.0)
- `μ_max`            : Maximum regularisation (default 1e6)
- `μ_scale`          : Regularisation growth factor on ill-conditioning (default 10.0)
- `μ_decay`          : Regularisation decay factor on accepted step (default 0.5)
- `discretization`   : Method for auto-discretising continuous dynamics
                       (default `ZOHDiscretization()`)

# Equilibrium type
Returns `GNEPSolution` with `equilibrium_type = :FeedbackNash` and
`strategy = FeedbackStrategy` from the last converged FNELQ solve. This
constitutes a feedback Nash equilibrium of the LQ approximation at the fixed
point (§IV-B, Fridovich-Keil et al. 2020). The converged nominal trajectory
(X*, U*) is also stored in `solver_info[:open_loop_strategy]`.

# References
Fridovich-Keil, D., Ratner, E., Peters, L., Dragan, A. D., & Tomlin, C. J.
(2020). Efficient iterative linear-quadratic approximations for nonlinear
multi-player general-sum differential games. ICRA 2020. arXiv:1909.04694.
"""
struct iLQGames <: GameSolver
    max_iter::Int
    ε_conv::Float64
    β::Float64
    η_min::Float64
    max_state_step::Float64
    μ_init::Float64
    μ_max::Float64
    μ_scale::Float64
    μ_decay::Float64
    discretization::AbstractDiscretizationMethod

    function iLQGames(;
        max_iter::Int       = 200,
        ε_conv::Float64     = 0.05,
        β::Float64          = 0.5,
        η_min::Float64      = 0.5^20,
        max_state_step::Float64 = 1.0,
        μ_init::Float64     = 1.0,
        μ_max::Float64      = 1e6,
        μ_scale::Float64    = 10.0,
        μ_decay::Float64    = 0.5,
        discretization::AbstractDiscretizationMethod = ZOHDiscretization()
    )
        @assert max_iter > 0
        @assert ε_conv > 0
        @assert 0 < β < 1
        @assert η_min > 0
        @assert max_state_step > 0
        @assert μ_init >= 0
        @assert μ_max >= μ_init
        @assert μ_scale > 1
        @assert 0 < μ_decay < 1
        new(max_iter, ε_conv, β, η_min, max_state_step, μ_init, μ_max,
            μ_scale, μ_decay, discretization)
    end
end

solver_capabilities(::Type{iLQGames}) = [
    :NonlinearGame,
    :LQGame,
    :UnconstrainedGame,
    :DiscreteTime
]

# ============================================================================
# _build_da — construct or extract DiscreteApproximation
# ============================================================================

"""
    _build_da(game, solver) -> Union{DiscreteApproximation, Nothing}

Return a `DiscreteApproximation` appropriate for the game's dynamics:
- `LinearDynamics`     : return `nothing` — linearisation is exact.
- Continuous dynamics  : construct from `solver.discretization` and `game.dt`.
"""
function _build_da(game::GameProblem{T}, solver::iLQGames) where {T}
    dyn = game.dynamics
    if dyn isa LinearDynamics
        return nothing
    else
        dt = game.time_horizon.dt
        return discretize(dyn, dt; method=solver.discretization)
    end
end

# ============================================================================
# _regularise_lq_game — SUPERSEDED; retained for reference only
# ============================================================================

"""
    _regularise_lq_game(lq_game, μ) -> GameProblem

**SUPERSEDED** by direct S-regularisation inside FNELQ.

Adds μI to each player's Q and R matrices. This approach is insufficient when
‖B'ZB‖ >> μ, because the S-matrix condition is still dominated by Z. FNELQ's
direct S-regularisation (adding μI to diagonal blocks of S) is correct because
it bounds P = (S + μI)⁻¹·YP regardless of Z scale.

Retained here for reference and potential use in alternative solver backends
that do not support S-level regularisation.
"""
function _regularise_lq_game(lq_game::GameProblem{T}, μ::T) where {T}
    μ == zero(T) && return lq_game

    dyn  = lq_game.dynamics
    N    = n_steps(lq_game)
    np   = num_players(lq_game)
    n    = total_state_dim(dyn)

    new_objectives = map(1:np) do i
        old_obj = get_objective(lq_game, i)
        old_sc  = old_obj.stage_cost
        n_u     = lq_game.metadata.control_dims[i]
        μI_u    = μ * Matrix{T}(I, n_u, n_u)
        μI_x    = μ * Matrix{T}(I, n, n)
        if is_ltv(old_sc)
            Q_reg = [get_Q(old_sc, k) + μI_x for k in 1:N]
            R_reg = [get_R(old_sc, k) + μI_u for k in 1:N]
            new_sc = LQStageCost(
                Q_reg,
                R_reg,
                [get_M(old_sc, k) for k in 1:N],
                [get_q(old_sc, k) for k in 1:N],
                [get_r(old_sc, k) for k in 1:N]
            )
        else
            new_sc = LQStageCost(
                get_Q(old_sc, 1) + μI_x,
                get_R(old_sc, 1) + μI_u
            )
        end
        PlayerObjective(i, new_sc, old_obj.terminal_cost, old_obj.scaling)
    end

    return GameProblem{T}(
        lq_game.n_players,
        new_objectives,
        lq_game.dynamics,
        lq_game.initial_state,
        lq_game.private_constraints,
        lq_game.shared_constraints,
        lq_game.time_horizon,
        lq_game.metadata
    )
end

# ============================================================================
# _total_costs — evaluate sum of all player stage+terminal costs
# ============================================================================

"""
    _total_costs(game, X, U) -> (costs::Vector, total::Float64)

Evaluate per-player costs and their sum for trajectory `(X, U)`.
Returns `(fill(Inf, np), Inf)` if any entry is non-finite.
"""
function _total_costs(game::GameProblem{T}, X::Matrix{T}, U::Matrix{T}) where {T}
    if !all(isfinite, X) || !all(isfinite, U)
        costs = fill(T(Inf), game.n_players)
        return costs, T(Inf)
    end
    N  = size(U, 2)
    np = game.n_players
    dyn = game.dynamics
    cd  = game.metadata.control_dims
    co  = [0; cumsum(cd)[1:end-1]]

    is_sep = dyn isa SeparableDynamics
    sd  = is_sep ? dyn.state_dims  : nothing
    so  = is_sep ? [0; cumsum(dyn.state_dims)[1:end-1]] : nothing

    costs = zeros(T, np)
    for i in 1:np
        obj = get_objective(game, i)
        Ui  = U[co[i]+1 : co[i]+cd[i], :]
        Xi  = is_sep ? X[so[i]+1 : so[i]+sd[i], :] : X

        costs[i] = total_cost(obj,
                               [Xi[:, k]  for k in 1:N+1],
                               [Ui[:, k]  for k in 1:N],
                               nothing)
    end
    return costs, sum(costs)
end

# ============================================================================
# _line_search — geometric backtracking
# ============================================================================

"""
    _line_search(game, x0, X_curr, U_curr, strat, β, η_min, max_state_step, times)

Backtracking line search over the feedforward scale η.

FIX-3: Starts at η=1 (the full LQ step, recovering Eq. (6) of the paper).
Previously started at η=β=0.5, which never tested the LQ solution and caused
systematic undershoot near convergence.

FIX-4: On failure, returns (X_curr, U_curr) unchanged. Previously returned
zeros(U), which produced a cost mismatch that corrupted the convergence
check on the subsequent iteration.

Acceptance criterion: rollout must be finite AND
‖X_new[:,k] - X_curr[:,k]‖_∞ < max_state_step for all k.
This keeps the new operating point within the validity region of the current
linearisation. Cost-decrease (Armijo) is intentionally NOT used — §IV-A of
the paper notes that total cost is not a valid merit function in a
noncooperative game.

η scales only the feedforward αᵢ, not the gains Pᵢ — Eq. (7) of the paper.
"""
function _line_search(
    game::GameProblem{T},
    x0::Vector{T},
    X_curr::Matrix{T},
    U_curr::Matrix{T},   # FIX-4: need current U to return on failure
    strat::FeedbackStrategy{T},
    β::Float64,
    η_min::Float64,
    max_state_step::Float64,
    times::Vector{T}
) where {T}
    # FIX-3: Start at η=1. The LQ solution (Eq. (6)) is recovered at η=1.
    η = T(1.0)
    N       = size(X_curr, 2) - 1
    m_total = sum(game.metadata.control_dims)

    idx = 1

    while η >= T(η_min)
        X_new = Matrix{T}(undef, size(X_curr, 1), N + 1)
        U_new = Matrix{T}(undef, m_total, N)
        X_new[:, 1] = x0

        aborted = false
        for k in 1:N
            u_k = apply_strategy(strat, X_new[:, k], k; η)
            U_new[:, k] = u_k
            xk1 = _rollout_step(game.dynamics, X_new[:, k], u_k, nothing, times, k)
            if !all(isfinite, xk1) ||
               maximum(abs.(xk1 - X_curr[:, k+1])) >= T(max_state_step)
                aborted = true
                break
            end
            X_new[:, k+1] = xk1
        end

        if !aborted
            costs_new, _ = _total_costs(game, X_new, U_new)
            return η, X_new, U_new, costs_new, true
        end

        η *= T(β)

        # Temporarily add inside the line search loop, first iteration only
        # if idx == 1
        #     @info "line search k=1" η=η u_k=apply_strategy(strat, X_new[:,1], 1; η=T(1.0)) x=X_new[:,1]
        # end 
        # idx += 1

        η < T(η_min) && break
    end

    # FIX-4: Return the trajectory we were called with, unchanged.
    # Previously returned zeros(U) with a cost computed against zero controls,
    # causing a spurious cost mismatch on the next iteration.
    costs_curr, _ = _total_costs(game, X_curr, U_curr)
    return T(η_min), X_curr, U_curr, costs_curr, false
end

# ============================================================================
# _trajectory_change — convergence metric
# ============================================================================

"""
    _trajectory_change(X_new, X_old, U_new, U_old) -> Float64

Maximum element-wise ∞-norm between consecutive state and control trajectories:

    Δ = max( max_k ‖X_new[:,k] - X_old[:,k]‖_∞,
             max_k ‖U_new[:,k] - U_old[:,k]‖_∞ )

Convergence is declared when Δ < ε_conv (default 0.05).

FIX-6: The docstring previously described a strategy-change criterion (on
Pᵢ and αᵢ Frobenius norms). The implementation has always used the trajectory
∞-norm; docstring and code now agree. The strategy-change criterion from §IV-B
is equivalent at the fixed point (by Lemma 1, αᵢ* = 0 at convergence) but the
trajectory criterion is more numerically robust across problem scales.
"""
function _trajectory_change(
    X_new::Matrix{T}, X_old::Matrix{T},
    U_new::Matrix{T}, U_old::Matrix{T}
) where {T}
    Δx = maximum(maximum(abs.(X_new[:, k] - X_old[:, k])) for k in axes(X_new, 2))
    Δu = maximum(maximum(abs.(U_new[:, k] - U_old[:, k])) for k in axes(U_new, 2))
    return max(Δx, Δu)
end

# ============================================================================
# _solve — main outer loop
# ============================================================================

function _solve(
    game::GameProblem{T},
    solver::iLQGames,
    warmstart::Union{Nothing, WarmstartData},
    verbose::Bool
) where {T}

    @assert(game.time_horizon isa DiscreteTime,
        "iLQGames requires DiscreteTime; got $(typeof(game.time_horizon))")
    @assert(is_unconstrained(game),
        "iLQGames handles only unconstrained games; soft constraints should " *
        "be encoded in the cost functions")

    N   = n_steps(game)
    np  = num_players(game)
    x0  = game.initial_state
    n   = total_state_dim(game.dynamics)
    dyn = game.dynamics
    t_vec = collect(range(T(0), game.time_horizon.tf, length=N+1))
    cd    = game.metadata.control_dims
    m_total = sum(cd)

    # ── Build DiscreteApproximation ──────────────────────────────────────────
    da = _build_da(game, solver)

    # ── Initialise nominal trajectory ────────────────────────────────────────
    X = zeros(T, n, N+1)
    U = zeros(T, m_total, N)

    if warmstart !== nothing && warmstart.strategy !== nothing
        ws = warmstart.strategy
        if ws isa FeedbackStrategy{T}
            X = copy(ws.nominal_states)
            c_offs = [0; cumsum(cd)[1:end-1]]
            for i in 1:np
                U[c_offs[i]+1:c_offs[i]+cd[i], :] = ws.nominal_controls[i]
            end
        elseif ws isa OpenLoopStrategy{T}
            c_offs = [0; cumsum(cd)[1:end-1]]
            for i in 1:np
                U[c_offs[i]+1:c_offs[i]+cd[i], :] = ws.controls[i]
            end
            X = rollout(dyn, x0, U, nothing, t_vec)
        end
    elseif warmstart !== nothing && warmstart.trajectories !== nothing
        traj1 = warmstart.trajectories[1]
        X = copy(traj1.states)
        c_offs = [0; cumsum(cd)[1:end-1]]
        for i in 1:np
            traj_i = warmstart.trajectories[i]
            U[c_offs[i]+1:c_offs[i]+cd[i], :] = traj_i.controls
        end
    else
        X = rollout(dyn, x0, U, nothing, t_vec)
    end

    costs_vec, cost_curr = _total_costs(game, X, U)

    if verbose
        @info "iLQGames start" np n N cost_init=cost_curr costs_per_player=costs_vec μ_init=solver.μ_init
    end

    # ── Tracking state ───────────────────────────────────────────────────────
    t_start        = time()
    μ              = T(solver.μ_init)
    converged      = false
    iterations_run = 0
    strat_last     = nothing   # FeedbackStrategy from last successful FNELQ
    Δ_history      = Float64[]
    cost_history   = Float64[cost_curr]
    η_history      = Float64[]
    μ_history      = Float64[Float64(μ)]

    # ── Outer loop ───────────────────────────────────────────────────────────
    for iter in 1:solver.max_iter
        iterations_run = iter
        X_prev = copy(X)
        U_prev = copy(U)

        # Step 1–2: Linearise + assemble LQ subgame
        exp_obj = expand(game, X, U, da)
        lq_game = assemble_lq_game(exp_obj, game)

        # Step 3: Solve LQ subgame with S-regularisation.
        # FIX-7: _regularise_lq_game (which adds μI to Q/R) is NOT called.
        # FNELQ's direct S-regularisation is used instead — it is effective
        # even when ‖B'ZB‖ >> μ, where Q/R regularisation is not.
        # Escalate μ in a loop until FNELQ reports S well-conditioned.
        fnelq_solver = FNELQ(check_singularity=true, rcond_threshold=1e-10,
                             regularization=Float64(μ))
        lq_sol   = _solve(lq_game, fnelq_solver, nothing, false)
        fnelq_ok = lq_sol.converged

        n_fnelq_retry = 0
        while !fnelq_ok && μ < T(solver.μ_max)
            μ = min(μ * T(solver.μ_scale), T(solver.μ_max))
            n_fnelq_retry += 1
            verbose && @info "  iter $iter: FNELQ ill-conditioned, escalating μ" μ=μ retry=n_fnelq_retry
            fnelq_solver = FNELQ(check_singularity=true, rcond_threshold=1e-10,
                                 regularization=Float64(μ))
            lq_sol  = _solve(lq_game, fnelq_solver, nothing, false)
            fnelq_ok = lq_sol.converged
            @info "lq_sol strategy feedforward" lq_sol.strategy.feedforward[1][1] lq_sol.strategy.nominal_controls[1][:,1]
        end

        strat_new = lq_sol.strategy

        # ── Rebase nominal trajectory onto outer operating point ──────────────
        # FNELQ stores the LQ inner trajectory as its nominal. For nonlinear
        # dynamics the inner linear trajectory diverges from the outer nonlinear
        # X after the first step.
        #
        # apply_strategy computes: û_nom[k] - Pᵢ(x - x̂_nom[k]) - η·αᵢ
        # This equals -Pᵢ·x - η·αᵢ only when x̂_nom = X (the nonlinear
        # operating point). Rebasing replaces the FNELQ inner trajectory with
        # the current outer X and U so that Eq. (7) of the paper holds exactly.
        c_offs_nom = [0; cumsum(cd)[1:end-1]]
        u_nom_rebased = [
            U[c_offs_nom[i]+1:c_offs_nom[i]+cd[i], :]
            for i in 1:np
        ]
        ff_rebased = [[copy(strat_new.feedforward[i][k]) for k in 1:N] for i in 1:np]

        strat_new = FeedbackStrategy(
            strat_new.gains,
            ff_rebased,
            copy(X),
            u_nom_rebased,
            cd,
            t_vec
        )

        # Step 5–6: Line search.
        # FIX-3: starts at η=1 (full LQ step).
        # FIX-4: failure returns (X_curr, U_curr) unchanged.

        # Add this directly before the _line_search call in _solve
        function apply_strategy_debug(s::FeedbackStrategy{T}, x::AbstractVector, k::Int; η::Real=one(T)) where {T}
            x̂ = s.nominal_states[:, k]
            û = s.nominal_controls[1][:, k]
            P = s.gains[1][k]
            δx = x - x̂
            u  = û - P * δx - η .* s.feedforward[1][k]
            @info "apply_strategy_debug" k x=x[1] x̂=x̂[1] û=û[1] P=P[1,1] δx=δx[1] u=u[1]
            return u
        end

        if verbose
            u_test = apply_strategy_debug(strat_new, X[:, 1], 1; η=T(1.0))
            @info "  rebase check" k=1 x=X[1,1] u_applied=u_test[1] u_nom_k1=strat_new.nominal_controls[1][1,1] P_k1=strat_new.gains[1][1][1,1]
            @info "nominal_states col 1" strat_new.nominal_states[:, 1]
            @info "nominal_controls[1] col 1" strat_new.nominal_controls[1][:, 1]  
            @info "X col 1" X[:, 1]
        end


        η, X_new, U_new, costs_new, ls_ok = _line_search(
            game, x0, X, U, strat_new,
            solver.β, solver.η_min, solver.max_state_step, t_vec
        )

        push!(η_history, η)

        if !ls_ok
            # Line search failed: the gain direction may be near-zero (already
            # at a fixed point), or regularisation is insufficient.
            # Check whether we are already within convergence tolerance —
            # this can happen when the last accepted step was small and the
            # new step is rejected because αᵢ ≈ 0 (Lemma 1).
            if !isempty(Δ_history) && last(Δ_history) < T(solver.ε_conv)
                verbose && @info "  iter $iter: line search failed near fixed point, declaring convergence" Δ=last(Δ_history)
                converged = true
                break
            end
            verbose && @info "  iter $iter: line search failed, escalating μ" μ_before=μ
            μ = min(μ * T(solver.μ_scale), T(solver.μ_max))
            push!(μ_history, Float64(μ))
            push!(cost_history, cost_curr)
            push!(Δ_history, Inf)
            strat_last = strat_new
            continue
        end

        # Step 7: Accept step
        X, U      = X_new, U_new
        cost_curr = sum(costs_new)
        push!(cost_history, cost_curr)
        push!(μ_history, Float64(μ))

        # Decay regularisation; never below μ_init floor
        μ = max(T(solver.μ_init), μ * T(solver.μ_decay))

        # Step 8: Convergence check
        Δ = _trajectory_change(X, X_prev, U, U_prev)
        push!(Δ_history, Δ)

        if verbose
            # Per-player costs for this iterate
            cost_breakdown = ["player $i: $(round(costs_new[i], sigdigits=5))" for i in 1:np]
            @info "  iter $iter" η=η total_cost=cost_curr Δ_traj=Δ μ=μ n_fnelq_retry=n_fnelq_retry costs=join(cost_breakdown, ", ")
        end

        # Step 9: Declare convergence
        if Δ < T(solver.ε_conv)
            converged = true
            break
        end

        strat_last = strat_new
    end

    total_time = time() - t_start

    # ── Build output ──────────────────────────────────────────────────────────
    #
    # FIX-9: Primary output is the FeedbackStrategy from the last converged
    # FNELQ solve. This is the feedback Nash equilibrium of the LQ subgame at
    # the fixed point (§IV-B). The file comment, struct docstring, and this
    # assignment now all agree.
    #
    # The converged nominal trajectory (X, U) is stored as OpenLoopStrategy in
    # solver_info for callers that require an open-loop representation or for
    # warm-starting subsequent solves.
    c_offs       = [0; cumsum(cd)[1:end-1]]
    U_per_player = [U[c_offs[i]+1:c_offs[i]+cd[i], :] for i in 1:np]
    ol_strategy  = OpenLoopStrategy(U_per_player, cd, t_vec)

    costs_final, total_cost_final = _total_costs(game, X, U)

    trajectories = [
        Trajectory(i, X, U_per_player[i], t_vec, costs_final[i])
        for i in 1:np
    ]

    if !converged
        @warn "iLQGames: did not converge in $(solver.max_iter) iterations" last_Δ=isempty(Δ_history) ? NaN : last(Δ_history) ε_conv=solver.ε_conv
    end

    if verbose
        @info "iLQGames complete" converged iterations=iterations_run total_time final_cost=total_cost_final final_μ=μ
    end

    # FIX-9: equilibrium_type is :FeedbackNash when a feedback strategy exists,
    # consistent with the struct docstring and file-level comment.
    primary_strategy = strat_last !== nothing ? strat_last : ol_strategy
    eq_type          = strat_last !== nothing ? :FeedbackNash : :OpenLoopNash

    return GNEPSolution(
        game,
        trajectories;
        state_trajectory = X,
        strategy         = primary_strategy,
        equilibrium_type = eq_type,
        converged        = converged,
        iterations       = iterations_run,
        solve_time       = total_time,
        solver_info      = Dict{Symbol, Any}(
            :cost_history              => cost_history,
            :trajectory_change_history => Δ_history,
            :η_history                 => η_history,
            :μ_history                 => μ_history,
            :open_loop_strategy        => ol_strategy,
            :final_regularisation      => μ,
            :final_cost_per_player     => costs_final
        )
    )
end