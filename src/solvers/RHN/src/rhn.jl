# ============================================================================
# solvers/RHN/src/rhn.jl
#
# Receding Horizon Nash (RHN) solver wrapper.
#
# Wraps any GameSolver to produce a closed-loop Nash strategy via the
# receding-horizon loop:
#
#   for t = 1, …, N_sim:
#     1. sub   = with_initial_state(game_win, X[:,t])   ← O(1) patch
#     2. sol   = solve(sub, inner_solver; warmstart=ws)  ← any solver
#     3. u_t   = apply_strategy(sol.strategy, X[:,t], 1) ← first action
#     4. X[:,t+1] = _rollout_step(dyn, X[:,t], u_t, …)  ← propagate
#     5. ws    = _shift_warmstart(sol)                   ← shift by 1 step
#
# The shift-and-warm-start step follows the TinyMPC / Mattingley et al.
# convention: drop the first element of the solution, repeat the last,
# and use this as the initial iterate for the next sub-problem solve.
#
# References:
#   Mattingley, Wang, Boyd (2011) — Receding Horizon Control.
#   Mahajan et al. (2026)        — Conic-TinyMPC (shift warm-start).
#   Laine et al. (2023)          — GFNE (generalized feedback Nash).
#
# Solvers package imports (available in scope):
#   GameSolver, solve, _solve, WarmstartData, solver_capabilities
#   + all DifferentialGamesBase exports
# ============================================================================

# ============================================================================
# RecedingHorizonNash solver struct
# ============================================================================

"""
    RecedingHorizonNash{S <: GameSolver} <: GameSolver

Receding-horizon wrapper for any Nash equilibrium solver.

At each simulation step the sub-game is solved over a fixed prediction
window, only the first optimal control is applied, and the previous
solution is shifted to warm-start the next solve.

# Fields
- `inner_solver`  : Any `GameSolver` that handles `GameProblem{T}`
  (FNELQ, iLQGames, ALGAMES, etc.).
- `warm_start`    : Whether to shift the previous solution as a warm-start
  for the next sub-problem (default `true`).
- `verbose_inner` : Whether to pass `verbose=true` to inner solves
  (default `false`).

# Construction
```julia
# Minimal — wraps any solver
solver = RecedingHorizonNash(FNELQ())
solver = RecedingHorizonNash(iLQGames(); warm_start=true, verbose_inner=false)
solver = RecedingHorizonNash(ALGAMES())
```

# Usage
```julia
prob = RecedingHorizonNashProblem(horizon_game, x0, 50)
sol  = solve(prob, RecedingHorizonNash(FNELQ()))
```

# Warm-start
After each sub-solve, the strategy is shifted by one step:
- `FeedbackStrategy`: gains, feedforward, and nominal trajectory are
  all shifted left by 1 step; the last element is repeated.
- `OpenLoopStrategy`: control sequences are shifted left by 1 step;
  the last control is repeated.
This provides a hot-start that is valid (feasible dynamics) and near
the new optimal point, reducing inner solve iterations significantly
for nonlinear solvers (iLQGames, ALGAMES).

For FNELQ the warm-start is structurally accepted but has no effect
on the solution since FNELQ solves exactly in one backward pass.
"""
struct RecedingHorizonNash{S <: GameSolver} <: GameSolver
    inner_solver ::S
    warm_start   ::Bool
    verbose_inner::Bool

    function RecedingHorizonNash(
        inner_solver ::S;
        warm_start   ::Bool = true,
        verbose_inner::Bool = false
    ) where {S <: GameSolver}
        new{S}(inner_solver, warm_start, verbose_inner)
    end
end

solver_capabilities(::Type{<:RecedingHorizonNash}) = [:RecedingHorizon, :FeedbackPolicies]

# ============================================================================
# Warm-start shift helpers
# ============================================================================

"""
    _shift_feedback_strategy(strat::FeedbackStrategy{T}) -> FeedbackStrategy{T}

Shift a feedback strategy forward by one step for warm-starting.

Implements the TinyMPC "shift and repeat last" heuristic:
- Gains Pᵢ(k): [P(2), P(3), …, P(N), P(N)]
- Feedforward αᵢ(k): similarly shifted
- Nominal states x̂:  [x̂₂, x̂₃, …, x̂_{N+1}, x̂_{N+1}]
- Nominal controls ûᵢ: [û₂, …, û_N, û_N]
- Times: unchanged (window-local [0, dt, 2dt, …])

At the new t=1, the repeated last element is the best available
guess for the final step of the new horizon — conservative but valid.
"""
function _shift_feedback_strategy(strat::FeedbackStrategy{T}) where {T}
    N  = n_steps(strat)
    np = strat.n_players

    # Gains and feedforward: shift left, repeat last
    new_gains = [
        [strat.gains[i][min(k+1, N)] for k in 1:N]
        for i in 1:np
    ]
    new_ff = [
        [strat.feedforward[i][min(k+1, N)] for k in 1:N]
        for i in 1:np
    ]

    # Nominal states: columns 2…N+1, then repeat N+1
    new_nom_x = hcat(
        strat.nominal_states[:, 2:N+1],       # columns 2 to N+1 (shift left)
        strat.nominal_states[:, N+1:N+1]       # repeat last state
    )

    # Nominal controls per player: columns 2…N, then repeat N
    new_nom_u = [
        if N > 1
            hcat(strat.nominal_controls[i][:, 2:N],
                 strat.nominal_controls[i][:, N:N])
        else
            copy(strat.nominal_controls[i])  # N=1: no shift possible
        end
        for i in 1:np
    ]

    return FeedbackStrategy(
        new_gains, new_ff,
        new_nom_x, new_nom_u,
        strat.control_dims,
        strat.times          # window-local times are unchanged
    )
end

"""
    _shift_open_loop_strategy(strat::OpenLoopStrategy{T}) -> OpenLoopStrategy{T}

Shift an open-loop strategy forward by one step for warm-starting.

Controls are shifted left by 1 step; the last control is repeated:
  ûᵢ(k) ← ûᵢ(k+1) for k = 1…N-1;  ûᵢ(N) ← ûᵢ(N).
"""
function _shift_open_loop_strategy(strat::OpenLoopStrategy{T}) where {T}
    N  = n_steps(strat)
    np = strat.n_players

    new_controls = [
        if N > 1
            hcat(strat.controls[i][:, 2:N],
                 strat.controls[i][:, N:N])
        else
            copy(strat.controls[i])
        end
        for i in 1:np
    ]

    return OpenLoopStrategy(new_controls, strat.control_dims, strat.times)
end

"""
    _shift_warmstart(sol::GNEPSolution{T}) -> WarmstartData{T}

Build a warm-start for the next sub-problem from a converged solution.

The strategy (if present) is shifted forward by one step.
Dual variables are forwarded as-is so constrained solvers (ALGAMES)
can continue from their previous multiplier estimates.
"""
function _shift_warmstart(sol::GNEPSolution{T}) where {T}
    shifted_strat = if sol.strategy isa FeedbackStrategy
        _shift_feedback_strategy(sol.strategy)
    elseif sol.strategy isa OpenLoopStrategy
        _shift_open_loop_strategy(sol.strategy)
    else
        nothing
    end

    # Carry dual variables (ALGAMES uses these for warm-starting multipliers)
    duals = get(sol.solver_info, :dual_variables, nothing)

    return WarmstartData{T}(nothing, shifted_strat, duals, Dict{Symbol, Any}())
end

# ============================================================================
# Cost accumulation helper
# ============================================================================

"""
    _compute_rh_costs(game, X, U) -> Vector{T}

Accumulate stage and terminal costs over the full closed-loop simulation.

For each player i:
  Jᵢ = Σ_{t=1}^{N_sim} scaling_i · ℓᵢ(X[:,t], U[i][:,t]) + scaling_i · φᵢ(X[:,end])

Uses the objectives from `game` (the horizon template game). For LQ
games, `ℓᵢ` is the quadratic stage cost and `φᵢ` is the quadratic
terminal cost. Works for any `AbstractStageCost`/`AbstractTerminalCost`
subtype via the `evaluate_stage_cost`/`evaluate_terminal_cost` interface.

`U[i]` must be (mᵢ × N_sim) with player i's own controls.
"""
function _compute_rh_costs(
    game::GameProblem{T},
    X   ::Matrix{T},
    U   ::Vector{Matrix{T}}
) where {T}
    np    = n_players(game)
    N_sim = size(X, 2) - 1
    costs = zeros(T, np)

    for i in 1:np
        obj = get_objective(game, i)
        sc  = obj.stage_cost
        tc  = obj.terminal_cost
        s   = obj.scaling

        # Stage cost sum
        for t in 1:N_sim
            costs[i] += evaluate_stage_cost(sc, X[:, t], U[i][:, t], nothing, t)
        end

        # Terminal cost at final state
        costs[i] += evaluate_terminal_cost(tc, X[:, N_sim+1], nothing)

        # Apply player scaling
        costs[i] *= s
    end

    return costs
end

# ============================================================================
# _solve — main receding-horizon loop
# ============================================================================

"""
    _solve(prob, solver, warmstart, verbose) -> RecedingHorizonNashSolution{T}

Receding-horizon Nash loop.

At each step t ∈ 1:n_sim_steps:
  1. Patch sub-problem: `sub = with_initial_state(horizon_game, X[:,t])`  [O(1)]
  2. Solve sub-problem: `sol = solve(sub, inner_solver; warmstart=ws)`
  3. Extract first control: `u_t = apply_strategy(sol.strategy, X[:,t], 1)`
     Fallback: reads from trajectory if strategy is nothing.
  4. Advance state: `X[:,t+1] = _rollout_step(dyn, X[:,t], u_t, …, 1)`
  5. Update warm-start: `ws = _shift_warmstart(sol)` (if warm_start=true)

The `warmstart` argument is the initial warm-start for step t=1 (typically nothing).
"""
function _solve(
    prob     ::RecedingHorizonNashProblem{T},
    solver   ::RecedingHorizonNash,
    warmstart::Union{Nothing, GNEPSolution, WarmstartData},
    verbose  ::Bool
) where {T}
    game   = prob.horizon_game
    np     = n_players(game)
    n      = state_dim(game)
    N_sim  = prob.n_sim_steps
    dyn    = game.dynamics

    # Control dimension bookkeeping
    cdims   = [control_dim(game, i) for i in 1:np]
    offsets = [0; cumsum(cdims)[1:end-1]]

    # Window-local time vector (reused every step — never changes)
    th = game.time_horizon
    @assert th isa DiscreteTime "RecedingHorizonNash requires a DiscreteTime horizon"
    times_w = collect(range(zero(T), T(th.tf), length=n_steps(game)+1))

    # Allocate output buffers
    X         = Matrix{T}(undef, n, N_sim+1)
    X[:, 1]   = prob.x0
    U_players = [Matrix{T}(undef, cdims[i], N_sim) for i in 1:np]

    # Normalise initial warm-start to WarmstartData
    ws::Union{Nothing, WarmstartData} = if warmstart isa GNEPSolution
        WarmstartData(warmstart)
    else
        warmstart  # WarmstartData or nothing
    end

    total_time    = 0.0
    n_converged   = 0
    all_converged = true

    verbose && @info "RHN: starting $(N_sim)-step simulation" inner=typeof(solver.inner_solver)

    for t in 1:N_sim
        # ── 1. Patch sub-problem initial condition ─────────────────────────────
        subgame = with_initial_state(game, X[:, t])

        # ── 2. Solve sub-problem ───────────────────────────────────────────────
        t0 = time()
        inner_sol = solve(
            subgame, solver.inner_solver;
            warmstart           = ws,
            verbose             = solver.verbose_inner,
            check_compatibility = false
        )
        total_time += time() - t0

        if inner_sol.converged
            n_converged += 1
        else
            all_converged = false
            verbose && @warn "RHN: inner solve did not converge at step $t"
        end

        # ── 3. Extract first joint control ────────────────────────────────────
        u_joint = if inner_sol.strategy !== nothing
            # Evaluate feedback law u = û(1) - P(1)·(x - x̂(1)) + η·α(1)
            # Since subgame was built with x̂(1) = X[:,t], δx = 0 at this step
            apply_strategy(inner_sol.strategy, X[:, t], 1)
        else
            # Fallback: first step of each player's stored trajectory
            vcat([get_trajectory(inner_sol, i).controls[:, 1] for i in 1:np]...)
        end

        # Store per-player controls
        for i in 1:np
            U_players[i][:, t] = u_joint[offsets[i]+1 : offsets[i]+cdims[i]]
        end

        # ── 4. Advance state ───────────────────────────────────────────────────
        # _rollout_step is exact for LinearDynamics; RK4 for nonlinear.
        X[:, t+1] = _rollout_step(dyn, X[:, t], u_joint, nothing, times_w, 1)

        # ── 5. Shift warm-start for next iteration ─────────────────────────────
        ws = solver.warm_start ? _shift_warmstart(inner_sol) : nothing

        if verbose && (t % max(1, N_sim ÷ 5) == 0 || t == 1 || t == N_sim)
            @info "  RHN t=$t/$N_sim" x_norm=norm(X[:, t+1]) converged=inner_sol.converged
        end
    end

    # ── Accumulate costs over the full closed-loop trajectory ─────────────────
    costs = _compute_rh_costs(game, X, U_players)

    verbose && @info "RHN complete" total_time=round(total_time, digits=3) converged=all_converged

    return RecedingHorizonNashSolution{T}(
        prob, X, U_players, costs, all_converged, total_time,
        Dict{Symbol, Any}(
            :n_steps_simulated => N_sim,
            :n_inner_converged => n_converged,
            :inner_solver_type => typeof(solver.inner_solver)
        )
    )
end
