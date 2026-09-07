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
    _shift_trajectories(trajs, game, dt) -> Vector{Trajectory{T}}

Shift each player's trajectory forward by one step for warm-starting, preserving the
node count exactly:
  states  (n_x × N+1): columns 2…N+1, then a GENUINE dynamics rollout from the new
                       terminal state/control pair (session correctness fix, see below).
  controls(n_u × N)  : columns 2…N,   then repeat the terminal control (the classical
                       Mattingley/Boyd MPC convention: hold the last control, per the file
                       header above -- this part was always correct).

Session correctness fix: the previous version repeated the OLD terminal STATE verbatim as
the new terminal state (`new_states = hcat(states[:,2:N+1], states[:,N+1:N+1])`), paired
with the repeated-hold control. This does NOT satisfy `x_new_{N+1} = f_d(x_new_N,
u_new_N)` in general -- it would require the trajectory to already be at a fixed point
under the held control, essentially never true for a real nonlinear system. Concretely:
new_states_prefix = [x₂,…,x_{N+1}] (shifted), new_controls = [û₂,…,û_N,û_N] (held); the OLD
code then appended x_{N+1} again as column N+2, but the dynamics-consistent value is
`f_d(x_{N+1}, û_N)`, not `x_{N+1}` itself. Confirmed on a real receding-horizon window
(solar-blocking DTMPC, window 97): this defect alone accounted for the ENTIRE dynamics-
infeasibility gap at the seed (χ(z̄) dropped from 0.584 to 4.5e-6 after correcting just
this one node), though it was not sufficient by itself to explain that window's remaining
failure (a separate hard-constraint issue) -- see the session's window-97 diagnostic.

`game` is the SAME (possibly tube/uncertainty-AUGMENTED) `GameProblem` the trajectories
were solved against. The genuine dynamics-rollout tail is only computable when `game.dynamics
isa SeparableDynamics`, i.e. `game.dynamics.player_dynamics[i]` gives a per-player function
fᵢ(xᵢ, uᵢ, p, t) operating on THAT PLAYER'S OWN state/control slice alone (physical rows AND
any augmented tube-channel rows alike -- no special-casing needed for augmented DTMPC
states, since the augmented game's own dynamics already propagate them correctly). This is
the actual dynamics representation used by every production RHN/DTMPC caller in this
repo (built via `PDGNEProblem`/`IGNEP.jl`). `dt` is the game's own timestep
(`game.time_horizon.dt`).

For any OTHER `DynamicsSpec` (e.g. `LinearDynamics`, which shares ONE joint state across
all players rather than stacking a private slice per player -- `x(k+1) = A x(k) + Σᵢ
Bᵢuᵢ(k)` genuinely needs every player's control, not just this trajectory's own) this
function falls back to the ORIGINAL naive terminal-state repeat. Do not attempt to unify
this without also redesigning the per-player `map` below to see every player's
trajectory jointly -- that redesign was explicitly out of scope for the session that
fixed the `SeparableDynamics` case (confirmed causally, see below) and this fallback is
what keeps the FNELQ+RHN (`LinearDynamics`) test suite passing unchanged.

Length preservation is REQUIRED, not cosmetic: FALCON's `_extract_flat_trajectory`
indexes `U[:, 1:N]` / `X[:, 1:N+1]` with no resampling and throws a `BoundsError` on a
short seed. `times` are window-local and therefore unchanged, matching the convention in
`_shift_feedback_strategy`/`_shift_open_loop_strategy` above.

The stale first column is harmless: the consumer overwrites node 1 with the true measured
initial state before solving.
"""
function _shift_trajectories(trajs::Vector{Trajectory{T}}, game::GameProblem{T}, dt::T) where {T}
    dyn = game.dynamics
    map(trajs) do tr
        N = size(tr.controls, 2)
        if N > 1
            new_states_prefix = tr.states[:, 2:N+1]
            new_controls      = hcat(tr.controls[:, 2:N], tr.controls[:, N:N])
            if dyn isa SeparableDynamics
                i = findfirst(obj -> obj.player_id == tr.player_id, game.objectives)
                fᵢ = dyn.player_dynamics[i]
                x_last = new_states_prefix[:, N]
                u_last = new_controls[:, N]
                x_new_tail = fᵢ(x_last, u_last, nothing, (N - 1) * dt)
                new_states = hcat(new_states_prefix, x_new_tail)
            else
                # LinearDynamics (shared joint state) and any other non-separable
                # DynamicsSpec: per-player rollout is not well-defined from a single
                # player's own Trajectory alone. Preserve the previously-tested behavior.
                new_states = hcat(new_states_prefix, new_states_prefix[:, N:N])
            end
        else
            new_states   = copy(tr.states)      # N=1: no shift possible
            new_controls = copy(tr.controls)
        end
        Trajectory{T}(tr.player_id, new_states, new_controls, copy(tr.times), tr.cost)
    end
end

"""
    _shift_warmstart(sol::GNEPSolution{T}) -> WarmstartData{T}

Build a warm-start for the next sub-problem from a converged solution.

The previous solution's TRAJECTORIES and strategy are both shifted forward by one step.
Dual variables are forwarded so constrained solvers (ALGAMES, FALCON) can continue from
their previous multiplier estimates.

NOTE (bug fix): this previously passed `nothing` for `trajectories` and looked duals up
under `:dual_variables` — a key FALCON does not emit (it emits `:λ` / `:dual_variables_λ`).
The result was that `warm_start=true` was a silent no-op: every receding-horizon window
fell through to FALCON's zero-control cold rollout (`_init_flat_trajectory`) with no dual
warm start, which is the entire point of a receding-horizon loop. Consumers all guard with
`isnothing`/length checks (see FALCON's `_try_load_dual!`), so supplying real data here is
backward compatible.
"""
function _shift_warmstart(sol::GNEPSolution{T}) where {T}
    shifted_strat = if sol.strategy isa FeedbackStrategy
        _shift_feedback_strategy(sol.strategy)
    elseif sol.strategy isa OpenLoopStrategy
        _shift_open_loop_strategy(sol.strategy)
    else
        nothing
    end

    shifted_trajs = isnothing(sol.trajectories) ? nothing :
        _shift_trajectories(sol.trajectories, sol.game, sol.game.time_horizon.dt)

    # Carry dual variables. Try the generic key first, then the keys FALCON actually
    # emits; wrap in a Dict{Symbol,Any} so downstream `haskey` lookups work uniformly.
    duals = get(sol.solver_info, :dual_variables, nothing)
    if isnothing(duals)
        d = Dict{Symbol, Any}()
        for k in (:λ, :dual_variables_λ, :μ)
            haskey(sol.solver_info, k) && (d[k] = sol.solver_info[k])
        end
        duals = isempty(d) ? nothing : d
    end

    return WarmstartData{T}(shifted_trajs, shifted_strat, duals, Dict{Symbol, Any}())
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
