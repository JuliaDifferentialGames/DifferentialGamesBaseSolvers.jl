# ============================================================================
# algames.jl
#
# ALGAMES: Augmented Lagrangian GAMEtheoretic Solver
#
# Reference: Le Cleac'h, Schwager, Manchester (2021),
#   "ALGAMES: A Fast Augmented Lagrangian Solver for Constrained Dynamic Games"
#   https://arxiv.org/abs/2104.08452
#
# ── Solution concept ─────────────────────────────────────────────────────────
#   Open-loop Nash equilibrium. With shared inequality constraints the solver
#   converges to a Normalized Nash Equilibrium (NNE, Rosen 1967) because
#   identical dual-ascent updates enforce equal multipliers on shared constraints
#   (Eqs. 20–23 and Section 6.4 of the paper).
#
# ── Algorithm (Algorithm 3) ──────────────────────────────────────────────────
#   Outer loop (up to outer_iter):
#     1. Newton inner loop: solve G(y) = 0 via δy = −H⁻¹G + backtracking.
#     2. Convergence check on ‖G_stat‖, ‖D‖, max(C).
#     3. Dual ascent on λ (Eq. 10).
#     4. Penalty increase ρ ← min(γρ, ρ_max) (Eq. 11).
#
#   G = [stationarity conditions for all players; dynamics residual]
#   H = ∂G/∂y computed via ForwardDiff (quasi-Newton, Section 4.5)
#
# ── Multiplier taxonomy ──────────────────────────────────────────────────────
#   μᵛₖ = Λ_dyn[i][:, k]   per-player dynamics dual    (solved in Newton)
#   λ                        shared constraint dual      (dual ascent, outer)
#   ρ                        shared penalty weight       (geometric schedule)
#
# ── File layout ──────────────────────────────────────────────────────────────
#   algames.jl  — ALGAMES struct, solver_capabilities, _solve, workspace
#   utils.jl    — residual/Jacobian assembly, line search, index helpers
# ============================================================================

# ── Namespace note ────────────────────────────────────────────────────────────
# This file is included by DifferentialGamesBaseSolvers. All names it uses
# (LinearAlgebra, ForwardDiff, Printf, DifferentialGamesBase types) must be
# imported at MODULE level in DifferentialGamesBaseSolvers.jl.
#
# Critically, `_solve` and `solver_capabilities` must be *imported* (not merely
# `using`-d) at module level so that definitions here add METHODS to DGB's
# existing generics rather than creating new module-local functions:
#
#   import DifferentialGamesBase: _solve, solver_capabilities
#
# This is already present in DifferentialGamesBaseSolvers.jl.  Do not add
# `using` or `import` statements to this file.

# ============================================================================
# ALGAMES struct
# ============================================================================

"""
    ALGAMES <: GameSolver

Augmented Lagrangian solver for open-loop (Generalized) Nash equilibria.

Default hyperparameters match the paper's recommended settings (ρ_init=1,
ρ_increase=10, Section 5.4).

# Keyword arguments
## Outer AL loop
- `outer_iter`   : Maximum AL iterations                          (50)
- `ρ_init`       : Initial penalty weight ρ⁽⁰⁾                   (1.0)
- `ρ_increase`   : Geometric multiplier γ > 1 (Eq. 11)            (10.0)
- `ρ_max`        : Cap on ρ                                       (1e8)

## Newton inner loop
- `inner_iter`   : Maximum Newton steps per outer iteration       (10)
- `reg`          : Tikhonov regularisation on H (ensures invertibility
                   when Nash equilibrium is non-unique, Section 6.1) (1e-6)

## Regularisation and line search
- `ls_iter`      : Maximum backtracks in Armijo line search           (20)
- `ls_beta`      : Sufficient-decrease fraction β                     (0.1)
- `ls_tau`       : Step contraction τ                                 (0.5)
- `reg`          : Base regularisation reg_0; actual reg per inner    (1e-3)
                   step l is reg_0 * l^4 (matches Algames.jl schedule)

## Convergence tolerances
- `tol_opt`      : ‖G_stationarity‖₁ / len                      (1e-4)
- `tol_dyn`      : ‖D‖₁ / len                                   (1e-4)
- `tol_con`      : max constraint violation                       (1e-3)

## Warm-start
- `reset_duals`  : If false, load λ, ρ, μ from WarmstartData     (true)
"""
@kwdef struct ALGAMES <: GameSolver
    outer_iter  ::Int     = 50
    ρ_init      ::Float64 = 1.0
    ρ_increase  ::Float64 = 10.0
    ρ_max       ::Float64 = 1e8
    inner_iter  ::Int     = 10
    reg         ::Float64 = 1e-3
    ls_iter     ::Int     = 20
    ls_beta     ::Float64 = 0.1
    ls_tau      ::Float64 = 0.5
    tol_opt     ::Float64 = 1e-4
    tol_dyn     ::Float64 = 1e-4
    tol_con     ::Float64 = 1e-3
    reset_duals ::Bool    = true
end

function solver_capabilities(::Type{ALGAMES})
    return [
        :LQGame,
        :NonlinearGame,
        :GNEP,
        :ConstrainedGame,
        :UnconstrainedGame,
        :DiscreteTime,
        :SeparableDynamics,
        :SharedConstraints,
    ]
end

# ============================================================================
# ALGAMESWorkspace
# ============================================================================

"""
    ALGAMESWorkspace{T}

Mutable primal-dual workspace for the Newton iterate.

## Primal
- `X`    : (n × N+1)    shared state; `X[:,1] = x₀` is fixed throughout
- `U`    : [mᵢ × N, …]  per-player controls

## Dynamics duals (μᵛ in the paper)
- `Λ_dyn` : [(n × N), …]  `Λ_dyn[i][:,k]` = μᵢₖ ∈ ℝⁿ

## Constraint duals (λ, ρ in the paper)
- `λ`    : (nc_step × N,) flattened; `λ[(k-1)*nc_step+j]` = λⱼ at step k
- `ρ`    : same layout as λ

## Cached dimensions
- `n`, `np`, `N`, `nc_step` : problem sizes
- `control_dims`    : [m₁, …, mₙₚ]
- `control_offsets` : zero-based start of each player's slice in joint u
"""
mutable struct ALGAMESWorkspace{T}
    X               ::Matrix{T}
    U               ::Vector{Matrix{T}}
    Λ_dyn           ::Vector{Matrix{T}}
    λ               ::Vector{T}
    ρ               ::Vector{T}
    n               ::Int
    np              ::Int
    N               ::Int
    nc_step         ::Int
    control_dims    ::Vector{Int}
    control_offsets ::Vector{Int}
end

# ============================================================================
# Workspace constructor
# ============================================================================

function ALGAMESWorkspace(
    game     ::GameProblem{T},
    solver   ::ALGAMES,
    warmstart::Union{Nothing, WarmstartData}
) where {T}
    n    = total_state_dim(game.dynamics)
    np   = game.n_players
    N    = n_steps(game)
    cdim = game.metadata.control_dims
    coff = game.metadata.control_offsets

    # ── Primal initialisation ────────────────────────────────────────────────
    # Default: zero controls, roll out dynamics from x₀ (Section 4.6 default).
    U = [zeros(T, cdim[i], N) for i in 1:np]

    if warmstart !== nothing && warmstart.trajectories !== nothing
        for i in 1:np
            idx = findfirst(t -> t.player_id == i, warmstart.trajectories)
            idx === nothing && continue
            tr  = warmstart.trajectories[idx]
            nu  = min(size(tr.controls, 2), N)
            U[i][:, 1:nu] .= tr.controls[:, 1:nu]
        end
    end

    X = _rollout_flat(game.dynamics, game.initial_state, U, N, game.time_horizon.dt)

    # ── Constraint dimension ─────────────────────────────────────────────────
    nc_step   = _count_nc_step(game)
    n_c_total = nc_step * N

    # ── Dual initialisation ──────────────────────────────────────────────────
    Λ_dyn = [zeros(T, n, N) for _ in 1:np]
    λ      = zeros(T, n_c_total)
    ρ      = fill(T(solver.ρ_init), n_c_total)

    if !solver.reset_duals && warmstart !== nothing
        dv = warmstart.dual_variables
        if dv !== nothing
            if haskey(dv, :λ) && length(dv[:λ]) == n_c_total
                λ .= convert(Vector{T}, dv[:λ])
            end
            if haskey(dv, :ρ) && length(dv[:ρ]) == n_c_total
                ρ .= convert(Vector{T}, dv[:ρ])
            end
            if haskey(dv, :Λ_dyn)
                prev = dv[:Λ_dyn]
                for i in 1:min(length(prev), np)
                    sz = min(size(prev[i], 2), N)
                    Λ_dyn[i][:, 1:sz] .= prev[i][:, 1:sz]
                end
            end
        end
    end

    return ALGAMESWorkspace{T}(
        X, U, Λ_dyn, λ, ρ, n, np, N, nc_step, cdim, coff
    )
end

# ============================================================================
# _solve  —  outer AL loop  (Algorithm 3)
# ============================================================================

# Concrete overload for GNEPSolution warmstart.
# DGB's forwarding method (_solve(::GameProblem, ::GameSolver, ::GNEPSolution, ::Bool))
# re-dispatches with a ::GameSolver-typed solver, which misses our concrete method.
# This overload intercepts the GNEPSolution case at the concrete ::ALGAMES level
# and converts to WarmstartData before calling the main method — matching the
# pattern used by iLQGames and other solvers in this package.
function _solve(
    game     ::GameProblem{T},
    solver   ::ALGAMES,
    warmstart::GNEPSolution{T},
    verbose  ::Bool
) where {T}
    _solve(game, solver, WarmstartData(warmstart), verbose)
end

"""
    _solve(game, solver::ALGAMES, warmstart, verbose) -> GNEPSolution{T}

Main entry point.  See module docstring for algorithm description.
"""
function _solve(
    game     ::GameProblem{T},
    solver   ::ALGAMES,
    warmstart::Union{Nothing, WarmstartData},
    verbose  ::Bool
) where {T}
    t_start = time()
    ws  = ALGAMESWorkspace(game, solver, warmstart)

    # G and y have the same length (see _G_dim == _y_dim in utils.jl).
    n_y = _y_dim(ws)
    buf = ALGAMESBuffers(T, n_y)

    converged = false
    outer_k   = 0
    opt_vio   = T(Inf)
    dyn_vio   = T(Inf)
    con_vio   = T(Inf)

    # Per-outer-iteration convergence history (for plotting)
    history_opt = T[]
    history_dyn = T[]
    history_con = T[]

    for k in 1:solver.outer_iter
        outer_k = k

        # ── Newton inner loop ─────────────────────────────────────────────
        # Matches reference Algames.jl solver_methods.jl:
        #   reg = reg_0 * l^4 (scheduled, reset each outer iteration)
        #   Line search with ρ_trial (fixed) — see _line_search
        #   Break on: (a) stationarity met, (b) Δ_step < Δ_min, (c) LS exhausted
        _build_residual!(buf, game, ws)
        Δ_min = T(1e-9)

        for l in 1:solver.inner_iter
            G_norm = norm(buf.G, 1)

            _build_jacobian!(buf, game, ws)

            reg_l = T(solver.reg) * T(l)^4
            _regularize!(buf.H, ws, Float64(reg_l))

            δy = _newton_step(buf.H, buf.G)
            α  = _line_search(game, ws, buf, G_norm, δy, solver)
            _apply_step!(ws, δy, T(α))

            _build_residual!(buf, game, ws)
            G_norm_new = norm(buf.G, 1)
            Δ_step     = α * norm(δy, 1) / n_y

            verbose && @printf("  [inner %2d] ‖G‖₁/n=%.3e  α=%.4f  reg=%.1e\n",
                               l, G_norm_new / n_y, α, reg_l)

            G_norm_new / n_y < solver.tol_opt && break   # (a) converged
            Δ_step < Δ_min && break                       # (b) step too small
            α ≤ solver.ls_tau^(solver.ls_iter - 1) && break  # (c) LS exhausted
        end

        # ── Convergence check ─────────────────────────────────────────────
        # buf.G is current — the inner loop ends with _build_residual!
        opt_vio = _stationarity_norm(buf, ws)
        dyn_vio = _dynamics_norm(buf, ws)
        con_vio = _constraint_violation(game, ws)

        verbose && @printf(
            "[outer %3d] opt=%.3e  dyn=%.3e  con=%.3e\n",
            k, opt_vio, dyn_vio, con_vio
        )

        push!(history_opt, opt_vio)
        push!(history_dyn, dyn_vio)
        push!(history_con, con_vio)

        if opt_vio < solver.tol_opt &&
           dyn_vio < solver.tol_dyn &&
           con_vio < solver.tol_con
            converged = true
            break
        end

        # ── Dual ascent (Eq. 10) ──────────────────────────────────────────
        _dual_ascent!(ws, game)

        # ── Penalty schedule (Eq. 11) ─────────────────────────────────────
        ws.ρ .= min.(ws.ρ .* T(solver.ρ_increase), T(solver.ρ_max))
    end

    return _assemble_solution(
        game, ws, converged, outer_k,
        time() - t_start, opt_vio, dyn_vio, con_vio,
        history_opt, history_dyn, history_con
    )
end

# ============================================================================
# _apply_step!
# ============================================================================

"""
    _apply_step!(ws, δy, α)

Update y ← y + α·δy by unpacking δy into (X, U, Λ_dyn).
Packing matches `_pack_y` in utils.jl:
  y = [x(2)…x(N+1) | U[1] | … | U[np] | Λ_dyn[1] | … | Λ_dyn[np]]
x(1) = x₀ is NEVER updated.
"""
function _apply_step!(ws::ALGAMESWorkspace{T}, δy::AbstractVector{T}, α::T) where {T}
    n, np, N = ws.n, ws.np, ws.N
    off = 0

    for k in 1:N
        ws.X[:, k+1] .+= α .* δy[off+1:off+n]
        off += n
    end
    for i in 1:np
        mi  = ws.control_dims[i]
        len = mi * N
        ws.U[i] .+= α .* reshape(δy[off+1:off+len], mi, N)
        off += len
    end
    for i in 1:np
        len = n * N
        ws.Λ_dyn[i] .+= α .* reshape(δy[off+1:off+len], n, N)
        off += len
    end

    @assert off == length(δy) "y packing mismatch: off=$off ≠ $(length(δy))"
end

# ============================================================================
# _dual_ascent!  (Eq. 10)
# ============================================================================

"""
    _dual_ascent!(ws, game)

Dual ascent on shared constraint multipliers λ (Eq. 10):
  Inequality j:  λ_j ← max(0, λ_j + ρ_j · C_j(X_k, U_k))
  Equality j:    λ_j ← λ_j + ρ_j · C_j(X_k, U_k)

λ and ρ are stored with per-timestep blocks of length nc_step.
This ensures the active-set logic in `_al_pen` can read the correct λ
for each timestep independently.

NNE property: because λ is initialised to zero and both players of a
shared constraint receive the SAME update, by induction λ^ν_k = λ^ω_k
for all shared constraints (Eqs. 20–23), yielding a Normalized Nash
Equilibrium (Section 6.4).
"""
function _dual_ascent!(ws::ALGAMESWorkspace{T}, game::GameProblem{T}) where {T}
    ws.nc_step == 0 && return
    U_flat = _concat_controls(ws)

    for k in 1:ws.N
        x_k = ws.X[:, k]
        u_k = U_flat[:, k]
        λ_k = view(ws.λ, (k-1)*ws.nc_step+1 : k*ws.nc_step)
        ρ_k = view(ws.ρ, (k-1)*ws.nc_step+1 : k*ws.nc_step)
        off = 0

        for c in Iterators.flatten((game.private_constraints,
                                    game.shared_constraints))
            cv = evaluate_constraint(c, x_k, u_k, nothing, k)
            nc = length(cv)
            for j in 1:nc
                g     = off + j
                new_λ = λ_k[g] + ρ_k[g] * cv[j]
                λ_k[g] = is_inequality(c) ? max(zero(T), new_λ) : new_λ
            end
            off += nc
        end

        @assert off == ws.nc_step "Constraint dim mismatch at k=$k"
    end
end

# ============================================================================
# _assemble_solution
# ============================================================================

function _assemble_solution(
    game        ::GameProblem{T},
    ws          ::ALGAMESWorkspace{T},
    converged   ::Bool,
    iters       ::Int,
    wall_time   ::Float64,
    opt_vio     ::T,
    dyn_vio     ::T,
    con_vio     ::T,
    history_opt ::Vector{T},
    history_dyn ::Vector{T},
    history_con ::Vector{T},
) where {T}
    th    = game.time_horizon
    N     = ws.N
    times = collect(range(T(0), th.tf, length = N + 1))

    U_flat = vcat(ws.U...)   # full joint control (m_total × N)

    trajectories = map(1:ws.np) do i
        obj = get_objective(game, i)

        # State list: LQStageCost with separable dynamics needs private slice (no offset fields).
        # NonlinearStageCost (CompositeCostTerm) uses player_slice internally → full joint state.
        if obj.stage_cost isa LQStageCost && length(game.metadata.state_dims) == ws.np
            sdim_i = game.metadata.state_dims[i]
            soff_i = game.metadata.state_offsets[i]
            srng_i = soff_i+1 : soff_i+sdim_i
            X_list = [ws.X[srng_i, k] for k in 1:N+1]
        else
            X_list = [ws.X[:, k] for k in 1:N+1]
        end

        # Control list: LQStageCost needs private control slice (R is mᵢ×mᵢ, no offset).
        # NonlinearStageCost (QuadraticControlCost) uses player_slice internally → full joint control.
        if obj.stage_cost isa LQStageCost
            U_list = [ws.U[i][:, k] for k in 1:N]
        else
            U_list = [U_flat[:, k] for k in 1:N]
        end

        c_i = total_cost(obj, X_list, U_list, nothing)
        Trajectory{T}(i, copy(ws.X), copy(ws.U[i]), times, c_i)
    end

    strategy = OpenLoopStrategy(
        [copy(ws.U[i]) for i in 1:ws.np],
        copy(ws.control_dims),
        times
    )

    solver_info = Dict{Symbol, Any}(
        :opt_vio => opt_vio,
        :dyn_vio => dyn_vio,
        :con_vio => con_vio,
        :history => Dict{Symbol, Any}(
            :opt => copy(history_opt),
            :dyn => copy(history_dyn),
            :con => copy(history_con),
        ),
        :dual_variables => Dict{Symbol, Any}(
            :λ     => copy(ws.λ),
            :ρ     => copy(ws.ρ),
            :Λ_dyn => [copy(L) for L in ws.Λ_dyn],
        )
    )

    return GNEPSolution(
        game, trajectories;
        state_trajectory = copy(ws.X),
        strategy         = strategy,
        equilibrium_type = :OpenLoopNash,
        converged        = converged,
        iterations       = iters,
        solve_time       = wall_time,
        solver_info      = solver_info,
    )
end