# ============================================================================
# libr.jl — Lexicographic Iterative Best Response (L-IBR)
#
# Reference: Miller, K. & Mitra, S. (2022).
#   "Multi-agent motion planning using differential games with lexicographic
#    preferences." IEEE CDC 2022, pp. 5751–5757.
#
# ── Game class ───────────────────────────────────────────────────────────────
#   Operates on LexicographicGameProblem{T}, which wraps a PDGNEProblem
#   (SeparableDynamics) and adds per-player lexicographic cost structure:
#     Jᵢ(z) = (Jᵢᶜᵒˡ, Jᵢᵖᵉʳ)  ∈ ℝ²  ordered lexicographically.
#
# ── Algorithm (Algorithm 1 of Miller & Mitra 2022) ───────────────────────────
#   Initialize: U[i] = 0 ∈ ℝ^{mᵢ × N}  for all i; rollout to get X[i].
#
#   for iter = 1, 2, ..., max_iter:
#     prev_U = copy(U)
#     for i = 1, ..., N:
#       z_{-i} = {(X[j], U[j]) : j ≠ i}            # fixed opponent trajectories
#
#       # Phase 1 — minimize collision cost Jᵢᶜᵒˡ
#       U[i] ← gradient_descent(Jᵢᶜᵒˡ(·, z_{-i}), U[i])
#       J*_col = Jᵢᶜᵒˡ(z_i, z_{-i})                # record Phase 1 optimum
#       X[i] ← rollout(fᵢ, x_i0, U[i])
#
#       # Phase 2 — minimize personal cost Jᵢᵖᵉʳ subject to Jᵢᶜᵒˡ ≤ J*_col + slack
#       U[i] ← gradient_descent(Jᵢᵖᵉʳ + ρ·max(0, Jᵢᶜᵒˡ−J*_col)², U[i])
#       X[i] ← rollout(fᵢ, x_i0, U[i])
#
#     converged? max_i ‖U[i] − prev_U[i]‖_∞ < tol → break
#
# ── Trajectory representation ────────────────────────────────────────────────
#   z_i is a Tuple{Matrix, Matrix} = (X_i, U_i) where:
#     X_i ∈ ℝ^{nᵢ × (N+1)} — private state trajectory (from Euler rollout)
#     U_i ∈ ℝ^{mᵢ × N}    — private control trajectory (optimization variable)
#   User cost functions must accept this representation.
#
# ── Integration ──────────────────────────────────────────────────────────────
#   Euler: x_{k+1} = x_k + dt · fᵢ(x_k, u_k, nothing, (k−1)·dt).
#   ForwardDiff-compatible — dual numbers propagate through rollout and costs.
#
# ── Namespace note ────────────────────────────────────────────────────────────
#   Included by DifferentialGamesBaseSolvers. All DGB types are available via
#   `using DifferentialGamesBase` in the parent module.
# ============================================================================

# ============================================================================
# LIBR solver struct
# ============================================================================

"""
    LIBR <: GameSolver

Lexicographic Iterative Best Response solver for `LexicographicGameProblem{T}`.

Implements Algorithm 1 of Miller & Mitra (2022). Each outer IBR iteration
updates every agent sequentially with two gradient-descent phases:

- **Phase 1**: minimize collision cost Jᵢᶜᵒˡ with all others fixed; record J*_col.
- **Phase 2**: minimize personal cost Jᵢᵖᵉʳ with a quadratic penalty enforcing
  Jᵢᶜᵒˡ ≤ J*_col + `col_slack`.

Converges to a Lexicographic Nash Equilibrium (L-NE), guaranteed to exist
as a pure strategy by Proposition 2 of Miller & Mitra (2022).

# Trajectory convention
User-supplied cost functions receive `z_i = (X_i, U_i)`:
  - `X_i :: Matrix{T}` of shape `(nᵢ, N+1)` — state trajectory
  - `U_i :: Matrix{T}` of shape `(mᵢ, N)`   — control trajectory

# Keyword Arguments
## IBR loop
- `max_iter`    : Maximum outer IBR iterations                         (100)
- `tol`         : Convergence: max change in any player's control       (1e-4)

## Per-phase gradient descent
- `inner_iter`  : Gradient steps per phase per player per IBR iteration (200)
- `step_size`   : Initial step size α₀ for backtracking line search     (0.01)
- `ls_beta`     : Step contraction factor for line search               (0.5)
- `ls_max_iter` : Maximum backtracking steps                            (15)
- `ls_armijo`   : Armijo sufficient-decrease constant c                 (1e-4)

## Phase 2 penalty
- `col_penalty` : Quadratic penalty weight ρ for Jᵢᶜᵒˡ > J*_col        (1e3)
- `col_slack`   : Absolute slack: Phase 2 target is Jᵢᶜᵒˡ ≤ J*_col + δ (1e-6)
"""
@kwdef struct LIBR <: GameSolver
    max_iter    ::Int     = 100
    tol         ::Float64 = 1e-4
    inner_iter  ::Int     = 200
    step_size   ::Float64 = 0.01
    ls_beta     ::Float64 = 0.5
    ls_max_iter ::Int     = 15
    ls_armijo   ::Float64 = 1e-4
    col_penalty ::Float64 = 1e3
    col_slack   ::Float64 = 1e-6
end

function solver_capabilities(::Type{LIBR})
    return [
        :LexicographicGame,
        :OpenLoopPolicies,
        :DiscreteTime,
        :SeparableDynamics,
    ]
end

# ============================================================================
# _libr_euler_rollout — ForwardDiff-compatible single-player Euler rollout
# ============================================================================

"""
    _libr_euler_rollout(f_i, x_i0, U_i, dt) -> Matrix

Euler-integrate player i's private continuous-time dynamics over N steps:
  x_{k+1} = x_k + dt · fᵢ(x_k, u_k, nothing, (k−1)·dt)

Accepts `AbstractVector` and `AbstractMatrix` for dual-number propagation
through ForwardDiff gradient computation.
"""
function _libr_euler_rollout(f_i, x_i0::AbstractVector{S}, U_i::AbstractMatrix, dt::Real) where {S}
    n_i = length(x_i0)
    N   = size(U_i, 2)
    TE  = promote_type(S, eltype(U_i))
    X_i = Matrix{TE}(undef, n_i, N + 1)
    X_i[:, 1] .= x_i0
    for k in 1:N
        ẋ = f_i(X_i[:, k], U_i[:, k], nothing, (k - 1) * dt)
        X_i[:, k+1] .= X_i[:, k] .+ dt .* ẋ
    end
    return X_i
end

# ============================================================================
# _libr_gradient_step — one gradient step with Armijo backtracking
# ============================================================================

"""
    _libr_gradient_step(cost_fn, u_vec, solver) -> (u_new, step_taken)

Compute the gradient of `cost_fn` at `u_vec` via ForwardDiff, then take an
Armijo-backtracked gradient step. Returns the new iterate and the step size.
"""
function _libr_gradient_step(cost_fn, u_vec::AbstractVector{T}, solver::LIBR) where {T}
    f0   = cost_fn(u_vec)
    grad = ForwardDiff.gradient(cost_fn, u_vec)
    g2   = dot(grad, grad)          # ‖∇f‖²

    iszero(g2) && return (u_vec, zero(T))

    α = T(solver.step_size)
    for _ in 1:solver.ls_max_iter
        u_trial = u_vec .- α .* grad
        if cost_fn(u_trial) ≤ f0 - T(solver.ls_armijo) * α * g2
            return (u_trial, α)
        end
        α *= T(solver.ls_beta)
    end
    # Accept smallest step even if Armijo not met — avoids stall
    return (u_vec .- α .* grad, α)
end

# ============================================================================
# _libr_phase — inner gradient descent for one phase
# ============================================================================

"""
    _libr_phase(game, i, f_i, x_i0, U_i_init, z_minus_i, dt, solver;
                phase=:collision, J_col_star=nothing)
    -> (U_i_opt, J_col_opt)

Run gradient descent for player i's Phase 1 (`:collision`) or Phase 2
(`:personal`) optimization.

**Phase 1** minimizes Jᵢᶜᵒˡ(z_i, z_{−i}) over U_i.
**Phase 2** minimizes Jᵢᵖᵉʳ(z_i, z_{−i}) + ρ·max(0, Jᵢᶜᵒˡ − J*_col + δ)²
  where J*_col = `J_col_star` and δ = `solver.col_slack`.

Returns `(U_i_opt, J_col_opt)`:
  - `U_i_opt`   : Optimized control trajectory for player i
  - `J_col_opt` : Collision cost at `U_i_opt` (same units, useful for bookkeeping)
"""
function _libr_phase(
    game       ::LexicographicGameProblem{T},
    i          ::Int,
    f_i        ,
    x_i0       ::AbstractVector{T},
    U_i_init   ::Matrix{T},
    z_minus_i  ::AbstractVector,
    dt         ::Real,
    solver     ::LIBR;
    phase      ::Symbol = :collision,
    J_col_star ::Union{Nothing, T} = nothing
) where {T}
    m_i = size(U_i_init, 1)
    N   = size(U_i_init, 2)
    ρ   = T(solver.col_penalty)
    δ   = T(solver.col_slack)

    u_vec = vec(copy(U_i_init))

    for _ in 1:solver.inner_iter
        cost_fn = if phase === :collision
            uv -> begin
                U_ad = reshape(uv, m_i, N)
                X_ad = _libr_euler_rollout(f_i, x_i0, U_ad, dt)
                collision_cost(game, i, (X_ad, U_ad), z_minus_i)
            end
        else # :personal
            uv -> begin
                U_ad  = reshape(uv, m_i, N)
                X_ad  = _libr_euler_rollout(f_i, x_i0, U_ad, dt)
                j_col = collision_cost(game, i, (X_ad, U_ad), z_minus_i)
                j_per = personal_cost(game, i, (X_ad, U_ad), z_minus_i)
                viol  = max(zero(eltype(uv)), j_col - (J_col_star + δ))
                j_per + ρ * viol^2
            end
        end

        u_new, _ = _libr_gradient_step(cost_fn, u_vec, solver)
        u_vec = u_new
    end

    U_opt   = reshape(u_vec, m_i, N)
    X_opt   = _libr_euler_rollout(f_i, x_i0, U_opt, dt)
    J_col   = collision_cost(game, i, (X_opt, U_opt), z_minus_i)

    return U_opt, T(J_col)
end

# ============================================================================
# _solve — main L-IBR loop
# ============================================================================

# GNEPSolution warmstart overload (mirrors ALGAMES pattern)
function _solve(
    game     ::LexicographicGameProblem{T},
    solver   ::LIBR,
    warmstart::GNEPSolution{T},
    verbose  ::Bool
) where {T}
    _solve(game, solver, WarmstartData(warmstart), verbose)
end

"""
    _solve(game::LexicographicGameProblem{T}, solver::LIBR, warmstart, verbose)
    -> GNEPSolution{T}

Main entry point for L-IBR. See module docstring for algorithm description.
"""
function _solve(
    game     ::LexicographicGameProblem{T},
    solver   ::LIBR,
    warmstart::Union{Nothing, WarmstartData},
    verbose  ::Bool
) where {T}
    t_start  = time()
    fwd      = game.forward_game
    N        = n_steps(fwd)
    dt       = fwd.time_horizon.dt
    np       = game.n_players
    dyn      = fwd.dynamics   # SeparableDynamics

    sdims  = dyn.state_dims
    soffs  = fwd.metadata.state_offsets
    cdims  = fwd.metadata.control_dims
    coffs  = fwd.metadata.control_offsets
    x0_all = fwd.initial_state   # joint initial state

    # ── Initialize controls ───────────────────────────────────────────────────
    U = Vector{Matrix{T}}(undef, np)
    for i in 1:np
        U[i] = zeros(T, cdims[i], N)
    end

    if warmstart !== nothing && warmstart.trajectories !== nothing
        for i in 1:np
            idx = findfirst(t -> t.player_id == i, warmstart.trajectories)
            idx === nothing && continue
            tr = warmstart.trajectories[idx]
            nc = min(size(tr.controls, 2), N)
            nr = min(size(tr.controls, 1), cdims[i])
            U[i][1:nr, 1:nc] .= tr.controls[1:nr, 1:nc]
        end
    end

    # ── Initial rollout for each player ──────────────────────────────────────
    X = Vector{Matrix{T}}(undef, np)
    for i in 1:np
        f_i   = dyn.player_dynamics[i]
        x_i0  = x0_all[soffs[i]+1 : soffs[i]+sdims[i]]
        X[i]  = _libr_euler_rollout(f_i, x_i0, U[i], dt)
    end

    converged  = false
    iters_done = 0
    Δ_hist     = T[]

    # ── Outer IBR loop ────────────────────────────────────────────────────────
    for iter in 1:solver.max_iter
        iters_done = iter
        U_prev     = [copy(Ui) for Ui in U]

        for i in 1:np
            f_i  = dyn.player_dynamics[i]
            x_i0 = x0_all[soffs[i]+1 : soffs[i]+sdims[i]]

            # Build z_{-i} from current (possibly updated) X and U
            z_minus_i = [(X[j], U[j]) for j in 1:np if j != i]

            # Phase 1: minimize collision cost
            U_new, J_col_star = _libr_phase(
                game, i, f_i, x_i0, U[i], z_minus_i, dt, solver;
                phase = :collision
            )
            U[i] = U_new
            X[i] = _libr_euler_rollout(f_i, x_i0, U[i], dt)

            # Refresh z_{-i} (X[i] changed, others unchanged within this sweep)
            z_minus_i = [(X[j], U[j]) for j in 1:np if j != i]

            # Phase 2: minimize personal cost s.t. J_col ≤ J_col_star + slack
            U_new2, _ = _libr_phase(
                game, i, f_i, x_i0, U[i], z_minus_i, dt, solver;
                phase       = :personal,
                J_col_star  = J_col_star
            )
            U[i] = U_new2
            X[i] = _libr_euler_rollout(f_i, x_i0, U[i], dt)
        end

        # ── Convergence check ─────────────────────────────────────────────────
        Δ = maximum(
            maximum(abs, U[i] .- U_prev[i]) for i in 1:np
        )
        push!(Δ_hist, Δ)

        verbose && @printf("[LIBR iter %3d] ΔU = %.3e\n", iter, Δ)

        if Δ < solver.tol
            converged = true
            break
        end
    end

    # ── Assemble solution ─────────────────────────────────────────────────────
    return _libr_assemble_solution(
        game, fwd, X, U, converged, iters_done, time() - t_start, Δ_hist
    )
end

# ============================================================================
# _libr_assemble_solution
# ============================================================================

function _libr_assemble_solution(
    game       ::LexicographicGameProblem{T},
    fwd        ::GameProblem{T},
    X          ::Vector{Matrix{T}},
    U          ::Vector{Matrix{T}},
    converged  ::Bool,
    iters      ::Int,
    wall_time  ::Float64,
    Δ_hist     ::Vector{T}
) where {T}
    np   = game.n_players
    N    = n_steps(fwd)
    th   = fwd.time_horizon
    times = collect(range(T(0), th.tf, length = N + 1))

    # Compute final lexicographic costs for all players
    J_col = Vector{T}(undef, np)
    J_per = Vector{T}(undef, np)
    for i in 1:np
        z_i  = (X[i], U[i])
        z_mi = [(X[j], U[j]) for j in 1:np if j != i]
        J_col[i] = collision_cost(game, i, z_i, z_mi)
        J_per[i] = personal_cost(game, i, z_i, z_mi)
    end

    trajectories = map(1:np) do i
        # Use J_col + J_per as scalar cost summary for the Trajectory record
        c_i = J_col[i] + J_per[i]
        Trajectory{T}(i, copy(X[i]), copy(U[i]), times, c_i)
    end

    strategy = OpenLoopStrategy(
        [copy(U[i]) for i in 1:np],
        fwd.metadata.control_dims,
        times
    )

    solver_info = Dict{Symbol, Any}(
        :J_col   => J_col,
        :J_per   => J_per,
        :Δ_hist  => Δ_hist,
    )

    return GNEPSolution(
        fwd, trajectories;
        strategy         = strategy,
        equilibrium_type = :Approximate,
        converged        = converged,
        iterations       = iters,
        solve_time       = wall_time,
        solver_info      = solver_info,
    )
end
