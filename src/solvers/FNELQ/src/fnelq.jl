"""
    FNELQ <: GameSolver

Solves N-player discrete-time LQ games for feedback Nash equilibrium via 
dynamic programming. Supports both LTI and LTV formulations via the unified
`LinearDynamics{T, AM, BM}` type — use `get_A(dyn, k)` / `get_B(dyn, i, k)`
accessors throughout; never access `.A` or `.B` directly in this file.

Based on Başar & Olsder (1999), Theorem 6.17 (discrete-time case).

# References
Başar, T., & Olsder, G. J. (1999). Dynamic Noncooperative Game Theory (2nd ed.). 
SIAM. Section 6.3, Theorem 6.17 (Equation 6.17a).
"""
struct FNELQ <: GameSolver
    check_singularity::Bool
    rcond_threshold::Float64

    function FNELQ(;
        check_singularity::Bool = true,
        rcond_threshold::Float64 = 1e-10
    )
        new(check_singularity, rcond_threshold)
    end
end

solver_capabilities(::Type{FNELQ}) = [
    :LQGame,
    :FeedbackPolicies,
    :UnconstrainedGame,
    :DiscreteTime
]

# ============================================================================
# extract_lq_matrices
# ============================================================================
#
# Returns the dynamics object and per-player cost accessors. Dispatches on
# dynamics type so callers get compile-time specialization.
#
# Returns:
#   dyn   — LinearDynamics (LTI or LTV); use get_A(dyn, k), get_B(dyn, i, k)
#   costs — Vector of PlayerObjective stage costs (LTI or LTV LQStageCost)
#           use get_Q(cost, k), get_R(cost, k), etc.
#   Qf    — Vector{Matrix{T}}, one per player (terminal cost, always static)
#   n     — state dimension (Int)
#   n_players — Int
#
# Validation of LQ structure is performed here so _solve can assume it.

"""
    extract_lq_matrices(game::GameProblem{T})

Extract dynamics and cost structure from an LQ `GameProblem`.

Returns `(dyn, stage_costs, Qf, n, n_players)` where:
- `dyn::LinearDynamics`          : use `get_A(dyn, k)`, `get_B(dyn, i, k)`
- `stage_costs::Vector`          : `LQStageCost` per player; use `get_Q(sc, k)` etc.
- `Qf::Vector{Matrix{T}}`        : terminal cost matrices (always static)
- `n::Int`                       : state dimension
- `n_players::Int`               : number of players

Dispatches on `LinearDynamics` subtype — LTI and LTV paths are resolved at
compile time via the `AM` type parameter; no runtime branching needed here.
"""
function extract_lq_matrices(game::GameProblem{T}) where {T}
    @assert game.dynamics isa LinearDynamics "extract_lq_matrices requires LinearDynamics"
    @assert all(obj.stage_cost isa LQStageCost for obj in game.objectives) "extract_lq_matrices requires LQStageCost objectives"

    dyn       = game.dynamics
    n         = dyn.state_dim
    n_players = game.n_players

    stage_costs = [game.objectives[i].stage_cost for i in 1:n_players]
    Qf          = [game.objectives[i].terminal_cost.Qf for i in 1:n_players]

    return dyn, stage_costs, Qf, n, n_players
end

# ============================================================================
# _solve
# ============================================================================

"""
    _solve(game::GameProblem{T}, solver::FNELQ, warmstart, verbose::Bool)

Solve discrete-time unconstrained LQ game for feedback Nash equilibrium.
Handles both LTI (`LinearDynamics{T, Matrix{T}, ...}`) and LTV
(`LinearDynamics{T, Vector{Matrix{T}}, ...}`) formulations transparently
via `get_A` / `get_B` / `get_Q` / `get_R` accessors.
"""
function _solve(
    game::GameProblem{T},
    solver::FNELQ,
    warmstart::Union{Nothing, WarmstartData},
    verbose::Bool
) where {T}

    @assert game.time_horizon isa DiscreteTime "FNELQ requires discrete-time formulation"
    @assert is_unconstrained(game) "FNELQ only handles unconstrained games"

    # ── Extract structure ────────────────────────────────────────────────────
    dyn, stage_costs, Qf, n, n_players = extract_lq_matrices(game)

    # Control dimensions — fixed across time even for LTV
    control_dims   = dyn.control_dims
    total_controls = sum(control_dims)
    control_ranges = [sum(control_dims[1:i-1])+1:sum(control_dims[1:i]) for i in 1:n_players]

    N  = n_steps(game)
    dt = game.time_horizon.dt

    # ── Cross-term / affine-term detection ───────────────────────────────────
    # For LTV costs, check the first timestep; if M/q/r are zero there and the
    # cost is LTI, that's representative. LTV detection is conservative: if any
    # player has a non-LTI cost we check all timesteps.
    function _has_nonzero_M(sc, k)
        M = get_M(sc, k)
        return norm(M) > eps(T)
    end
    function _has_affine(sc, k)
        norm(get_q(sc, k)) > eps(T) || norm(get_r(sc, k)) > eps(T)
    end

    # Check representative timesteps (k=1 for LTI, all k for LTV)
    ks_to_check = is_ltv(dyn) ? (1:N) : (1:1)

    if any(_has_nonzero_M(stage_costs[i], k) for i in 1:n_players, k in ks_to_check)
        @warn "FNELQ: Cross terms M ≠ 0 detected; not yet implemented in backward pass, ignoring"
    end
    has_affine_terms = any(_has_affine(stage_costs[i], k) for i in 1:n_players, k in ks_to_check)

    verbose && @info "FNELQ: Solving $(n_players)-player discrete LQ game" n_states=n time_steps=N dt=dt control_dims=control_dims is_ltv=is_ltv(dyn)

    t_start = time()

    # ── Initialise cost-to-go ────────────────────────────────────────────────
    Z = [copy(Qf[i]) for i in 1:n_players]   # Zᵢ(N) = Qfᵢ
    ζ = has_affine_terms ? [zeros(T, n) for _ in 1:n_players] : nothing

    # ── Storage ──────────────────────────────────────────────────────────────
    P_history = [Matrix{T}[] for _ in 1:n_players]
    α_history = has_affine_terms ? [Vector{T}[] for _ in 1:n_players] : nothing
    Z_history = [Matrix{T}[] for _ in 1:n_players]
    for i in 1:n_players
        push!(Z_history[i], copy(Z[i]))
    end
    rcond_history = Float64[]

    # ── Backward pass ────────────────────────────────────────────────────────
    for k in N:-1:1

        # Pull time-varying matrices at step k via accessors.
        # For LTI dynamics/costs the accessor ignores k (zero overhead).
        A_k = get_A(dyn, k)
        # B_k[i] = get_B(dyn, i, k) — fetched per-player inside the loop below

        # Per-player cost matrices at step k
        Q_k = [get_Q(stage_costs[i], k) for i in 1:n_players]
        R_k = [get_R(stage_costs[i], k) for i in 1:n_players]
        q_k = has_affine_terms ? [get_q(stage_costs[i], k) for i in 1:n_players] : nothing
        r_k = has_affine_terms ? [get_r(stage_costs[i], k) for i in 1:n_players] : nothing

        # Concatenated B at step k: [B₁(k) | ... | Bₙ(k)]
        B_k = get_B_concatenated(dyn, k)   # [n × total_controls]

        # ── Build coupled S matrix and RHS ───────────────────────────────────
        S  = zeros(T, total_controls, total_controls)
        YP = zeros(T, total_controls, n)
        Yα = has_affine_terms ? zeros(T, total_controls) : nothing

        for i in 1:n_players
            rng_i  = control_ranges[i]
            Bi_k   = get_B(dyn, i, k)      # [n × mᵢ] at step k
            Ri_k   = R_k[i]
            Zi     = Z[i]
            BiZi   = Bi_k' * Zi            # [mᵢ × n]

            for j in 1:n_players
                rng_j = control_ranges[j]
                Bj_k  = get_B(dyn, j, k)
                # Diagonal: Rᵢ + BᵢᵀZᵢBᵢ; off-diagonal: BᵢᵀZᵢBⱼ
                S[rng_i, rng_j] = (i == j ? Ri_k : zeros(T, length(rng_i), length(rng_j))) + BiZi * Bj_k
            end

            YP[rng_i, :] = BiZi * A_k     # BᵢᵀZᵢA

            if has_affine_terms
                Yα[rng_i] = Bi_k' * ζ[i] + r_k[i]
            end
        end

        # ── Conditioning check ───────────────────────────────────────────────
        if solver.check_singularity
            S_rcond = 1.0 / cond(S)
            push!(rcond_history, S_rcond)
            if S_rcond < solver.rcond_threshold
                @warn "FNELQ: S matrix poorly conditioned at k=$k" rcond=S_rcond
            end
        end

        # ── Solve for feedback gains: S·P = YP ───────────────────────────────
        P   = S \ YP                                    # [total_controls × n]
        P_k = [P[control_ranges[i], :] for i in 1:n_players]
        for i in 1:n_players
            pushfirst!(P_history[i], copy(P_k[i]))
        end

        # ── Affine gains ──────────────────────────────────────────────────────
        α_k_vals = nothing
        β        = nothing
        if has_affine_terms
            α      = S \ Yα
            α_k_vals = [α[control_ranges[i]] for i in 1:n_players]
            for i in 1:n_players
                pushfirst!(α_history[i], copy(α_k_vals[i]))
            end
            β = -B_k * α
        end

        # ── Cost-to-go update: Zᵢ(k) = FᵀZᵢ(k+1)F + Qᵢ(k) + PᵢᵀRᵢ(k)Pᵢ ──
        F = A_k - B_k * P   # Closed-loop dynamics at step k

        for i in 1:n_players
            Pi   = P_k[i]
            Ri_k = R_k[i]
            PRi  = Pi' * Ri_k

            Z[i] = F' * Z[i] * F + Q_k[i] + PRi * Pi
            Z[i] = (Z[i] + Z[i]') / 2      # Enforce symmetry against floating-point drift
            pushfirst!(Z_history[i], copy(Z[i]))

            # ζᵢ(k) = Fᵀ(ζᵢ(k+1) + Zᵢ(k+1)β) + qᵢ(k) + PᵢᵀRᵢ(k)αᵢ - Pᵢᵀrᵢ(k)
            # (Başar & Olsder Eq. 6.17c)
            if has_affine_terms
                αi   = α_k_vals[i]
                ζ[i] = F' * (ζ[i] + Z[i] * β) + q_k[i] + PRi * αi - Pi' * r_k[i]
            end
        end

        verbose && (k % 10 == 0 || k == N) &&
            @info "  Timestep $(N-k+1)/$(N)" rcond=get(rcond_history, length(rcond_history), NaN)
    end

    backward_time = time() - t_start
    verbose && @info "FNELQ: Backward pass complete" time=backward_time

    # ── Forward pass ─────────────────────────────────────────────────────────
    t_forward_start = time()

    x_traj = zeros(T, n, N+1)
    x_traj[:, 1] = game.initial_state
    u_trajs = [zeros(T, control_dims[i], N) for i in 1:n_players]

    for k in 1:N
        x_k = x_traj[:, k]
        A_k = get_A(dyn, k)

        u_k = Vector{Vector{T}}(undef, n_players)
        for i in 1:n_players
            u_k[i] = -P_history[i][k] * x_k
            if has_affine_terms
                u_k[i] -= α_history[i][k]
            end
            u_trajs[i][:, k] = u_k[i]
        end

        # x(k+1) = A(k)·x(k) + Σᵢ Bᵢ(k)·uᵢ(k)
        u_total = sum(get_B(dyn, i, k) * u_k[i] for i in 1:n_players)
        x_traj[:, k+1] = A_k * x_k + u_total
    end

    forward_time = time() - t_forward_start
    total_time   = time() - t_start
    verbose && @info "FNELQ: Forward pass complete" time=forward_time total_time=total_time

    # ── Build trajectories ───────────────────────────────────────────────────
    t_vec = range(T(0), game.time_horizon.tf, length=N+1)

    # Compute per-player costs before building Trajectory (cost field required)
    costs = zeros(T, n_players)
    for i in 1:n_players
        sc = stage_costs[i]
        for k in 1:N
            xk = x_traj[:, k]
            uk = u_trajs[i][:, k]
            Qk = get_Q(sc, k); Rk = get_R(sc, k)
            costs[i] += xk' * Qk * xk + uk' * Rk * uk
            if has_affine_terms
                costs[i] += 2 * get_q(sc, k)' * xk + 2 * get_r(sc, k)' * uk
            end
        end
        xf = x_traj[:, N+1]
        costs[i] += xf' * Qf[i] * xf
    end

    trajectories = Trajectory{T}[]
    for i in 1:n_players
        traj = Trajectory(
            i,
            x_traj,
            u_trajs[i],
            collect(t_vec),
            costs[i]
        )
        push!(trajectories, traj)
    end

    # ── Pack solution ────────────────────────────────────────────────────────
    solver_info = Dict{Symbol, Any}(
        :feedback_gains      => P_history,
        :cost_to_go_matrices => Z_history,
        :costs               => costs,
        :backward_pass_time  => backward_time,
        :forward_pass_time   => forward_time,
        :rcond_history       => rcond_history
    )
    has_affine_terms && (solver_info[:affine_gains] = α_history)

    converged = true
    if solver.check_singularity && !isempty(rcond_history) && any(rcond_history .< solver.rcond_threshold)
        converged = false
        @warn "FNELQ: Solution may be unreliable due to ill-conditioned S matrices"
    end

    return GameSolution(
        game,
        trajectories;
        equilibrium_type = :FeedbackNash,
        converged        = converged,
        iterations       = 1,
        solve_time       = total_time,
        solver_info      = solver_info
    )
end