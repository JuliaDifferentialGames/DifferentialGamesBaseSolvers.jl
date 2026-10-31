# ============================================================================
# FNELQ — Finite-horizon Nash Equilibrium LQ game solver
# (diagnostic prints added in correct positions — NOT inside player loop)
# ============================================================================

"""
    FNELQ(; check_singularity=true, rcond_threshold=1e-10, regularization=0.0)

Feedback Nash Equilibrium LQ (FNELQ) solver.

Computes the exact, closed-form feedback Nash equilibrium for finite-horizon discrete-time
LQ games by solving the coupled backward Riccati recursion in a single pass. No iteration
required; the solution is exact up to floating-point precision.

Applicable to `LQGameProblem`, `LTVLQGameProblem`, and `PDGNEProblem` with LQ costs.

# Options
- `check_singularity`: warn when the Nash gain matrix is ill-conditioned
- `rcond_threshold`: reciprocal condition number below which a warning is issued
- `regularization`: Tikhonov regularization added to the gain system
"""
struct FNELQ <: GameSolver
    check_singularity::Bool
    rcond_threshold::Float64
    regularization::Float64

    function FNELQ(;
        check_singularity::Bool  = true,
        rcond_threshold::Float64 = 1e-10,
        regularization::Float64  = 0.0
    )
        new(check_singularity, rcond_threshold, regularization)
    end
end

solver_capabilities(::Type{FNELQ}) = [
    :LQGame,
    :FeedbackPolicies,
    :UnconstrainedGame,
    :DiscreteTime
]

function extract_lq_matrices(game::GameProblem{T}) where {T}
    @assert(game.dynamics isa LinearDynamics,
        "FNELQ requires LinearDynamics; got $(typeof(game.dynamics))")
    @assert(all(obj.stage_cost isa Union{LQStageCost, JointLQStageCost} for obj in game.objectives),
        "FNELQ requires LQStageCost or JointLQStageCost for all objectives")
    @assert(all(obj.terminal_cost isa LQTerminalCost for obj in game.objectives),
        "FNELQ requires LQTerminalCost for all objectives")

    dyn          = game.dynamics
    np           = game.n_players
    stage_costs  = [get_objective(game, i).stage_cost for i in 1:np]
    tc           = [get_objective(game, i).terminal_cost for i in 1:np]
    Qf           = [tc[i].Qf for i in 1:np]
    qf_affine    = [tc[i].qf for i in 1:np]

    return dyn, stage_costs, Qf, qf_affine, total_state_dim(dyn), np
end

# Player i's stage cost as joint-control matrices (Q, R, M, q, r) with R: m×m,
# M: n×m, r: m. An own-control `LQStageCost` is embedded in player i's control
# block; a `JointLQStageCost` is used as is.
function _joint_cost_mats(sc::JointLQStageCost, i, k, ranges, n, m, ::Type{T}) where {T}
    return get_Q(sc, k), get_R(sc, k), get_M(sc, k), get_q(sc, k), get_r(sc, k)
end

function _joint_cost_mats(sc::LQStageCost, i, k, ranges, n, m, ::Type{T}) where {T}
    ri = ranges[i]
    R = zeros(T, m, m); R[ri, ri] = get_R(sc, k)
    M = zeros(T, n, m); M[:, ri] = get_M(sc, k)
    r = zeros(T, m);    r[ri]    = get_r(sc, k)
    return get_Q(sc, k), R, M, get_q(sc, k), r
end

# Stage cost of player i given the joint control (own-control costs take u^i only).
_stage_cost_value(sc::JointLQStageCost, i, ranges, x, u_joint, k) =
    evaluate_stage_cost(sc, x, u_joint, nothing, k)
_stage_cost_value(sc::LQStageCost, i, ranges, x, u_joint, k) =
    evaluate_stage_cost(sc, x, u_joint[ranges[i]], nothing, k)
_stage_cost_value(sc::AbstractStageCost, i, ranges, x, u_joint, k) =
    evaluate_stage_cost(sc, x, u_joint[ranges[i]], nothing, k)

# Result returned when a stage's stacked Nash system cannot be solved: zero
# policies, `converged = false`, and an explicit status instead of a false equilibrium.
function _failed_fnelq_solution(game::GameProblem{T}, solver, status, k, t_vec, control_dims, elapsed) where {T}
    np = game.n_players; n = length(game.initial_state); N = length(t_vec) - 1
    x_traj = repeat(game.initial_state, 1, N + 1)
    u_trajs = [zeros(T, control_dims[i], N) for i in 1:np]
    gains = [[zeros(T, control_dims[i], n) for _ in 1:N] for i in 1:np]
    ff = [[zeros(T, control_dims[i]) for _ in 1:N] for i in 1:np]
    strategy = FeedbackStrategy(gains, ff, x_traj, [zeros(T, control_dims[i], N) for i in 1:np], control_dims, t_vec)
    trajs = [Trajectory(i, x_traj, u_trajs[i], t_vec, T(NaN)) for i in 1:np]
    @warn "FNELQ: stage Nash system is singular; no equilibrium computed" k
    return GNEPSolution(game, trajs; state_trajectory = x_traj, strategy = strategy,
        equilibrium_type = :FeedbackNash, converged = false, iterations = 1, solve_time = elapsed,
        solver_info = Dict{Symbol, Any}(:status => status, :failed_at_step => k,
            :rcond_history => Float64[], :min_rcond => 0.0))
end

function _solve(
    game::GameProblem{T},
    solver::FNELQ,
    warmstart::Union{Nothing, WarmstartData},
    verbose::Bool
) where {T}

    dyn, stage_costs, Qf, qf_affine, n, np = extract_lq_matrices(game)

    control_dims   = dyn.control_dims
    total_m        = sum(control_dims)
    c_offs         = [0; cumsum(control_dims)[1:end-1]]
    control_ranges = [c_offs[i]+1 : c_offs[i]+control_dims[i] for i in 1:np]

    N     = n_steps(game)
    t_vec = collect(range(T(0), game.time_horizon.tf, length=N+1))

    m_tot = total_m
    mats(i, k) = _joint_cost_mats(stage_costs[i], i, k, control_ranges, n, m_tot, T)
    t_start = time()

    Z = [((Qf[i] + Qf[i]') / 2) for i in 1:np]
    ζ = [copy(qf_affine[i]) for i in 1:np]

    P_history  = [Matrix{T}[] for _ in 1:np]
    α_history  = [Vector{T}[] for _ in 1:np]
    Z_history  = [[copy(Z[i])] for i in 1:np]
    rcond_hist = Float64[]
    min_own_eig = Inf   # smallest eigenvalue of any player's own-control Hessian block

    max_P_norm  = 0.0
    max_Z_norm  = 0.0
    min_rcond   = Inf

    for k in N:-1:1
        Ak  = get_A(dyn, k)
        Bk  = get_B_concatenated(dyn, k)
        cm  = [mats(i, k) for i in 1:np]
        Qk  = [c[1] for c in cm]; Rk = [c[2] for c in cm]; Mk = [c[3] for c in cm]
        qk  = [c[4] for c in cm]; rk = [c[5] for c in cm]

        # Stationarity of player i in its own control block: row block i of
        #   (R_i + BᵀZ_iB) u + (M_iᵀ + BᵀZ_iA) x + (r_i + Bᵀζ_i) = 0.
        S  = zeros(T, m_tot, m_tot)
        YP = zeros(T, m_tot, n)
        Yα = zeros(T, m_tot)
        for i in 1:np
            ri = control_ranges[i]
            BZ = Bk' * Z[i]
            S[ri, :]  = (Rk[i] + BZ * Bk)[ri, :]
            YP[ri, :] = (Mk[i]' + BZ * Ak)[ri, :]
            Yα[ri]    = (rk[i] + Bk' * ζ[i])[ri]
            # Second-order condition: own-control block must be positive definite.
            Sii = Symmetric(S[ri, ri])
            min_own_eig = min(min_own_eig, minimum(eigvals(Sii)))
        end

        if solver.regularization > zero(T)
            μ_s = T(solver.regularization)
            for j in 1:m_tot
                S[j, j] += μ_s
            end
        end

        rc = NaN
        if solver.check_singularity
            rc = 1.0 / cond(S)
            push!(rcond_hist, rc)
            min_rcond = min(min_rcond, rc)
            if rc < solver.rcond_threshold
                @warn "FNELQ: S ill-conditioned at k=$k" rcond=rc μ=solver.regularization
            end
        end

        local P_full, α_full
        try
            P_full = S \ YP
            α_full = S \ Yα
        catch err
            err isa LinearAlgebra.SingularException || rethrow()
            return _failed_fnelq_solution(game, solver, :singular_nash_system, k, t_vec, control_dims, time() - t_start)
        end
        if !all(isfinite, P_full) || !all(isfinite, α_full)
            return _failed_fnelq_solution(game, solver, :singular_nash_system, k, t_vec, control_dims, time() - t_start)
        end
        β_affine = -Bk * α_full
        Pk = [P_full[control_ranges[i], :] for i in 1:np]
        αk = [α_full[control_ranges[i]] for i in 1:np]
        for i in 1:np
            pushfirst!(P_history[i], Pk[i])
            pushfirst!(α_history[i], αk[i])
            max_P_norm = max(max_P_norm, norm(Pk[i]))
        end

        F = Ak - Bk * P_full
        for i in 1:np
            Z_old_i = Z[i]
            # Substituting u = -P x - α into the stage cost and next value function.
            Z[i] = F' * Z_old_i * F + Qk[i] + P_full' * Rk[i] * P_full -
                   Mk[i] * P_full - P_full' * Mk[i]'
            Z[i] = (Z[i] + Z[i]') / 2
            pushfirst!(Z_history[i], Z[i])
            max_Z_norm = max(max_Z_norm, norm(Z[i]))
            ζ[i] = F' * (ζ[i] + Z_old_i * β_affine) + qk[i] +
                   P_full' * Rk[i] * α_full - P_full' * rk[i] - Mk[i] * α_full
        end

        if verbose && (k % max(1, N ÷ 5) == 0 || k == 1 || k == N)
            @info "  FNELQ k=$k" rcond=rc Z_norm=max_Z_norm P_norm=max_P_norm
        end
    end

    backward_time = time() - t_start

    t_fwd  = time()
    x_traj = zeros(T, n, N+1)
    x_traj[:, 1] = game.initial_state
    u_trajs = [zeros(T, control_dims[i], N) for i in 1:np]

    for k in 1:N
        xk = x_traj[:, k]
        B_u_sum = zeros(T, n)
        for i in 1:np
            ui_opt = -P_history[i][k] * xk - α_history[i][k]
            u_trajs[i][:, k] = ui_opt
            B_u_sum += get_B(dyn, i, k) * ui_opt
        end
        x_traj[:, k+1] = get_A(dyn, k) * xk + B_u_sum
    end

    forward_time = time() - t_fwd

    costs = zeros(T, np)
    for i in 1:np
        sc = stage_costs[i]
        for k in 1:N
            u_joint = vcat((u_trajs[j][:, k] for j in 1:np)...)
            costs[i] += _stage_cost_value(sc, i, control_ranges, x_traj[:, k], u_joint, k)
        end
        costs[i] += evaluate_terminal_cost(get_objective(game, i).terminal_cost,
                                            x_traj[:, N+1], nothing)
        costs[i] *= get_objective(game, i).scaling
    end

    u_nom_mat = [
        hcat([-P_history[i][k] * x_traj[:, k] for k in 1:N]...)
        for i in 1:np
    ]

    ff = α_history

    # ─────────────────────────────────────────────────────────────────────────

    strategy = FeedbackStrategy(P_history, ff, x_traj, u_nom_mat, control_dims, t_vec)

    trajectories = [Trajectory(i, x_traj, u_trajs[i], t_vec, costs[i]) for i in 1:np]

    total_time = time() - t_start
    ill_conditioned = solver.check_singularity && !isempty(rcond_hist) &&
                      any(rcond_hist .< solver.rcond_threshold)
    non_nash = !(min_own_eig > 0)
    status = ill_conditioned ? :ill_conditioned_nash_system :
             non_nash        ? :second_order_violation : :solved
    converged = status === :solved

    if ill_conditioned
        @warn "FNELQ: solution may be unreliable (S ill-conditioned)" min_rcond threshold=solver.rcond_threshold
    elseif non_nash
        @warn "FNELQ: a player's own-control Hessian block is not positive definite; the stationary point is not a Nash equilibrium" min_own_eig
    end

    return GNEPSolution(
        game,
        trajectories;
        state_trajectory = x_traj,
        strategy         = strategy,
        equilibrium_type = :FeedbackNash,
        converged        = converged,
        iterations       = 1,
        solve_time       = total_time,
        solver_info      = Dict{Symbol, Any}(
            :cost_to_go_matrices  => Z_history,
            :rcond_history        => rcond_hist,
            :min_rcond            => min_rcond,
            :status               => status,
            :min_own_hessian_eig  => min_own_eig,
            :second_order_ok      => min_own_eig > 0,
            :max_Z_norm           => max_Z_norm,
            :max_P_norm           => max_P_norm,
            :backward_pass_time   => backward_time,
            :forward_pass_time    => forward_time
        )
    )
end