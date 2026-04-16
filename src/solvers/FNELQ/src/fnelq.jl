# ============================================================================
# FNELQ — Finite-horizon Nash Equilibrium LQ game solver
# (diagnostic prints added in correct positions — NOT inside player loop)
# ============================================================================

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
    @assert(all(obj.stage_cost isa LQStageCost for obj in game.objectives),
        "FNELQ requires LQStageCost for all objectives")
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

    ks = is_ltv(dyn) ? (1:N) : (1:1)

    has_M = any(
        norm(get_M(stage_costs[i], k)) > eps(T)
        for i in 1:np, k in ks
    )
    has_affine = any(
        norm(get_q(stage_costs[i], k)) > eps(T) ||
        norm(get_r(stage_costs[i], k)) > eps(T) ||
        norm(qf_affine[i]) > eps(T)
        for i in 1:np, k in ks
    )

    t_start = time()

    Z = [((Qf[i] + Qf[i]') / 2) for i in 1:np]
    ζ = has_affine ? [copy(qf_affine[i]) for i in 1:np] : nothing

    P_history  = [Matrix{T}[] for _ in 1:np]
    α_history  = has_affine ? [Vector{T}[] for _ in 1:np] : nothing
    Z_history  = [[copy(Z[i])] for i in 1:np]
    rcond_hist = Float64[]

    max_P_norm  = 0.0
    max_Z_norm  = 0.0
    min_rcond   = Inf

    for k in N:-1:1
        Ak  = get_A(dyn, k)
        Bk  = get_B_concatenated(dyn, k)
        Qk  = [get_Q(stage_costs[i], k) for i in 1:np]
        Rk  = [get_R(stage_costs[i], k) for i in 1:np]
        Mk  = has_M      ? [get_M(stage_costs[i], k) for i in 1:np] : nothing
        qk  = has_affine ? [get_q(stage_costs[i], k) for i in 1:np] : nothing
        rk  = has_affine ? [get_r(stage_costs[i], k) for i in 1:np] : nothing

        S  = zeros(T, total_m, total_m)
        YP = zeros(T, total_m, n)
        Yα = has_affine ? zeros(T, total_m) : nothing

        # ── Player loop: assemble S, YP, Yα only ─────────────────────────────
        # NOTE: NO solve, NO pushfirst! inside this loop.
        for i in 1:np
            ri   = control_ranges[i]
            Bi   = get_B(dyn, i, k)
            BiZi = Bi' * Z[i]

            for j in 1:np
                rj = control_ranges[j]
                Bj = get_B(dyn, j, k)
                S[ri, rj] = (i == j ? Rk[i] : zeros(T, length(ri), length(rj))) +
                             BiZi * Bj
            end

            YP[ri, :] = BiZi * Ak
            has_M && (YP[ri, :] += Mk[i]')
            has_affine && (Yα[ri] = -(Bi' * ζ[i] + rk[i]))
        end

        if solver.regularization > zero(T)
            μ_s = T(solver.regularization)
            for i in 1:np
                ri = control_ranges[i]
                for j in ri
                    S[j, j] += μ_s
                end
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

        P_full = S \ YP
        Pk = [P_full[control_ranges[i], :] for i in 1:np]
        for i in 1:np
            pushfirst!(P_history[i], Pk[i])
            max_P_norm = max(max_P_norm, norm(Pk[i]))
        end

        # ── Affine solve: exactly once per timestep, outside player loop ──────
        αk       = nothing
        β_affine = nothing
        if has_affine
            α_full   = S \ Yα
            β_affine = -Bk * α_full
            αk = [α_full[control_ranges[i]] for i in 1:np]

            for i in 1:np
                pushfirst!(α_history[i], αk[i])
            end
        end

        F = Ak - Bk * P_full

        for i in 1:np
            Pi = Pk[i]
            Ri = Rk[i]

            # FIX-1: save Z_old before overwriting
            Z_old_i = Z[i]

            Z[i] = F' * Z[i] * F + Qk[i] + Pi' * Ri * Pi
            if has_M
                Mi   = Mk[i]
                Z[i] -= Mi * Pi + Pi' * Mi'
            end
            Z[i] = (Z[i] + Z[i]') / 2
            pushfirst!(Z_history[i], Z[i])
            max_Z_norm = max(max_Z_norm, norm(Z[i]))

            if has_affine
                ζ[i] = F' * (ζ[i] + Z_old_i * β_affine) +
                        qk[i] + Pi' * Ri * αk[i] - Pi' * rk[i]
            end
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
            ui_opt = -P_history[i][k] * xk
            has_affine && (ui_opt -= α_history[i][k])
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
            costs[i] += evaluate_stage_cost(sc, x_traj[:, k], u_trajs[i][:, k], nothing, k)
        end
        costs[i] += evaluate_terminal_cost(get_objective(game, i).terminal_cost,
                                            x_traj[:, N+1], nothing)
        costs[i] *= get_objective(game, i).scaling
    end

    u_nom_mat = [
        hcat([-P_history[i][k] * x_traj[:, k] for k in 1:N]...)
        for i in 1:np
    ]

    ff = has_affine ? α_history :
                      [[zeros(T, control_dims[i]) for _ in 1:N] for i in 1:np]

    # ─────────────────────────────────────────────────────────────────────────

    strategy = FeedbackStrategy(P_history, ff, x_traj, u_nom_mat, control_dims, t_vec)

    trajectories = [Trajectory(i, x_traj, u_trajs[i], t_vec, costs[i]) for i in 1:np]

    total_time = time() - t_start
    converged  = !(solver.check_singularity &&
                   !isempty(rcond_hist) &&
                   any(rcond_hist .< solver.rcond_threshold))

    if !converged
        @warn "FNELQ: solution may be unreliable (S ill-conditioned)" min_rcond threshold=solver.rcond_threshold
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
            :max_Z_norm           => max_Z_norm,
            :max_P_norm           => max_P_norm,
            :backward_pass_time   => backward_time,
            :forward_pass_time    => forward_time
        )
    )
end