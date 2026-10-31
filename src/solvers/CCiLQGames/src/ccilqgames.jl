# ============================================================================
# CCiLQGames — Chance-Constrained iterative LQ Stochastic Games
#
# Follows the nested structure of Zhong, Shimizu & Chen, "Chance-Constrained
# Iterative Linear-Quadratic Stochastic Games", RA-L 2023:
#
#   outer loop  augmented-Lagrangian multiplier / penalty updates
#   inner loop  nominal stochastic rollout → linearisation → Gaussian chance
#               residuals → local LQ game (solved by FNELQ) → line search
#
# GUARANTEES (see the `solver_info` metadata):
#   * Each inner step solves a local LQ game exactly; no global equilibrium or
#     convergence theorem is claimed for the nonlinear stochastic game.
#   * Chance residuals are the Gaussian-quantile surrogate
#         c = g(ψ̄) + Φ⁻¹(1-ε) · sqrt(∇gᵀ Σψ ∇g)
#     exact for affine g under the Gaussian law, a first-order approximation
#     otherwise (`guarantee_classes`). Probability statements are only as
#     strong as that class; the Monte Carlo summary is an estimate.
#   * `covariance_mode = :frozen_covariance`: within each local game the
#     covariance-dependent term σ is held fixed (its derivative w.r.t. the
#     nominal and the policy is neglected), and the covariance trace term of
#     the expected cost (½tr(QΣ)) is dropped. Stationarity is therefore that of
#     the frozen-covariance surrogate, NOT of the full stochastic game.
#   * Perfect state information only; constraints are imposed on the physical
#     state/control at stages k = 1…N (not on the terminal state).
# ============================================================================

"""
    CCiLQGames(; max_iter=20, max_inner_iter=30, constraint_tol=1e-4,
                 stationarity_tol=1e-3, ρ0=10.0, ρ_growth=5.0, ρ_max=1e4,
                 regularization=1e-6, covariance_mode=:frozen_covariance,
                 evaluator=GaussianQuantileEvaluator(), line_search_β=0.5,
                 η_min=1e-3, max_state_step=1e3, mc_samples=2000, seed=0)

Solve a [`StochasticGameProblem`](@ref) with additive-or-state-dependent Gaussian
noise and private/shared chance constraints for a locally stationary feedback
Nash policy of the Gaussian-surrogate game. Returns a `GNEPSolution` whose
`solver_info` carries the status, residual histories, covariance trajectories,
risk allocation, guarantee classes, multipliers, and (if `mc_samples > 0`) an
empirical constraint-violation summary.

`covariance_mode = :consistent_covariance` is reserved and not implemented.
"""
struct CCiLQGames <: GameSolver
    max_iter::Int
    max_inner_iter::Int
    constraint_tol::Float64
    stationarity_tol::Float64
    ρ0::Float64
    ρ_growth::Float64
    ρ_max::Float64
    regularization::Float64
    covariance_mode::Symbol
    evaluator::AbstractChanceEvaluator
    line_search_β::Float64
    η_min::Float64
    max_state_step::Float64
    mc_samples::Int
    seed::Int
end

function CCiLQGames(;
    max_iter::Int = 20, max_inner_iter::Int = 30,
    constraint_tol::Real = 1e-4, stationarity_tol::Real = 1e-3,
    ρ0::Real = 10.0, ρ_growth::Real = 5.0, ρ_max::Real = 1e4,
    regularization::Real = 1e-6, covariance_mode::Symbol = :frozen_covariance,
    evaluator::AbstractChanceEvaluator = GaussianQuantileEvaluator(),
    line_search_β::Real = 0.5, η_min::Real = 1e-3, max_state_step::Real = 1e3,
    mc_samples::Int = 2000, seed::Int = 0
)
    covariance_mode in (:frozen_covariance, :consistent_covariance) ||
        throw(ArgumentError("covariance_mode must be :frozen_covariance or :consistent_covariance"))
    CCiLQGames(max_iter, max_inner_iter, Float64(constraint_tol), Float64(stationarity_tol),
        Float64(ρ0), Float64(ρ_growth), Float64(ρ_max), Float64(regularization), covariance_mode,
        evaluator, Float64(line_search_β), Float64(η_min), Float64(max_state_step), mc_samples, seed)
end

solver_capabilities(::Type{CCiLQGames}) =
    [:Stochastic, :ChanceConstrained, :FeedbackPolicies, :DiscreteTime, :NonlinearGame, :LQGame]

# ─── problem context ─────────────────────────────────────────────────────────

# One block of constraint events (a ChanceConstraint or a deterministic inequality).
struct _Spec
    kind::Symbol                 # :chance | :det
    players::Vector{Int}         # players whose augmented cost carries it
    inner                        # underlying (non-chance) constraint
    dim::Int
    eps::Matrix{Float64}         # dim × N risk levels (unused for :det)
    guarantee::Symbol
end

struct _Ctx{T, SD, EV}
    game::GameProblem{T}
    sg::StochasticGameProblem{T}
    sd::SD
    n::Int; m::Int; np::Int; N::Int; dt::T
    ranges::Vector{UnitRange{Int}}
    linear_drift::Bool
    specs::Vector{_Spec}
    evaluator::EV
end

_tk(ctx::_Ctx, k) = ctx.linear_drift ? k : (k - 1) * ctx.dt
_W(ctx::_Ctx, k, x, u) = ctx.sg.process_covariances === nothing ?
    process_covariance(ctx.sd, x, u, nothing, _tk(ctx, k)) : ctx.sg.process_covariances[k]

function _build_specs(sg::StochasticGameProblem{T}, N::Int, evaluator) where {T}
    g = sg.nominal_game
    ccs = all_chance_constraints(sg)
    budget = sg.joint_risk === nothing ? nothing :
        allocate_horizon_risk(sg.risk_allocation, ccs, N, sg.joint_risk)
    specs = _Spec[]
    for (ci, cc) in enumerate(ccs)
        inner = cc.constraint
        dim = constraint_output_dim(inner)
        eps = budget === nothing ?
            repeat(reshape(Float64.(allocate_risk(sg.risk_allocation, cc)), dim, 1), 1, N) :
            Float64.(budget[ci])
        players = is_private(cc) ? [get_player(cc)] : collect(get_players(cc))
        push!(specs, _Spec(:chance, players, inner, dim, eps, guarantee_class(cc, evaluator)))
    end
    for c in Iterators.flatten((g.private_constraints, g.shared_constraints))
        c isa ChanceConstraint && continue
        is_inequality(c) || error("CCiLQGames supports inequality constraints only; got $(typeof(c))")
        players = is_private(c) ? [get_player(c)] : collect(get_players(c))
        push!(specs, _Spec(:det, players, c, constraint_output_dim(c), zeros(constraint_output_dim(c), N), :deterministic))
    end
    return specs
end

# ─── dynamics: nominal rollout, linearisation, covariance ────────────────────

_step_mean(ctx, x, u, k) = mean_step(ctx.sd, x, u, nothing, _tk(ctx, k))

function _rollout_policy(ctx, X̄, Ū, P, α, η, x0)
    X = zeros(ctx.n, ctx.N + 1); U = zeros(ctx.m, ctx.N)
    X[:, 1] = x0
    for k in 1:ctx.N
        δx = X[:, k] - X̄[:, k]
        U[:, k] = Ū[:, k] - P[k] * δx - η * α[k]
        X[:, k+1] = _step_mean(ctx, X[:, k], U[:, k], k)
        (all(isfinite, X[:, k+1]) && maximum(abs.(X[:, k+1] - X̄[:, k+1])) < 1e6) || return nothing
    end
    return X, U
end

function _linearize(ctx, X, U)
    As = Matrix{Float64}[]; Bs = Matrix{Float64}[]
    for k in 1:ctx.N
        Jx, Ju = jacobian(ctx.sd, X[:, k], U[:, k], nothing, _tk(ctx, k))
        push!(As, Matrix{Float64}(Jx)); push!(Bs, Matrix{Float64}(Ju))
    end
    return As, Bs
end

# Physical-state/control covariances under u_k = ū_k - K_k (x_k - x̄_k), Σ₀ = initial covariance.
function _moments(ctx, As, Bs, K, X, U)
    Σ = Matrix(ctx.sg.initial_covariance)
    Σx = Matrix{Float64}[]; Σu = Matrix{Float64}[]; Σxu = Matrix{Float64}[]
    for k in 1:ctx.N
        push!(Σx, copy(Σ)); push!(Σu, Matrix(Symmetric(K[k] * Σ * K[k]')))
        push!(Σxu, -Σ * K[k]')
        F = As[k] - Bs[k] * K[k]
        Σ = Matrix(Symmetric(F * Σ * F' + _W(ctx, k, X[:, k], U[:, k])))
    end
    push!(Σx, Σ)
    return Σx, Σu, Σxu
end

# ─── chance residual events ──────────────────────────────────────────────────

function _events(ctx, spec::_Spec, x, u, k, Σx, Σu, Σxu)
    c = Vector{Float64}(evaluate_constraint(spec.inner, x, u, nothing, k))
    Jx, Ju = constraint_jacobian(spec.inner, x, u, nothing, k)
    res = copy(c)
    if spec.kind === :chance
        for i in 1:spec.dim
            gx = Jx[i, :]; gu = Ju[i, :]
            v = dot(gx, Σx * gx) + 2 * dot(gx, Σxu * gu) + dot(gu, Σu * gu)
            res[i] = chance_residual(ctx.evaluator, c[i], v, spec.eps[i, k])
        end
    end
    return res, Matrix{Float64}(Jx), Matrix{Float64}(Ju), c
end

function _feasibility(ctx, X, U, mom)
    Σx, Σu, Σxu = mom
    r = 0.0
    for spec in ctx.specs, k in 1:ctx.N
        res = _events(ctx, spec, X[:, k], U[:, k], k, Σx[k], Σu[k], Σxu[k])[1]
        r = max(r, maximum(max.(res, 0.0)))
    end
    return r
end

# ─── local LQ game ───────────────────────────────────────────────────────────

function _psd_clip(H::AbstractMatrix, floor::Real)
    S = Symmetric((H + H') / 2)
    F = eigen(S)
    return Matrix(Symmetric(F.vectors * Diagonal(max.(F.values, floor)) * F.vectors'))
end

function _local_game(ctx, solver, X, U, As, Bs, mom, λ, ρ)
    g = ctx.game; n, m, np, N = ctx.n, ctx.m, ctx.np, ctx.N
    Σx, Σu, Σxu = mom
    Qs = [Matrix{Float64}[] for _ in 1:np]; Rs = [Matrix{Float64}[] for _ in 1:np]
    Ms = [Matrix{Float64}[] for _ in 1:np]; qs = [Vector{Float64}[] for _ in 1:np]
    rs = [Vector{Float64}[] for _ in 1:np]
    for k in 1:N
        x = X[:, k]; u = U[:, k]
        ev = [_events(ctx, spec, x, u, k, Σx[k], Σu[k], Σxu[k]) for spec in ctx.specs]
        for i in 1:np
            obj = get_objective(g, i)
            f = z -> obj.scaling * _stage_cost_value(obj.stage_cost, i, ctx.ranges, z[1:n], z[n+1:end], k)
            z0 = vcat(x, u)
            gz = ForwardDiff.gradient(f, z0); Hz = ForwardDiff.hessian(f, z0)
            for (si, spec) in enumerate(ctx.specs)
                i in spec.players || continue
                res, Jx, Ju, _ = ev[si]
                for mm in 1:spec.dim
                    s = λ[si][mm, k] + ρ[si][mm, k] * res[mm]
                    s > 0 || continue
                    ∇ = vcat(Jx[mm, :], Ju[mm, :])
                    gz = gz + s * ∇
                    Hz = Hz + ρ[si][mm, k] * (∇ * ∇')
                end
            end
            Q = _psd_clip(Hz[1:n, 1:n], 0.0)
            R = Matrix((Hz[n+1:end, n+1:end] + Hz[n+1:end, n+1:end]') / 2)
            ri = ctx.ranges[i]
            R[ri, ri] = _psd_clip(R[ri, ri], 0.0) + solver.regularization * I
            R = Matrix(Symmetric(R))
            push!(Qs[i], Q); push!(Rs[i], R); push!(Ms[i], Matrix(Hz[1:n, n+1:end]))
            push!(qs[i], gz[1:n]); push!(rs[i], gz[n+1:end])
        end
    end
    objectives = map(1:np) do i
        obj = get_objective(g, i)
        φ = x -> obj.scaling * evaluate_terminal_cost(obj.terminal_cost, x, nothing)
        xN = X[:, N+1]
        Qf = _psd_clip(ForwardDiff.hessian(φ, xN), 0.0)
        qf = ForwardDiff.gradient(φ, xN)
        PlayerObjective(i, JointLQStageCost(Qs[i], Rs[i]; M = Ms[i], q = qs[i], r = rs[i]),
                        LQTerminalCost(Qf, Vector{Float64}(qf), 0.0))
    end
    Bi = [[Matrix(Bs[k][:, ctx.ranges[i]]) for k in 1:N] for i in 1:np]
    dyn = LinearDynamics(As, Bi)
    cd = [length(r) for r in ctx.ranges]
    meta = GameMetadata([n], cd, [0], [first(r) - 1 for r in ctx.ranges],
        CouplingGraph(sparse(trues(np, np)), Vector{Int}[], nothing), false, nothing)
    return GameProblem{Float64}(np, objectives, dyn, zeros(n),
        AbstractPrivateConstraint[], AbstractSharedConstraint[], DiscreteTime(N * ctx.dt, ctx.dt), meta)
end

# ─── costs, validation ───────────────────────────────────────────────────────

function _nominal_costs(ctx, X, U)
    c = zeros(ctx.np)
    for i in 1:ctx.np
        obj = get_objective(ctx.game, i)
        for k in 1:ctx.N
            c[i] += _stage_cost_value(obj.stage_cost, i, ctx.ranges, X[:, k], U[:, k], k)
        end
        c[i] += evaluate_terminal_cost(obj.terminal_cost, X[:, ctx.N+1], nothing)
        c[i] *= obj.scaling
    end
    return c
end

function _monte_carlo(ctx, solver, X, U, K)
    chance = [s for s in ctx.specs if s.kind === :chance]
    (isempty(chance) || solver.mc_samples <= 0) && return nothing
    rng = Random.MersenneTwister(solver.seed)
    L0 = DifferentialGamesBase._psd_factor(Matrix(ctx.sg.initial_covariance))
    n_events = sum(s.dim for s in chance) * ctx.N
    G = zeros(n_events, solver.mc_samples)
    for s in 1:solver.mc_samples
        x = X[:, 1] + L0 * randn(rng, size(L0, 2)); row = 0
        for k in 1:ctx.N
            u = U[:, k] - K[k] * (x - X[:, k])
            for spec in chance
                c = evaluate_constraint(spec.inner, x, u, nothing, k)
                G[row+1:row+spec.dim, s] = c; row += spec.dim
            end
            x = sample_step(rng, ctx.sd, x, u, nothing, _tk(ctx, k))
        end
    end
    return empirical_violation_rates(G)
end

# ─── _solve ──────────────────────────────────────────────────────────────────

function _solve(sg::StochasticGameProblem{T}, solver::CCiLQGames, warmstart, verbose::Bool) where {T}
    t_start = time()
    solver.covariance_mode === :consistent_covariance && error(
        "CCiLQGames: covariance_mode = :consistent_covariance is not implemented; use :frozen_covariance.")
    sg.information isa PerfectStateInformation || error(
        "CCiLQGames supports PerfectStateInformation only; got $(typeof(sg.information)).")
    sde = sg.stochastic_dynamics
    sde isa ItoSDEDynamics || error("CCiLQGames requires ItoSDEDynamics stochastic dynamics.")
    g = sg.nominal_game
    g isa GameProblem || error("CCiLQGames requires a GameProblem nominal game.")

    N = n_steps(g); dt = g.time_horizon.dt
    n = total_state_dim(sde); m = total_control_dim(sde); np = num_players(g)
    cd = g.metadata.control_dims
    offs = [0; cumsum(cd)[1:end-1]]
    ranges = [offs[i]+1:offs[i]+cd[i] for i in 1:np]
    sd = discretize_stochastic(sde, dt)
    specs = _build_specs(sg, N, solver.evaluator)
    ctx = _Ctx(g, sg, sd, n, m, np, N, dt, ranges, sde.drift isa LinearDynamics, specs, solver.evaluator)

    x0 = Vector{Float64}(g.initial_state)
    zerosP() = [zeros(m, n) for _ in 1:N]
    K = zerosP()
    X = zeros(n, N + 1); U = zeros(m, N); X[:, 1] = x0
    for k in 1:N
        X[:, k+1] = _step_mean(ctx, X[:, k], U[:, k], k)
    end
    λ = [zeros(s.dim, N) for s in specs]
    ρ = [fill(solver.ρ0, s.dim, N) for s in specs]

    hist = (r_stat = Float64[], r_feas = Float64[], r_step = Float64[], η = Float64[],
            cost = Vector{Vector{Float64}}(), outer_of_inner = Int[])
    outer_hist = (r_feas = Float64[], ρ_max = Float64[])
    status = :max_iterations
    P_last = nothing
    prev_feas = Inf
    r_stat = Inf

    for outer in 1:solver.max_iter
        η_cap = 1.0; prev_stat = Inf     # damping: shrinks while stationarity stalls (breaks 2-cycles)
        for inner in 1:solver.max_inner_iter
            As, Bs = _linearize(ctx, X, U)
            mom = _moments(ctx, As, Bs, K, X, U)
            lq = _local_game(ctx, solver, X, U, As, Bs, mom, λ, ρ)
            det = _solve(lq, FNELQ(check_singularity = true, rcond_threshold = 1e-12), nothing, false)
            st = det.solver_info[:status]
            if st !== :solved
                status = st === :singular_nash_system ? :singular_nash_system : :ill_conditioned_nash_system
                st === :second_order_violation && (status = :second_order_violation)
                @warn "CCiLQGames: local LQ game failed" status outer inner
                return _ccil_solution(ctx, solver, g, X, U, K, mom, specs, λ, ρ, hist, outer_hist, status, t_start, 0.0, r_stat, outer)
            end
            P = [vcat((det.strategy.gains[i][k] for i in 1:np)...) for k in 1:N]
            α = [vcat((det.strategy.feedforward[i][k] for i in 1:np)...) for k in 1:N]
            r_stat = maximum(maximum(abs, a) for a in α)
            r_stat > 0.9 * prev_stat && (η_cap = max(0.02, η_cap / 2))
            prev_stat = r_stat

            # Backtracking on a feasible-and-bounded step (general-sum games have no scalar merit).
            η = η_cap; new = nothing
            while η >= solver.η_min
                cand = _rollout_policy(ctx, X, U, P, α, η, x0)
                if cand !== nothing && maximum(abs.(cand[1] - X)) <= solver.max_state_step
                    new = cand; break
                end
                η *= solver.line_search_β
            end
            if new === nothing
                status = :line_search_failed
                return _ccil_solution(ctx, solver, g, X, U, K, mom, specs, λ, ρ, hist, outer_hist, status, t_start, 0.0, r_stat, outer)
            end
            Xn, Un = new
            r_step = max(maximum(abs.(Xn - X)), maximum(abs.(Un - U)))
            X, U, K = Xn, Un, P
            As, Bs = _linearize(ctx, X, U)
            mom = _moments(ctx, As, Bs, K, X, U)
            push!(hist.r_stat, r_stat); push!(hist.r_step, r_step); push!(hist.η, η)
            push!(hist.r_feas, _feasibility(ctx, X, U, mom)); push!(hist.cost, _nominal_costs(ctx, X, U))
            push!(hist.outer_of_inner, outer)
            verbose && @info "CCiLQGames" outer inner r_stat r_feas = hist.r_feas[end] η
            r_stat <= solver.stationarity_tol && break
        end

        As, Bs = _linearize(ctx, X, U)
        mom = _moments(ctx, As, Bs, K, X, U)
        r_feas = _feasibility(ctx, X, U, mom)
        push!(outer_hist.r_feas, r_feas); push!(outer_hist.ρ_max, maximum(maximum(r) for r in ρ; init = 0.0))
        if r_feas <= solver.constraint_tol
            if r_stat <= solver.stationarity_tol
                status = :locally_stationary_surrogate; break
            end
        end
        # Multiplier and penalty updates (inequality augmented Lagrangian).
        Σx, Σu, Σxu = mom
        for (si, spec) in enumerate(specs), k in 1:N
            res = _events(ctx, spec, X[:, k], U[:, k], k, Σx[k], Σu[k], Σxu[k])[1]
            for mm in 1:spec.dim
                λ[si][mm, k] = max(0.0, λ[si][mm, k] + ρ[si][mm, k] * res[mm])
                if res[mm] > solver.constraint_tol && r_feas > 0.25 * prev_feas
                    ρ[si][mm, k] = min(solver.ρ_max, solver.ρ_growth * ρ[si][mm, k])
                end
            end
        end
        prev_feas = r_feas
        if outer == solver.max_iter
            status = r_feas <= solver.constraint_tol ? :surrogate_feasible : :max_iterations
        end
    end
    As, Bs = _linearize(ctx, X, U)
    mom = _moments(ctx, As, Bs, K, X, U)
    return _ccil_solution(ctx, solver, g, X, U, K, mom, specs, λ, ρ, hist, outer_hist, status, t_start,
                          _feasibility(ctx, X, U, mom), r_stat, length(outer_hist.r_feas))
end

function _ccil_solution(ctx, solver, g, X, U, K, mom, specs, λ, ρ, hist, outer_hist, status, t_start, r_feas, r_stat, outer)
    np, N = ctx.np, ctx.N
    cd = [length(r) for r in ctx.ranges]
    t_vec = collect(range(0.0, N * ctx.dt, length = N + 1))
    U_pp = [U[ctx.ranges[i], :] for i in 1:np]
    gains = [[K[k][ctx.ranges[i], :] for k in 1:N] for i in 1:np]
    ff = [[zeros(cd[i]) for _ in 1:N] for i in 1:np]
    strategy = FeedbackStrategy(gains, ff, X, U_pp, cd, t_vec)
    costs = _nominal_costs(ctx, X, U)
    trajs = [Trajectory(i, X, U_pp[i], t_vec, costs[i]) for i in 1:np]
    Σx, Σu, Σxu = mom

    chance_specs = [(si, s) for (si, s) in enumerate(specs) if s.kind === :chance]
    ccs = all_chance_constraints(ctx.sg)
    checks = isempty(ccs) ? ChanceCheck[] :
        check_chance_constraints(ccs, X[:, 1:N], U, Σx[1:N], Σu, Σxu;
                                 risk = [s.eps for (_, s) in chance_specs], evaluator = ctx.evaluator)
    mc = status in (:singular_nash_system, :ill_conditioned_nash_system, :second_order_violation, :line_search_failed) ?
        nothing : _monte_carlo(ctx, solver, X, U, K)
    converged = status in (:locally_stationary_surrogate, :surrogate_feasible)

    info = Dict{Symbol, Any}(
        :status => status,
        :covariance_mode => solver.covariance_mode,
        :frozen_terms => solver.covariance_mode === :frozen_covariance ?
            [:sigma_derivative_wrt_nominal_and_policy, :expected_cost_trace_term] : Symbol[],
        :r_stat => r_stat, :r_feas => r_feas,
        :residual_history => hist, :outer_history => outer_hist,
        :outer_iterations => outer,
        :state_covariances => Σx, :control_covariances => Σu, :state_control_covariances => Σxu,
        :multipliers => λ, :penalties => ρ,
        :risk_allocation => [s.eps for (_, s) in chance_specs],
        :joint_risk_bound => isempty(chance_specs) ? 0.0 : sum(sum(s.eps) for (_, s) in chance_specs),
        :guarantee_classes => [s.guarantee for (_, s) in chance_specs],
        :chance_checks => checks,
        :analytic_joint_violation_bound => joint_violation_bound(checks),
        :monte_carlo => mc,
        :expected_costs => costs)
    return GNEPSolution(g, trajs; state_trajectory = X, strategy = strategy,
        equilibrium_type = :FeedbackNash, converged = converged,
        iterations = length(hist.r_stat), solve_time = time() - t_start, solver_info = info)
end
