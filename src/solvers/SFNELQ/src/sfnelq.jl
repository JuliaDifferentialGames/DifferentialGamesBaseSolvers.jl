# ============================================================================
# SFNELQ — Stochastic Feedback Nash Equilibrium, Linear-Quadratic.
#
# Ported from Alessandro-Work/sgnep/{SGNEP.jl,SFNELQ.jl}, reconciled with
# DifferentialGamesBase.jl's reworked StochasticGameProblem (which now wraps
# an ItoSDEDynamics rather than a standalone NoiseModel — see
# DifferentialGamesBase.jl's docs/src/api/stochastic.md).
#
# SCOPE OF THIS PORT: additive noise only (StochasticGameProblem built with
# `ConstantDiffusion`, i.e. `is_certainty_equivalent(sgame) == true`).
# Alessandro's original SFNELQ also handled state-/control-multiplicative
# noise (NoiseModel's `C`/`D` fields); DifferentialGamesBase.jl's current
# ItoSDEDynamics-based StochasticGameProblem has no multiplicative-noise
# representation yet, so that part of the math is not ported — it can be
# re-added once Base exposes it (see Base's docs/src/api/stochastic.md
# roadmap). Under additive-only noise, certainty equivalence ALWAYS holds:
# the feedback Nash gains are identical to the deterministic FNELQ gains,
# and the noise only shifts each player's expected cost by a computable
# constant, `Σ_k ½ tr(P_{k+1}^i W)`. SFNELQ therefore delegates the gain
# computation to FNELQ and adds this correction.
#
# References:
#   Başar & Olsder, Dynamic Noncooperative Game Theory, 2nd ed., SIAM 1999,
#   Ch. 6 (feedback Nash for discrete-time stochastic LQ games).
# ============================================================================

"""
    SFNELQ(; check_singularity=true, rcond_threshold=1e-10, regularization=0.0)

Stochastic feedback Nash equilibrium solver for additive-noise linear-
quadratic games (`StochasticGameProblem` with `is_certainty_equivalent ==
true`). Delegates gain computation to `FNELQ` (certainty equivalence) and
reports the expected cost correction due to the additive noise.

# Usage
```julia
sol = solve(sgame, SFNELQ())
sol.solver_info[:noise_cost]            # per-player expected cost from noise
sol.solver_info[:certainty_equivalent]  # always true in this round's scope
```
"""
struct SFNELQ <: GameSolver
    check_singularity::Bool
    rcond_threshold::Float64
    regularization::Float64
end

SFNELQ(;
    check_singularity::Bool = true,
    rcond_threshold::Real = 1e-10,
    regularization::Real = 0.0
) = SFNELQ(check_singularity, Float64(rcond_threshold), Float64(regularization))

get_covariance_trajectory(sol) = sol.solver_info[:state_covariances]
get_expected_costs(sol) = sol.solver_info[:expected_costs]
get_mean_trajectory(sol) = sol.state_trajectory

solver_capabilities(::Type{SFNELQ}) = [:LQGame, :Stochastic, :FeedbackPolicies, :DiscreteTime]

"""
    additive_noise_cost(sgame::StochasticGameProblem, Z_history) -> Vector

Per-player expected cost contributed by additive process noise,
`Σ_{k=1}^{N} ½ tr(P_{k+1}^i W)`, where `P^i` is the (FNELQ) value-function
Hessian sequence (`Z_history[i]`, length `N+1`) and `W = get_W(sgame)`.
"""
function additive_noise_cost(sgame::StochasticGameProblem{T}, Z_history::Vector{Vector{Matrix{T}}}) where {T}
    g = sgame.nominal_game
    np = num_players(g)
    N = n_steps(g)
    out = zeros(T, np)
    for i in 1:np
        # Uncertain initial state: E[½x₀ᵀZ₀x₀] = ½μ₀ᵀZ₀μ₀ + ½tr(Z₀Σ₀).
        out[i] += T(0.5) * tr(Z_history[i][1] * sgame.initial_covariance)
        for k in 1:N
            out[i] += T(0.5) * tr(Z_history[i][k+1] * get_W(sgame, k))
        end
        # FNELQ scales each player's deterministic cost by the objective scaling;
        # the noise contribution must be scaled identically to stay an expected cost.
        out[i] *= get_objective(g, i).scaling
    end
    return out
end

"""
    gaussian_trajectory_moments(sgame, strategy) -> NamedTuple

Mean-independent second moments of the physical state and joint control under
the affine feedback `u_k = -P_k x_k - α_k`:

- `state_covariances[k]`        = Cov(x_k),  k = 1,…,N+1 (`Σ₀ = initial_covariance`)
- `control_covariances[k]`      = Cov(u_k) = P Σ Pᵀ,  k = 1,…,N
- `state_control_covariances[k]`= Cov(x_k, u_k) = -Σ Pᵀ

`P_k` is the joint gain stacked over players. These are covariances of the
*physical* state under perfect state observation.
"""
function gaussian_trajectory_moments(sgame::StochasticGameProblem{T}, strategy) where {T}
    g = sgame.nominal_game
    dyn = g.dynamics
    np = num_players(g)
    N = n_steps(g)
    Σ = copy(sgame.initial_covariance)
    Σs = [copy(Σ)]; Σu = Matrix{T}[]; Σxu = Matrix{T}[]
    for k in 1:N
        P = vcat((strategy.gains[i][k] for i in 1:np)...)
        F = get_A(dyn, k) - get_B_concatenated(dyn, k) * P
        push!(Σu, Symmetric(P * Σ * P') |> Matrix)
        push!(Σxu, -Σ * P')
        Σ = Symmetric(F * Σ * F' + get_W(sgame, k)) |> Matrix
        push!(Σs, copy(Σ))
    end
    return (state_covariances = Σs, control_covariances = Σu, state_control_covariances = Σxu)
end

"""
    shared_observation_moments(sgame, strategy) -> NamedTuple

Moments of the physical state, the common estimate, and the joint control for
`u_k = -P_k x̂_k - α_k` under `SharedObservation` (see
[`common_information_moments`](@ref)). The physical-state covariance
`state_covariances` differs from the filter's `posterior_covariances`.
"""
function shared_observation_moments(sgame::StochasticGameProblem{T}, strategy) where {T}
    g = sgame.nominal_game; dyn = g.dynamics
    np = num_players(g); N = n_steps(g)
    obs = sgame.information
    A = [Matrix{T}(get_A(dyn, k)) for k in 1:N]
    B = [Matrix{T}(get_B_concatenated(dyn, k)) for k in 1:N]
    P = [vcat((strategy.gains[i][k] for i in 1:np)...) for k in 1:N]
    W = [Matrix{T}(get_W(sgame, k)) for k in 1:N]
    n = size(A[1], 1)
    H = [observation_matrix(obs, zeros(n), nothing, k) for k in 1:N+1]
    return common_information_moments(A, B, P, W, H, obs.noise_cov, sgame.initial_covariance)
end

"""
    expected_costs_from_moments(sgame, det_sol, Σx, Σu, Σxu) -> Vector

Exact expected player costs of the linear policy: the deterministic cost along
the mean trajectory plus the covariance terms
`½tr(QΣx) + ½tr(RΣu) + ⟨M, Cov(x,u)⟩` per stage and `½tr(Qf Σx_{N+1})`,
scaled by each objective's `scaling`. `Σx` has `N+1` entries, `Σu`/`Σxu` have `N`.
"""
function expected_costs_from_moments(sgame::StochasticGameProblem{T}, det_sol, Σx, Σu, Σxu) where {T}
    g = sgame.nominal_game; dyn = g.dynamics
    np = num_players(g); N = n_steps(g)
    n = total_state_dim(dyn); cd = dyn.control_dims; m = sum(cd)
    offs = [0; cumsum(cd)[1:end-1]]
    ranges = [offs[i]+1:offs[i]+cd[i] for i in 1:np]
    out = zeros(T, np)
    for i in 1:np
        obj = get_objective(g, i)
        extra = zero(T)
        for k in 1:N
            Q, R, M, _, _ = _joint_cost_mats(obj.stage_cost, i, k, ranges, n, m, T)
            extra += T(0.5) * tr(Q * Σx[k]) + T(0.5) * tr(R * Σu[k]) + sum(M .* Σxu[k])
        end
        extra += T(0.5) * tr(obj.terminal_cost.Qf * Σx[N+1])
        out[i] = det_sol.trajectories[i].cost + extra * obj.scaling
    end
    return out
end

function _solve(
    sgame::StochasticGameProblem{T},
    solver::SFNELQ,
    warmstart::Union{Nothing, GNEPSolution, WarmstartData},
    verbose::Bool
) where {T}
    is_certainty_equivalent(sgame) || error("""
        SFNELQ (this round's port) only supports additive (certainty-equivalent) noise.
        sgame.stochastic_dynamics has multiplicative (state/control-dependent) diffusion —
        this path is not yet implemented, see DifferentialGamesBase.jl docs/src/api/stochastic.md.
        """)

    info = sgame.information
    (info isa PerfectStateInformation || info isa SharedObservation) || error("""
        SFNELQ supports PerfectStateInformation and SharedObservation (common linear-Gaussian
        measurements) only; got $(typeof(info)). The separation of estimation and feedback-Nash
        control is not valid for private or asymmetric information.
        """)
    if info isa SharedObservation
        n = total_state_dim(sgame.nominal_game.dynamics)
        is_affine_observation(info, n) || error("""
            SFNELQ requires an affine common observation y = Hx + c + v (unit noise gain);
            a nonlinear observation model needs an EKF-based solver (e.g. CCiLQGames).
            """)
    end
    isempty(all_chance_constraints(sgame)) || error("""
        SFNELQ solves unconstrained games; the problem declares chance constraints.
        Use a chance-constrained solver (e.g. CCiLQGames) instead.
        """)

    g = sgame.nominal_game
    verbose && @info("SFNELQ: delegating gains to FNELQ (certainty equivalence)",
        players = num_players(g), steps = n_steps(g))

    det_sol = _solve(g, FNELQ(
        check_singularity = solver.check_singularity,
        rcond_threshold = solver.rcond_threshold,
        regularization = solver.regularization), warmstart, verbose)

    if !haskey(det_sol.solver_info, :cost_to_go_matrices)
        # FNELQ could not solve a stage game: propagate its failure status unchanged.
        return det_sol
    end
    Z_history = det_sol.solver_info[:cost_to_go_matrices]
    noise_cost = additive_noise_cost(sgame, Z_history)

    trajectories = [
        Trajectory(t.player_id, t.states, t.controls, t.times, t.cost + noise_cost[t.player_id])
        for t in det_sol.trajectories
    ]

    solver_info = copy(det_sol.solver_info)
    solver_info[:noise_cost] = noise_cost
    solver_info[:certainty_equivalent] = true
    solver_info[:expected_costs] = [t.cost for t in trajectories]
    if info isa SharedObservation
        # Separation: gains are the FNELQ gains applied to the common estimate x̂. The expected
        # cost is recomputed exactly from the joint law of (x, x̂) rather than from the
        # perfect-information value-function correction.
        mom = shared_observation_moments(sgame, det_sol.strategy)
        ec = expected_costs_from_moments(sgame, det_sol, mom.state_covariances,
                                         mom.control_covariances, mom.state_control_covariances)
        trajectories = [Trajectory(t.player_id, t.states, t.controls, t.times, ec[t.player_id])
                        for t in det_sol.trajectories]
        solver_info[:noise_cost] = ec .- [t.cost for t in det_sol.trajectories]
        solver_info[:expected_costs] = ec
        solver_info[:information] = :shared_observation
        solver_info[:separation_theorem_used] = true
        merge!(solver_info, Dict{Symbol, Any}(pairs(mom)))
    else
        merge!(solver_info, Dict{Symbol, Any}(pairs(gaussian_trajectory_moments(sgame, det_sol.strategy))))
    end

    return GNEPSolution(
        det_sol.game, trajectories;
        state_trajectory = det_sol.state_trajectory,
        strategy = det_sol.strategy,
        equilibrium_type = det_sol.equilibrium_type,
        converged = det_sol.converged,
        iterations = det_sol.iterations,
        solve_time = det_sol.solve_time,
        solver_info = solver_info
    )
end
