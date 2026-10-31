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
    W = get_W(sgame)
    out = zeros(T, np)
    for i in 1:np, k in 1:N
        out[i] += T(0.5) * tr(Z_history[i][k+1] * W)
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

    g = sgame.nominal_game
    verbose && @info("SFNELQ: delegating gains to FNELQ (certainty equivalence)",
        players = num_players(g), steps = n_steps(g))

    det_sol = _solve(g, FNELQ(
        check_singularity = solver.check_singularity,
        rcond_threshold = solver.rcond_threshold,
        regularization = solver.regularization), warmstart, verbose)

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
