# ============================================================================
# DGSQP — Dynamic Game Sequential Quadratic Programming
#
# Julia implementation of the solver of
#
#   E. L. Zhu and F. Borrelli, "A Sequential Quadratic Programming Approach to
#   the Solution of Open-Loop Generalized Nash Equilibria", ICRA 2023
#   (arXiv:2203.16478), reference implementation https://github.com/zhu-edward/DGSQP
#   (MIT licence, Copyright (c) 2022 Edward Zhu).
#
# The method finds a local open-loop generalized Nash equilibrium (GNE) of a
# discrete-time game with nonlinear dynamics and constraints by solving the KKT
# system  ∇_{uⁱ}𝓛ⁱ(u, λ) = 0,  C(u) ≤ 0,  λ ≥ 0,  λᵀC(u) = 0  (shared multiplier λ,
# i.e. a normalized GNE) with a *single* convex QP per iteration:
#
#     min_p  ½ pᵀ B p + hᵀ p      s.t.  C + G p ≤ 0
#
# where h = [∇_{uⁱ}Jⁱ]ᵢ, G = ∂C/∂u and B is the symmetric positive semidefinite
# projection of the (non-symmetric) game Jacobian L = ∂[∇_{uⁱ}𝓛ⁱ]ᵢ/∂u. Steps are
# accepted with a merit function built on the KKT residual and a non-monotone
# (watchdog) line search.
#
# This file follows the reference implementation (DGSQP/solvers/DGSQP.py) line by
# line where the paper leaves a choice open: the merit function (`stat_l1`), the
# choice of μ, the monotone and watchdog line searches, the stopping tests and the
# least-squares dual initialisation. Differences are listed in DGSQP.md.
# ============================================================================

import Clarabel
using SparseArrays

# ============================================================================
# Solver options
# ============================================================================

"""
    DGSQP(; kwargs...) <: GameSolver

Dynamic Game SQP solver for open-loop generalized Nash equilibria (Zhu & Borrelli, ICRA 2023).
Defaults are those of the reference implementation's `DGSQPParams`;
[`dgsqp_paper_options`](@ref) returns the settings used for the paper's experiments.

# Options
- `sqp_iters = 50`           : maximum SQP iterations
- `p_tol = 1e-3`, `d_tol = 1e-3` : primal feasibility tolerance and stationarity/complementarity tolerance;
  converged when `max(0, maxᵢCᵢ) < p_tol`, `‖λ∘C‖∞ < d_tol` and `‖∇𝓛‖∞ < d_tol`
- `reg = 1e-3`               : regularization `ε` added to the projected Hessian, `B = proj_{⪰0}((L+Lᵀ)/2) + εI`
- `beta = 0.25`, `tau = 0.5` : sufficient-decrease constant and backtracking factor of the line search
- `line_search_iters = 50`   : backtracking steps
- `nonmono_ls = false`       : use the non-monotone (watchdog) line search of Sec. IV-B
- `watchdog_iters = 5`       : relaxed steps `t̂` the watchdog allows before insisting on a decrease
- `merit_function = :stat_l1`: `:stat_l1` (the paper's merit function, Eq. 15) or `:stat` (no ℓ¹ term; `μ = 0`)
- `rel_tol_iters = 3`        : consecutive iterations of tiny steps (`‖Δu‖ < p_tol/2`, `‖Δλ‖ < d_tol/2`) that
  terminate the solve when the iterate is primal feasible ("relative-tolerance" convergence)
- `divergence_tol = 1e5`     : declare divergence when `‖∇𝓛‖∞` exceeds this
- `merit_max = 1e6`          : the watchdog abandons relaxed steps when the merit exceeds this
- `time_limit = Inf`         : wall-clock limit in seconds
- `terminal_constraints = false` : for a `GameProblem`, also enforce state-only constraint rows at the terminal state
- `divergence_recovery = false` : when the stationarity residual blows up (`> divergence_tol`) after a relaxed step, return to the last
  checkpoint with a ten times larger regularization and a ten times smaller step limit instead of stopping (at most `max_recoveries` times).
  An extension of the reference, which stops with status `:diverged`; only for `step_acceptance = :nms`.
- `elastic_mode = false`      : when a QP is infeasible, solve an *elastic* QP instead, `min ½pᵀBp + hᵀp + ½ρ‖t‖² + η𝟙ᵀt` s.t. `C + Gp ≤ t`, `t ≥ 0`
  (`elastic_rho`, `elastic_eta`). The reference code defines this QP but its solve loop never uses it; here it removes `:qp_fail`
  stops caused by infeasible linearized constraints.
- `step_acceptance = :line_search` : `:nms` selects the step-acceptance scheme of the racing paper
  ([`dgsqp_racing_options`](@ref)): a decaying regularization `ε_q = reg·reg_decay^q`, "d-steps" (relaxed full steps while
  the step norm stays below a decaying limit `δ`) and "m-steps" every `nms_frequency` iterations that must decrease the
  merit function below `(1−merit_decrease)·max(memory)` of the last `nms_memory_size` values, with a return to the last
  checkpoint and an Armijo line search otherwise

The returned `GNEPSolution` has `equilibrium_type = :OpenLoopNash`, `solver_info[:status]` ∈
`:conv_abs_tol`, `:conv_rel_tol`, `:max_it`, `:diverged`, `:qp_fail`, `:time_limit`, and
the KKT residuals in `solver_info[:cond]`.
"""
Base.@kwdef struct DGSQP <: GameSolver
    sqp_iters            ::Int     = 50
    p_tol                ::Float64 = 1e-3
    d_tol                ::Float64 = 1e-3
    reg                  ::Float64 = 1e-3
    beta                 ::Float64 = 0.25
    tau                  ::Float64 = 0.5
    line_search_iters    ::Int     = 50
    nonmono_ls           ::Bool    = false
    watchdog_iters       ::Int     = 5
    merit_function       ::Symbol  = :stat_l1
    rel_tol_iters        ::Int     = 3
    divergence_tol       ::Float64 = 1e5
    merit_max            ::Float64 = 1e6
    time_limit           ::Float64 = Inf
    terminal_constraints ::Bool    = false
    # ── step acceptance of the racing paper (arXiv:2404.00186; `DGSQP_v2_frenet_approx.py`) ──
    step_acceptance      ::Symbol  = :line_search   # :line_search (ICRA 2023 solver) or :nms
    reg_decay            ::Float64 = 1.0            # ε_q = reg · reg_decay^q
    nms_frequency        ::Int     = 5
    nms_memory_size      ::Int     = 3
    merit_decrease       ::Float64 = 0.01
    delta_decay          ::Float64 = 0.95
    delta0_factor        ::Float64 = 20.0
    # ── extension (not in the reference solve loop): elastic QP when the QP is infeasible ──
    divergence_recovery  ::Bool    = false
    max_recoveries       ::Int     = 20
    elastic_mode         ::Bool    = false
    elastic_rho          ::Float64 = 1e3
    elastic_eta          ::Float64 = 1e3
end

"""
    dgsqp_paper_options(; kwargs...) -> DGSQP

The settings used for the curved-track experiment of the paper
(`scripts/DGSQP_ALGAMES_monte_carlo_curve.py` in the reference repository):
`reg = 0`, `nonmono_ls = true`, `beta = 0.01`, `tau = 0.5`, 50 SQP iterations, tolerances `1e-3`.
Keyword arguments override.
"""
dgsqp_paper_options(; kwargs...) =
    DGSQP(; reg = 0.0, nonmono_ls = true, beta = 0.01, tau = 0.5, sqp_iters = 50,
            p_tol = 1e-3, d_tol = 1e-3, line_search_iters = 50, kwargs...)

"""
    dgsqp_racing_options(; kwargs...) -> DGSQP

The settings of the racing paper (Zhu & Borrelli, arXiv:2404.00186; `DGSQPV2Params` of the reference code):
non-monotone step acceptance (`step_acceptance = :nms`), regularization `ε₀ = 100` decaying by `0.95` per
iteration, tolerances `1e-4`, up to 500 iterations, merit decrease `0.01`, checkpoint every 5 iterations, memory 3.
"""
dgsqp_racing_options(; kwargs...) =
    DGSQP(; step_acceptance = :nms, reg = 1e2, reg_decay = 0.95, p_tol = 1e-4, d_tol = 1e-4, sqp_iters = 500,
            line_search_iters = 50, tau = 0.5, merit_decrease = 0.01, nms_frequency = 5, nms_memory_size = 3,
            delta_decay = 0.95, merit_function = :stat_l1, divergence_tol = 1e10, kwargs...)

"""
    dgsqp_tuned_options(; kwargs...) -> DGSQP

The configuration tuned for the Austin racing benchmark of the FALCON paper (see the manual): the racing paper's non-monotone
step acceptance with a larger initial regularization (`reg = 1000`, decaying by `0.9`), a d-step limit decaying by `0.95`, checkpoints
every 10 iterations with a merit memory of 10, tolerances `p_tol = 1e-4`, `d_tol = 2e-2`, and the two extensions `elastic_mode` and
`divergence_recovery`. It converges on all 50 benchmark initial conditions (the paper's DG-SQP settings: 33 of 50), in about 390 iterations.
Keyword arguments override.
"""
dgsqp_tuned_options(; kwargs...) =
    DGSQP(; step_acceptance = :nms, nms_frequency = 10, nms_memory_size = 10, line_search_iters = 20, sqp_iters = 500,
            p_tol = 1e-4, d_tol = 2e-2, reg = 1e3, reg_decay = 0.9, delta_decay = 0.95, merit_decrease = 0.01, tau = 0.5,
            merit_function = :stat_l1, divergence_tol = 1e10, elastic_mode = true, divergence_recovery = true, kwargs...)

function solver_capabilities(::Type{DGSQP})
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
# Functional open-loop game
# ============================================================================

"""
    OpenLoopGame(control_dims, N, costs, constraints)

A finite-horizon open-loop game written directly in the players' control sequences. It is the
interface `DGSQP` works on; a `GameProblem` is converted to it by `solve`.

- `control_dims::Vector{Int}` : control dimension `mᵢ` of each player
- `N::Int`                    : horizon (number of control steps)
- `costs::Vector{Function}`   : `costs[i](w) -> Jⁱ`, the cost of player `i` as a function of the stacked decision vector
- `constraints::Function`     : `constraints(w) -> Vector` of all constraints `C(w) ≤ 0` (private and shared), or `nothing` for an unconstrained game

The decision vector `w = [vec(U¹); …; vec(U^M)]` stacks the players' controls, player by player, with
`Uⁱ ∈ ℝ^{mᵢ×N}` flattened column by column (time-major within a player). See [`control_ranges`](@ref),
[`split_controls`](@ref) and [`join_controls`](@ref). The functions must accept `ForwardDiff` dual numbers (they are
differentiated twice) — write them generically in the element type.

Because a player's cost is an arbitrary function of the *joint* decision vector, anything that depends on the
previous control (the input-rate costs and constraints of the paper) can be written directly in `w`.
"""
abstract type AbstractOpenLoopGame end

struct OpenLoopGame <: AbstractOpenLoopGame
    control_dims ::Vector{Int}
    N            ::Int
    costs        ::Vector{Function}
    constraints  ::Union{Nothing, Function}
end

"""
    control_ranges(game::OpenLoopGame) -> Vector{UnitRange{Int}}

Index range of each player's block of the decision vector.
"""
function control_ranges(g::AbstractOpenLoopGame)
    off = 0
    map(g.control_dims) do m
        r = off+1:off+m*g.N
        off += m * g.N
        r
    end
end

"""
    split_controls(game::OpenLoopGame, w) -> Vector{Matrix}

Per-player control matrices `Uⁱ` (`mᵢ × N`) of a decision vector.
"""
split_controls(g::AbstractOpenLoopGame, w::AbstractVector) =
    [reshape(w[r], g.control_dims[i], g.N) for (i, r) in enumerate(control_ranges(g))]

"""
    join_controls(U::Vector{<:AbstractMatrix}) -> Vector

Stack per-player control matrices into a decision vector (inverse of [`split_controls`](@ref)).
"""
join_controls(U::Vector{<:AbstractMatrix}) = reduce(vcat, vec.(U))

"""
    DGSQPResult

Result of [`dgsqp_solve`](@ref): `w` (decision vector), `λ` (multipliers, one per constraint row),
`converged`, `status`, `iterations`, `qp_solves`, `cond` (final `(p_feas, comp, stat)` residuals),
`history` (per-iteration residuals), `time`.
"""
struct DGSQPResult
    w          ::Vector{Float64}
    λ          ::Vector{Float64}
    converged  ::Bool
    status     ::Symbol
    iterations ::Int
    qp_solves  ::Int
    cond       ::NamedTuple{(:p_feas, :comp, :stat), NTuple{3, Float64}}
    history    ::Vector{NamedTuple{(:p_feas, :comp, :stat, :qp_solves), NTuple{4, Float64}}}
    time       ::Float64
end

# ============================================================================
# Derivatives
# ============================================================================

# Replace the block `rng` of `w` by `wi`, keeping the most general element type so
# that nested ForwardDiff duals propagate.
function _embed(w::AbstractVector, rng::UnitRange{Int}, wi::AbstractVector)
    S = promote_type(eltype(w), eltype(wi))
    W = Vector{S}(undef, length(w))
    W .= w
    W[rng] .= wi
    return W
end

_block_gradient(f, w, rng) = ForwardDiff.gradient(wi -> f(_embed(w, rng, wi)), w[rng])

"""
Evaluate `(Q, h, G, g)` at `(w, λ)`: the game Jacobian `Q = L_q` (non-symmetric), the stacked own-block
cost gradients `h`, the constraint Jacobian `G` and the constraint values `g`. With `hessian = false`
returns `(h, G, g)`.
"""
function _evaluate(game::OpenLoopGame, w::Vector{Float64}, λ::Vector{Float64}; hessian::Bool = true)
    ranges = control_ranges(game)
    M = length(ranges)
    has_c = game.constraints !== nothing
    g = has_c ? Float64.(game.constraints(w)) : Float64[]
    G = if has_c && !isempty(g)
        ForwardDiff.jacobian(game.constraints, w)
    else
        zeros(0, length(w))
    end
    h = reduce(vcat, [_block_gradient(game.costs[i], w, ranges[i]) for i in 1:M])
    hessian || return h, G, g

    Lag(i) = W -> game.costs[i](W) + (has_c && !isempty(g) ? dot(λ, game.constraints(W)) : zero(eltype(W)))
    F = ww -> reduce(vcat, [_block_gradient(Lag(i), ww, ranges[i]) for i in 1:M])
    Q = ForwardDiff.jacobian(F, w)
    return Q, h, G, g
end


# ============================================================================
# Structured exact derivatives for separable rollouts
# ============================================================================

"""
    SeparableRolloutGame(control_dims, state_dims, N, x0, steps, cost, constraints; u_prev0)

An open-loop game in which every player `i` has its own state `xⁱ` evolving by a discrete map
`xⁱ_{k+1} = steps[i](xⁱ_k, uⁱ_k)` (separable dynamics, as in a PD-GNEP) and everything else couples the
players only *at a node*. The cost and constraint functions are written per node:

- `cost(i, k, X, U, Up) -> Real` for `k = 0, …, N` (the stage cost at `k < N` and the terminal cost at `k = N`)
- `constraints(k, X, U, Up) -> Vector` for `k = 0, …, N` (rows that do not apply at a node are simply left out), meaning `≤ 0`

where `X`, `U`, `Up` are vectors with one entry per player: the states `xⁱ_k`, the controls `uⁱ_k` (zeros at
`k = N`) and the previous controls `uⁱ_{k−1}` (`u_prev0` at `k = 0`). The functions must be generic in the number type.

Compared with the generic [`OpenLoopGame`](@ref) the derivatives are assembled exactly from per-step
Jacobians and a second-order adjoint recursion instead of differentiating the whole rollout twice, so
the cost of the game Jacobian grows like `N²` small operations rather than `N²` full rollouts. This is what makes
long-horizon problems (hundreds of decision variables with an expensive step map) tractable.
"""
struct SeparableRolloutGame <: AbstractOpenLoopGame
    control_dims ::Vector{Int}
    state_dims   ::Vector{Int}
    N            ::Int
    x0           ::Vector{Vector{Float64}}
    steps        ::Vector{Function}
    cost         ::Function
    constraints  ::Function
    u_prev0      ::Vector{Vector{Float64}}
end

SeparableRolloutGame(control_dims, state_dims, N, x0, steps, cost, constraints;
                     u_prev0 = [zeros(m) for m in control_dims]) =
    SeparableRolloutGame(collect(Int, control_dims), collect(Int, state_dims), N,
                         [Vector{Float64}(x) for x in x0], collect(Function, steps), cost, constraints,
                         [Vector{Float64}(u) for u in u_prev0])

# Split a feature vector z = [x¹;…;x^M; u¹;…;u^M; up¹;…;up^M] into per-player pieces.
function _split_features(g::SeparableRolloutGame, z::AbstractVector)
    M = length(g.control_dims)
    o = 0
    X = [z[(o += g.state_dims[i]) - g.state_dims[i] + 1:o] for i in 1:M]
    U = [z[(o += g.control_dims[i]) - g.control_dims[i] + 1:o] for i in 1:M]
    Up = [z[(o += g.control_dims[i]) - g.control_dims[i] + 1:o] for i in 1:M]
    return X, U, Up
end

"""
Per-step Jacobians of the players' maps, the states and the sensitivities `Sⁱ_k = ∂xⁱ_k/∂uⁱ`.
"""
function _rollout_sensitivities(g::SeparableRolloutGame, w::AbstractVector; second_order::Bool = false)
    M, N = length(g.control_dims), g.N
    ranges = control_ranges(g)
    X  = [Vector{Vector{Float64}}(undef, N + 1) for _ in 1:M]
    As = [Vector{Matrix{Float64}}(undef, N) for _ in 1:M]
    S  = [Vector{Matrix{Float64}}(undef, N + 1) for _ in 1:M]
    for i in 1:M
        n, m = g.state_dims[i], g.control_dims[i]
        Ui = reshape(w[ranges[i]], m, N)
        X[i][1] = g.x0[i]
        S[i][1] = zeros(n, m * N)
        for k in 0:N-1
            x, u = X[i][k+1], Ui[:, k+1]
            y = vcat(x, u)
            J = ForwardDiff.jacobian(yy -> g.steps[i](yy[1:n], yy[n+1:end]), y)
            As[i][k+1] = J[:, 1:n]
            X[i][k+2] = g.steps[i](x, u)
            Sn = As[i][k+1] * S[i][k+1]
            Sn[:, k*m+1:(k+1)*m] .+= J[:, n+1:end]
            S[i][k+2] = Sn
        end
    end
    return X, As, S
end

function _evaluate(game::SeparableRolloutGame, w::Vector{Float64}, λ::Vector{Float64}; hessian::Bool = true)
    M, N = length(game.control_dims), game.N
    ranges = control_ranges(game)
    nw = length(w)
    nxs = sum(game.state_dims); mt = sum(game.control_dims); nz = nxs + 2mt
    xoff = cumsum([0; game.state_dims]); uoff = cumsum([0; game.control_dims])
    Uc = [reshape(w[ranges[i]], game.control_dims[i], N) for i in 1:M]
    X, As, S = _rollout_sensitivities(game, w)

    # Feature vectors and their Jacobians with respect to w at every node
    zs  = Vector{Vector{Float64}}(undef, N + 1)
    Jz  = Vector{Matrix{Float64}}(undef, N + 1)
    for k in 0:N
        z = Float64[]
        for i in 1:M; append!(z, X[i][k+1]); end
        for i in 1:M; append!(z, k < N ? Uc[i][:, k+1] : zeros(game.control_dims[i])); end
        for i in 1:M; append!(z, k > 0 ? Uc[i][:, k] : game.u_prev0[i]); end
        zs[k+1] = z
        J = zeros(nz, nw)
        for i in 1:M
            J[xoff[i]+1:xoff[i+1], ranges[i]] .= S[i][k+1]
            m = game.control_dims[i]
            if k < N
                for r in 1:m; J[nxs + uoff[i] + r, ranges[i][k*m + r]] = 1.0; end
            end
            if k > 0
                for r in 1:m; J[nxs + mt + uoff[i] + r, ranges[i][(k-1)*m + r]] = 1.0; end
            end
        end
        Jz[k+1] = J
    end

    cfun(k) = zz -> (t = _split_features(game, zz); game.constraints(k, t[1], t[2], t[3]))
    costfun(i, k) = zz -> (t = _split_features(game, zz); game.cost(i, k, t[1], t[2], t[3]))

    # Constraints
    g = Float64[]; Grows = Matrix{Float64}[]; Gz = Vector{Matrix{Float64}}(undef, N + 1)
    rows = Vector{UnitRange{Int}}(undef, N + 1)
    for k in 0:N
        ck = Float64.(cfun(k)(zs[k+1]))
        rows[k+1] = length(g)+1:length(g)+length(ck)
        append!(g, ck)
        if isempty(ck)
            Gz[k+1] = zeros(0, nz)
        else
            Gz[k+1] = ForwardDiff.jacobian(cfun(k), zs[k+1])
            push!(Grows, Gz[k+1] * Jz[k+1])
        end
    end
    G = isempty(Grows) ? zeros(0, nw) : reduce(vcat, Grows)

    # Cost gradients h
    h = zeros(nw)
    gz = [Vector{Vector{Float64}}(undef, N + 1) for _ in 1:M]
    for i in 1:M, k in 0:N
        gzi = ForwardDiff.gradient(costfun(i, k), zs[k+1])
        gz[i][k+1] = gzi
        h[ranges[i]] .+= Jz[k+1][:, ranges[i]]' * gzi
    end
    hessian || return h, G, g

    Q = zeros(nw, nw)
    for i in 1:M
        n, m = game.state_dims[i], game.control_dims[i]
        d = Vector{Vector{Float64}}(undef, N + 1)            # direct ∂𝓛ⁱ/∂xⁱ_k
        for k in 0:N
            λk = isempty(rows[k+1]) ? Float64[] : λ[rows[k+1]]
            Lf = zz -> costfun(i, k)(zz) + (isempty(λk) ? zero(eltype(zz)) : dot(λk, cfun(k)(zz)))
            Hz = ForwardDiff.hessian(Lf, zs[k+1])
            Q[ranges[i], :] .+= Jz[k+1][:, ranges[i]]' * Hz * Jz[k+1]
            gl = gz[i][k+1] .+ (isempty(λk) ? 0.0 : Gz[k+1]' * λk)
            d[k+1] = gl[xoff[i]+1:xoff[i+1]]
        end
        # second-order dynamics term: Σ_k T_kᵀ (p_{k+1}·∇²f) T_k with the adjoint p
        p = d[N+1]
        for k in N-1:-1:0
            y = vcat(X[i][k+1], Uc[i][:, k+1])
            Hs = ForwardDiff.hessian(yy -> dot(p, game.steps[i](yy[1:n], yy[n+1:end])), y)
            T = zeros(n + m, m * N)
            T[1:n, :] .= S[i][k+1]
            for r in 1:m; T[n + r, k*m + r] = 1.0; end
            Q[ranges[i], ranges[i]] .+= T' * Hs * T
            p = d[k+1] .+ As[i][k+1]' * p
        end
    end
    return Q, h, G, g
end


"""
    OpenLoopGame(game::SeparableRolloutGame) -> OpenLoopGame

The same game with generic (whole-rollout) derivatives. Slow for long horizons; used to cross-check the structured
derivatives and for `best_response_gap`. Constraint rows are ordered node by node, as in the structured game.
"""
function OpenLoopGame(g::SeparableRolloutGame)
    M, N = length(g.control_dims), g.N
    ranges = control_ranges(g)
    function nodes(w)
        S = eltype(w)
        Uc = [reshape(w[ranges[i]], g.control_dims[i], N) for i in 1:M]
        Xs = [Vector{Vector{S}}(undef, N + 1) for _ in 1:M]
        for i in 1:M
            Xs[i][1] = S.(g.x0[i])
            for k in 1:N
                Xs[i][k+1] = g.steps[i](Xs[i][k], Uc[i][:, k])
            end
        end
        return Uc, Xs, S
    end
    function args(Uc, Xs, S, k)
        X = [Xs[i][k+1] for i in 1:M]
        U = [k < N ? Uc[i][:, k+1] : zeros(S, g.control_dims[i]) for i in 1:M]
        Up = [k > 0 ? Uc[i][:, k] : S.(g.u_prev0[i]) for i in 1:M]
        return X, U, Up
    end
    cost(i) = function (w)
        Uc, Xs, S = nodes(w)
        sum(g.cost(i, k, args(Uc, Xs, S, k)...) for k in 0:N)
    end
    cons = function (w)
        Uc, Xs, S = nodes(w)
        reduce(vcat, [Vector{S}(g.constraints(k, args(Uc, Xs, S, k)...)) for k in 0:N])
    end
    return OpenLoopGame(g.control_dims, N, Function[cost(i) for i in 1:M], cons)
end

# ============================================================================
# QP step
# ============================================================================

"""
    _nearest_psd(A) -> Matrix

Symmetric positive semidefinite projection used for `B_q`: the symmetric part of `A` with its
negative eigenvalues replaced by `1e-10` (as in the reference implementation).
"""
function _nearest_psd(A::AbstractMatrix)
    B = (A + A') / 2
    s, U = eigen(Symmetric(B))
    s = max.(s, 1e-10)
    C = U * Diagonal(s) * U'
    return (C + C') / 2
end

# Elastic QP: min ½pᵀBp + hᵀp + ½ρ‖t‖² + η𝟙ᵀt  s.t.  g + G p ≤ t, t ≥ 0. Returns (p, λ̂) (multipliers of the first block) or `nothing`.
function _solve_elastic_qp(B::Matrix{Float64}, h::Vector{Float64}, G::Matrix{Float64}, g::Vector{Float64}, ρ::Float64, η::Float64)
    n, nc = length(h), length(g)
    nc == 0 && return nothing
    try
        P = sparse(triu(cat(B, Matrix(ρ * I, nc, nc); dims = (1, 2))))
        q = vcat(h, fill(η, nc))
        A = vcat(hcat(sparse(G), -sparse(1.0I, nc, nc)), hcat(spzeros(nc, n), -sparse(1.0I, nc, nc)))
        b = vcat(-g, zeros(nc))
        solver = Clarabel.Solver()
        Clarabel.setup!(solver, P, q, A, b, [Clarabel.NonnegativeConeT(2nc)], Clarabel.Settings(verbose = false))
        r = Clarabel.solve!(solver)
        (r.status == Clarabel.SOLVED || r.status == Clarabel.ALMOST_SOLVED) || return nothing
        (all(isfinite, r.x) && all(isfinite, r.z)) || return nothing
        return Vector{Float64}(r.x[1:n]), Vector{Float64}(r.z[1:nc])
    catch
        return nothing
    end
end

# Solve  min ½ pᵀBp + hᵀp  s.t.  g + G p ≤ 0.  Returns (p, λ̂) or `nothing` when the QP fails.
function _solve_qp(B::Matrix{Float64}, h::Vector{Float64}, G::Matrix{Float64}, g::Vector{Float64})
    nc = length(g)
    if nc == 0
        p = try
            -(Symmetric(B) \ h)
        catch
            return nothing
        end
        return p, Float64[]
    end
    try
        P = sparse(triu(B))
        solver = Clarabel.Solver()
        settings = Clarabel.Settings(verbose = false)
        Clarabel.setup!(solver, P, h, sparse(G), -g, [Clarabel.NonnegativeConeT(nc)], settings)
        r = Clarabel.solve!(solver)
        (r.status == Clarabel.SOLVED || r.status == Clarabel.ALMOST_SOLVED) || return nothing
        (all(isfinite, r.x) && all(isfinite, r.z)) || return nothing
        return Vector{Float64}(r.x), Vector{Float64}(r.z)
    catch
        return nothing
    end
end

# QP step with the optional elastic fallback of `solver`.
function _qp_step(B::Matrix{Float64}, h::Vector{Float64}, G::Matrix{Float64}, g::Vector{Float64}, solver)
    qp = _solve_qp(B, h, G, g)
    qp === nothing && solver.elastic_mode && (qp = _solve_elastic_qp(B, h, G, g, solver.elastic_rho, solver.elastic_eta))
    return qp
end

# ============================================================================
# Merit function (Eq. 15, reference implementation's `stat_l1`)
# ============================================================================

# φ = ½(‖h + Gᵀλ‖² + (λᵀg)²) + μ Σ(g − s)
_merit(l, s, h, G, g, μ) =
    0.5 * (sum(abs2, h .+ G' * l) + dot(l, g)^2) + μ * sum(g .- s)

# Directional derivative of φ along (du, dl) with s moving along ds.
function _dmerit(du, l, dl, Q, h, G, g, s, μ)
    d = h .+ G' * l
    dstat = dot(d, Q * du .+ G' * dl) + dot(l, g) * (dot(l, G * du) + dot(dl, g))
    return dstat - μ * sum(g .- s)
end

# Penalty parameter μ (reference `_get_mu`, ρ = 0.5).
function _get_mu(du, l, dl, s, Q, h, G, g, solver::DGSQP)
    solver.merit_function === :stat && return 0.0
    vio = sum(g .- s)
    dstat = _dmerit(du, l, dl, Q, h, G, g, s, 0.0)
    ρ = 0.5
    if dstat < 0 && vio > 0
        return -dstat / ((1 - ρ) * vio)
    elseif dstat >= 0 && vio > 0
        return dstat / ((1 - ρ) * vio)
    end
    return 0.0
end

# ============================================================================
# Line searches
# ============================================================================

# Monotone backtracking line search on the merit function (reference `_line_search_3`).
# Returns the last trial point even when the search is exhausted, as the reference does.
function _line_search(game, solver, w, du, l, dl, s, ds, Q, h, G, g, μ)
    φ  = _merit(l, s, h, G, g, μ)
    dφ = _dmerit(du, l, dl, Q, h, G, g, s, μ)
    α = 1.0
    w_t, l_t, φ_t = w, l, φ
    for _ in 1:solver.line_search_iters
        w_t = w .+ α .* du
        l_t = l .+ α .* dl
        s_t = s .+ α .* ds
        h_t, G_t, g_t = _evaluate(game, w_t, l_t; hessian = false)
        φ_t = _merit(l_t, s_t, h_t, G_t, g_t, μ)
        φ_t <= φ + solver.beta * α * dφ && break
        α *= solver.tau
    end
    return w_t, l_t, φ_t
end

# Compute the QP step at (w, l) from freshly evaluated derivatives.
# Returns (w+..., stuff) as a named tuple or `nothing` when the QP fails.
function _step_at(game, solver, w, l)
    Q, h, G, g = _evaluate(game, w, l)
    B = _nearest_psd(Q) + solver.reg * I
    qp = _qp_step(Matrix{Float64}(B), h, G, g, solver)
    qp === nothing && return nothing
    du, l_hat = qp
    s  = min.(0.0, g)
    ds = g .+ G * du .- s
    return (; Q, h, G, g, du, l_hat, dl = l_hat .- l, s, ds)
end

# Non-monotone (watchdog) line search (reference `_watchdog_line_search_4`).
# `st` is the step computed at the current iterate. Returns (w, l, n_qp).
function _watchdog(game, solver, w_k, l_k, st, μ)
    n_qp = 0
    φ_k  = _merit(l_k, st.s, st.h, st.G, st.g, μ)
    dφ_k = _dmerit(st.du, l_k, st.dl, st.Q, st.h, st.G, st.g, st.s, μ)

    # Relaxed (full) step
    w_1 = w_k .+ st.du
    l_1 = l_k .+ st.dl
    s_1 = st.s .+ st.ds
    h_1, G_1, g_1 = _evaluate(game, w_1, l_1; hessian = false)
    φ_1 = _merit(l_1, s_1, h_1, G_1, g_1, μ)
    φ_1 <= φ_k + solver.beta * dφ_k && return w_1, l_1, n_qp

    fail = false
    w_t, l_t = w_1, l_1
    w_t1, l_t1, φ_t1 = w_1, l_1, φ_1
    t0 = time()
    for _ in 1:solver.watchdog_iters
        s_t = _step_at(game, solver, w_t, l_t)
        n_qp += 1
        if s_t === nothing
            fail = true
            break
        end
        w_t1 = w_t .+ s_t.du
        l_t1 = s_t.l_hat
        sN   = s_t.s .+ s_t.ds
        h_n, G_n, g_n = _evaluate(game, w_t1, l_t1; hessian = false)
        φ_t1 = _merit(l_t1, sN, h_n, G_n, g_n, μ)
        φ_t1 > solver.merit_max && break
        if φ_t1 <= φ_k + solver.beta * dφ_k
            return w_t1, l_t1, n_qp
        end
        w_t, l_t = w_t1, l_t1
        if time() - t0 > solver.time_limit
            fail = true
            break
        end
    end

    # Insist on a decrease of the merit function
    if !fail
        s_t = _step_at(game, solver, w_t, l_t)
        n_qp += 1
        if s_t === nothing
            fail = true
        else
            μ_t = μ
            w_t1, l_t1, φ_t1 = _line_search(game, solver, w_t, s_t.du, l_t, s_t.dl, s_t.s, s_t.ds,
                                            s_t.Q, s_t.h, s_t.G, s_t.g, μ_t)
            if φ_t1 <= φ_k + solver.beta * dφ_k
                return w_t1, l_t1, n_qp
            elseif φ_t1 > φ_k
                fail = true
            else
                s_n = _step_at(game, solver, w_t1, l_t1)
                if s_n === nothing
                    w_f, l_f, _ = _line_search(game, solver, w_k, st.du, l_k, st.dl, st.s, st.ds,
                                               st.Q, st.h, st.G, st.g, μ)
                    return w_f, l_f, n_qp
                end
                n_qp += 1
                w_2, l_2, _ = _line_search(game, solver, w_t1, s_n.du, l_t1, s_n.dl, s_n.s, s_n.ds,
                                           s_n.Q, s_n.h, s_n.G, s_n.g, μ)
                return w_2, l_2, n_qp
            end
        end
    end

    # No decrease found: go back to a line search along the step at iteration k
    w_f, l_f, _ = _line_search(game, solver, w_k, st.du, l_k, st.dl, st.s, st.ds,
                               st.Q, st.h, st.G, st.g, μ)
    return w_f, l_f, n_qp
end


# ============================================================================
# Racing-paper variant: decaying regularization and non-monotone step acceptance
# ============================================================================

_phi2(l, s, h, G, μ) = 0.5 * sum(abs2, h .+ G' * l) + μ * sum(s)

# Port of `DGSQP_v2_frenet_approx.py`. Differences from the ICRA solver: the merit function has no complementarity
# term, the slack is the constraint violation `max(0, g)`, the regularization decays, and steps are accepted by the
# d-step / m-step scheme with a merit memory and checkpoints. Convergence wins over the iteration limit when both
# hold on the same iteration (the reference lets the later test overwrite the status). The "elastic mode" QP of the
# reference, used when the QP is infeasible, is not implemented.
function _dgsqp_nms(game::AbstractOpenLoopGame, w0, λ0, solver::DGSQP, verbose::Bool)
    t_start = time()
    w = Vector{Float64}(w0)
    h0, G0, g0 = _evaluate(game, w, Float64[]; hessian = false)
    nc = length(g0)
    l = λ0 !== nothing ? Vector{Float64}(λ0) :
        (nc == 0 ? Float64[] : max.(0.0, -(pinv(G0 * G0') * (G0 * h0))))
    w_prev, l_prev = copy(w), copy(l)
    memory = [_phi2(l, max.(0.0, g0), h0, G0, 1.0)]
    reg = solver.reg
    ck_counter, ck_index, ck_delta, ck_reg = 0, 0, 0.0, reg
    delta = 0.0
    data = NamedTuple[]
    history = NamedTuple{(:p_feas, :comp, :stat, :qp_solves), NTuple{4, Float64}}[]
    cond = (p_feas = Inf, comp = Inf, stat = Inf)
    converged = false; status = :max_it
    rel_its = 0; qp_solves = 0; it = 0; m_it = 0; recoveries = 0

    # merit decrease along a (possibly reloaded) step: Armijo backtracking, as the reference `line_search`
    function line_search(u, l_, du, dl, μ, phi_ref, dphi)
        α = 1.0
        u_t, l_t = u, l_
        s_t = Float64[]; h_t = Float64[]; G_t = zeros(0, length(u)); g_t = Float64[]
        for _ in 1:solver.line_search_iters
            u_t = u .+ α .* du; l_t = l_ .+ α .* dl
            h_t, G_t, g_t = _evaluate(game, u_t, l_t; hessian = false)
            s_t = max.(0.0, g_t)
            _phi2(l_t, s_t, h_t, G_t, μ) <= phi_ref + solver.merit_decrease * α * dphi && break
            α *= solver.tau
        end
        return u_t, l_t, α, _phi2(l_t, s_t, h_t, G_t, 1.0)
    end

    while true
        Q, h, G, g = _evaluate(game, w, l)
        d = h .+ G' * l
        p_feas = isempty(g) ? 0.0 : max(0.0, maximum(g))
        comp   = isempty(g) ? 0.0 : norm(g .* l, Inf)
        stat   = norm(d, Inf)
        cond = (p_feas = p_feas, comp = comp, stat = stat)
        verbose && @printf("DGSQP-NMS it %3d | p feas %.3e | comp %.3e | stat %.3e | reg %.2e\n", it, p_feas, comp, stat, reg)
        push!(history, (p_feas = p_feas, comp = comp, stat = stat, qp_solves = Float64(qp_solves)))

        done = false
        if (stat > solver.divergence_tol || !isfinite(stat)) && solver.divergence_recovery && recoveries < solver.max_recoveries && !isempty(data)
            ck = data[min(ck_index, length(data) - 1) + 1]
            w, l = copy(ck.u), copy(ck.l)
            w_prev, l_prev = copy(w), copy(l)
            reg = min(max(ck_reg, reg) * 10, 1e8)
            delta = max(ck_delta, 1e-8) * 0.1
            ck_delta, ck_reg = delta, reg
            ck_counter = 0
            recoveries += 1
            verbose && @printf("DGSQP-NMS divergence: restart from checkpoint %d with reg %.1e\n", ck_index, reg)
            continue
        end
        if stat > solver.divergence_tol
            status, done = :diverged, true
        end
        if m_it >= solver.sqp_iters
            status, done = :max_it, true
        end
        if time() - t_start > solver.time_limit
            status, done = :time_limit, true
        end
        if p_feas < solver.p_tol && comp < solver.d_tol && stat < solver.d_tol
            converged, status, done = true, :conv_abs_tol, true
        end
        done && break

        B = _nearest_psd(Q) + reg * I
        qp = _qp_step(Matrix{Float64}(B), h, G, g, solver)
        qp_solves += 1
        d_step = false; m_step = false
        u, lc = w, l
        du = dl = s = ds = Float64[]; μ = 0.0
        phi_ref = dphi = 0.0
        if qp === nothing
            if isempty(data)
                status = :qp_fail
                break
            end
            # QP failed: m-step from the last checkpoint
            ck = data[min(ck_index, length(data) - 1) + 1]
            u, lc, du, dl, s, ds, μ = ck.u, ck.l, ck.du, ck.dl, ck.s, ck.ds, ck.μ
            phi_ref, dphi = ck.phi_ref, ck.dphi
            m_step = true
            this = ck
        else
            du, l_hat = qp
            dl = l_hat .- l
            if it == 0
                delta = solver.delta0_factor * norm(vcat(du, dl))
                ck_delta = delta
            end
            if ck_counter >= solver.nms_frequency
                m_step = true
            elseif norm(vcat(du, dl)) < delta
                d_step = true
            else
                m_step = true
            end
            s  = max.(0.0, g)
            ds = max.(0.0, g .+ G * du) .- s
            dφc = dot(d, Q * du .+ G' * dl)
            μ = 0.0
            if solver.merit_function === :stat_l1
                ssum = sum(s)
                if ssum > 0
                    μ = dφc < 0 ? -dφc / (0.5 * ssum) : dφc / (0.5 * ssum)
                end
            end
            phi_ref = _phi2(l, s, h, G, μ)
            dphi = dφc - μ * sum(s)
            this = (u = copy(w), l = copy(l), du = du, dl = dl, s = s, ds = ds, μ = μ, phi_ref = phi_ref, dphi = dphi)
        end

        if d_step
            w = w .+ du; l = l .+ dl
            delta *= solver.delta_decay
            ck_counter += 1
        end
        if m_step
            m_it += 1
            u_t = u .+ du; l_t = lc .+ dl
            h_t, G_t, g_t = _evaluate(game, u_t, l_t; hessian = false)
            s_t = max.(0.0, g_t)
            φ = _phi2(l_t, s_t, h_t, G_t, μ)
            R = (1 - solver.merit_decrease) * maximum(memory)
            if φ <= R
                w, l = u_t, l_t
            else
                if ck_index <= length(data) - 1 && qp !== nothing
                    # watchdog: return to the last checkpoint and enforce a decrease from there
                    ck = data[ck_index + 1]
                    u, lc, du, dl, μ = ck.u, ck.l, ck.du, ck.dl, ck.μ
                    phi_ref, dphi = ck.phi_ref, ck.dphi
                    this = ck
                    delta = ck_delta
                    reg = ck_reg
                end
                w, l, _, φ = line_search(u, lc, du, dl, μ, phi_ref, dphi)
            end
            if norm(w .- w_prev) < solver.p_tol && norm(l .- l_prev) < solver.d_tol
                rel_its += 1
                if rel_its >= solver.rel_tol_iters && p_feas < solver.p_tol && comp < solver.d_tol && stat < solver.d_tol
                    converged, status = true, :conv_rel_tol
                    push!(data, this)
                    break
                end
            else
                rel_its = 0
            end
            w_prev, l_prev = copy(w), copy(l)
            reg *= solver.reg_decay
            push!(memory, φ); length(memory) > solver.nms_memory_size && popfirst!(memory)
            ck_counter = 0; ck_delta = delta; ck_reg = reg; ck_index = it + 1
        end
        push!(data, this)
        it += 1
        if m_it >= solver.sqp_iters + 1 && false
            break
        end
    end
    return DGSQPResult(w, l, converged, status, it, qp_solves, cond, history, time() - t_start)
end

# ============================================================================
# Main loop (Algorithm 1)
# ============================================================================

"""
    dgsqp_solve(game::OpenLoopGame, w0; solver = DGSQP(), λ0 = nothing, verbose = false) -> DGSQPResult

Run Algorithm 1 of Zhu & Borrelli on an [`OpenLoopGame`](@ref) from the initial control guess `w0`.
`λ0` defaults to the least-squares estimate `max(0, −(GGᵀ)⁺ G h)` of the paper (Alg. 1, line 1).
"""
function dgsqp_solve(game::AbstractOpenLoopGame, w0::AbstractVector;
                     solver::DGSQP = DGSQP(), λ0 = nothing, verbose::Bool = false)
    solver.step_acceptance === :nms && return _dgsqp_nms(game, w0, λ0, solver, verbose)
    t_start = time()
    w = Vector{Float64}(w0)
    h0, G0, g0 = _evaluate(game, w, Float64[]; hessian = false)
    nc = length(g0)
    l = if λ0 !== nothing
        Vector{Float64}(λ0)
    elseif nc == 0
        Float64[]
    else
        max.(0.0, -(pinv(G0 * G0') * (G0 * h0)))
    end

    history = NamedTuple{(:p_feas, :comp, :stat, :qp_solves), NTuple{4, Float64}}[]
    cond = (p_feas = Inf, comp = Inf, stat = Inf)
    converged = false
    status = :max_it
    rel_its = 0
    qp_solves = 0
    it = 0

    while true
        Q, h, G, g = _evaluate(game, w, l)
        d = h .+ G' * l
        p_feas = isempty(g) ? 0.0 : max(0.0, maximum(g))
        comp   = isempty(g) ? 0.0 : norm(g .* l, Inf)
        stat   = norm(d, Inf)
        cond = (p_feas = p_feas, comp = comp, stat = stat)
        verbose && @printf("DGSQP it %3d | p feas %.3e | comp %.3e | stat %.3e\n", it, p_feas, comp, stat)

        if stat > solver.divergence_tol
            status = :diverged
            push!(history, (p_feas = p_feas, comp = comp, stat = stat, qp_solves = Float64(qp_solves)))
            break
        end
        if p_feas < solver.p_tol && comp < solver.d_tol && stat < solver.d_tol
            converged, status = true, :conv_abs_tol
            push!(history, (p_feas = p_feas, comp = comp, stat = stat, qp_solves = Float64(qp_solves)))
            break
        end

        B = _nearest_psd(Q) + solver.reg * I
        qp = _qp_step(Matrix{Float64}(B), h, G, g, solver)
        qp_solves += 1
        if qp === nothing
            status = :qp_fail
            push!(history, (p_feas = p_feas, comp = comp, stat = stat, qp_solves = Float64(qp_solves)))
            break
        end
        du, l_hat = qp
        dl = l_hat .- l
        s  = min.(0.0, g)
        ds = g .+ G * du .- s
        st = (; Q, h, G, g, du, l_hat, dl, s, ds)
        μ  = _get_mu(du, l, dl, s, Q, h, G, g, solver)

        w_prev, l_prev = w, l
        if solver.nonmono_ls
            w, l, n_qp = _watchdog(game, solver, w, l, st, μ)
            qp_solves += n_qp
        else
            w, l, _ = _line_search(game, solver, w, du, l, dl, s, ds, Q, h, G, g, μ)
        end
        push!(history, (p_feas = p_feas, comp = comp, stat = stat, qp_solves = Float64(qp_solves)))

        # Relative-tolerance convergence
        if norm(w .- w_prev) < solver.p_tol / 2 && norm(l .- l_prev) < solver.d_tol / 2
            rel_its += 1
            if rel_its >= solver.rel_tol_iters && p_feas < solver.p_tol
                converged, status = true, :conv_rel_tol
                break
            end
        else
            rel_its = 0
        end

        it += 1
        if it >= solver.sqp_iters
            status = :max_it
            break
        end
        if time() - t_start > solver.time_limit
            status = :time_limit
            break
        end
    end

    return DGSQPResult(w, l, converged, status, it, qp_solves, cond, history, time() - t_start)
end

# ============================================================================
# GameProblem adapter
# ============================================================================

"""
    OpenLoopGame(game::GameProblem; terminal_constraints = false) -> OpenLoopGame

Express a Base `GameProblem` in the players' control sequences: the state trajectory is the
rollout of the discrete dynamics (RK4 for continuous dynamics, the exact map for `LinearDynamics`,
as in `ALGAMES`), the cost of player `i` is `total_cost` of the rollout, and the constraints are all
private and shared constraints evaluated at `(x_k, u_k)` for `k = 1, …, N` (with `x_1 = x₀`). Cost and
constraint functions follow the same argument conventions as `ALGAMES`.
"""
function OpenLoopGame(game::GameProblem{T}; terminal_constraints::Bool = false) where {T}
    np   = game.n_players
    N    = n_steps(game)
    cdim = [game.metadata.control_dims[i] for i in 1:np]
    coff = [game.metadata.control_offsets[i] for i in 1:np]
    mtot = sum(cdim)
    n    = total_state_dim(game.dynamics)
    dt   = game.time_horizon.dt
    x0   = Vector{T}(game.initial_state)
    separable = length(game.metadata.state_dims) == np

    # Rollout of the joint state from a decision vector (generic in the element type)
    function rollout(w)
        S = eltype(w)
        U = Matrix{S}(undef, mtot, N)
        off = 0
        for i in 1:np
            m = cdim[i]
            Ui = reshape(w[off+1:off+m*N], m, N)
            U[coff[i]+1:coff[i]+m, :] .= Ui
            off += m * N
        end
        X = Vector{Vector{S}}(undef, N + 1)
        X[1] = S.(x0)
        for k in 1:N
            X[k+1] = _discrete_step(game.dynamics, X[k], U[:, k], nothing, k, dt)
        end
        return X, U
    end

    cost_fn(i) = function (w)
        X, U = rollout(w)
        obj = get_objective(game, i)
        if obj.stage_cost isa LQStageCost
            Xl = separable ?
                [X[k][game.metadata.state_offsets[i]+1:game.metadata.state_offsets[i]+game.metadata.state_dims[i]] for k in 1:N+1] :
                X
            Ul = [U[coff[i]+1:coff[i]+cdim[i], k] for k in 1:N]
            return total_cost(obj, Xl, Ul, nothing)
        else
            return total_cost(obj, X, [U[:, k] for k in 1:N], nothing)
        end
    end

    cons = collect(Iterators.flatten((game.private_constraints, game.shared_constraints)))
    term_mask = (terminal_constraints && !isempty(cons)) ? _terminal_mask(game, n, mtot) : Bool[]
    constraints = isempty(cons) ? nothing : function (w)
        X, U = rollout(w)
        S = eltype(w)
        parts = Vector{S}[]
        for k in 1:N
            for c in cons
                push!(parts, Vector{S}(evaluate_constraint(c, X[k], U[:, k], nothing, k)))
            end
        end
        if !isempty(term_mask)
            xT = X[N+1]; uT = zeros(S, mtot)
            Cs = reduce(vcat, [Vector{S}(evaluate_constraint(c, xT, uT, nothing, N + 1)) for c in cons])
            push!(parts, Cs[term_mask])
        end
        return reduce(vcat, parts)
    end

    return OpenLoopGame(cdim, N, Function[cost_fn(i) for i in 1:np], constraints)
end

# ============================================================================
# _solve — Base interface
# ============================================================================

_solve(game::GameProblem{T}, solver::DGSQP, warmstart::GNEPSolution{T}, verbose::Bool) where {T} =
    _solve(game, solver, WarmstartData(warmstart), verbose)

function _solve(game::GameProblem{T}, solver::DGSQP,
                warmstart::Union{Nothing, WarmstartData}, verbose::Bool) where {T}
    olg = OpenLoopGame(game; terminal_constraints = solver.terminal_constraints)
    np, N = game.n_players, olg.N

    U0 = [zeros(Float64, olg.control_dims[i], N) for i in 1:np]
    if warmstart !== nothing && warmstart.trajectories !== nothing
        for i in 1:np
            idx = findfirst(t -> t.player_id == i, warmstart.trajectories)
            idx === nothing && continue
            Ui = warmstart.trajectories[idx].controls
            nu = min(size(Ui, 2), N)
            U0[i][:, 1:nu] .= Ui[:, 1:nu]
        end
    end
    λ0 = nothing
    if warmstart !== nothing && warmstart.dual_variables !== nothing &&
       haskey(warmstart.dual_variables, :λ_dgsqp)
        λ0 = warmstart.dual_variables[:λ_dgsqp]
    end

    res = dgsqp_solve(olg, join_controls(U0); solver = solver, λ0 = λ0, verbose = verbose)
    return _assemble_dgsqp_solution(game, olg, res)
end

function _assemble_dgsqp_solution(game::GameProblem{T}, olg::OpenLoopGame, res::DGSQPResult) where {T}
    np, N = game.n_players, olg.N
    cdim = olg.control_dims
    U = split_controls(olg, res.w)
    dt = game.time_horizon.dt
    n = total_state_dim(game.dynamics)
    X = Matrix{T}(undef, n, N + 1)
    X[:, 1] .= game.initial_state
    Uflat = vcat([Matrix{T}(u) for u in U]...)
    for k in 1:N
        X[:, k+1] .= _discrete_step(game.dynamics, X[:, k], Uflat[:, k], nothing, k, dt)
    end
    times = collect(range(T(0), game.time_horizon.tf, length = N + 1))
    trajectories = map(1:np) do i
        c_i = olg.costs[i](res.w)
        Trajectory{T}(i, copy(X), Matrix{T}(U[i]), times, T(c_i))
    end
    strategy = OpenLoopStrategy([Matrix{T}(u) for u in U], copy(cdim), times)
    info = Dict{Symbol, Any}(
        :status    => res.status,
        :cond      => res.cond,
        :qp_solves => res.qp_solves,
        :history   => res.history,
        :duals     => Dict{Symbol, Any}(:λ_dgsqp => copy(res.λ)),
        :dual_variables => Dict{Symbol, Any}(:λ_dgsqp => copy(res.λ)),
    )
    return GNEPSolution(game, trajectories;
        state_trajectory = X, strategy = strategy, equilibrium_type = :OpenLoopNash,
        converged = res.converged, iterations = res.iterations, solve_time = res.time,
        solver_info = info)
end

# ============================================================================
# Equilibrium check: unilateral deviations
# ============================================================================

"""
    best_response_gap(game::OpenLoopGame, w, i; solver) -> (gap, w_br)

Equilibrium certificate for player `i`: hold the other players' controls in `w` fixed, solve player
`i`'s own constrained problem (an open-loop best response, from `w`'s own block, with `DGSQP` itself as the
single-player SQP solver) and return `gap = Jⁱ(w) − Jⁱ(w_br) ≥ 0` together with the decision vector with the
improved block. At a local generalized Nash equilibrium no unilateral deviation lowers the cost, so every
player's gap is zero (up to the solver tolerance); a positive gap exposes a point that is not an equilibrium
even when the KKT residuals look small. Local: it finds the nearest local best response.
"""
function best_response_gap(game::OpenLoopGame, w::AbstractVector, i::Int;
                           solver::DGSQP = DGSQP(; nonmono_ls = true, reg = 0.0, beta = 0.01,
                                                   p_tol = 1e-7, d_tol = 1e-7, sqp_iters = 100))
    r = control_ranges(game)[i]
    w = Vector{Float64}(w)
    sub = OpenLoopGame([game.control_dims[i]], game.N,
                       Function[wi -> game.costs[i](_embed(w, r, wi))],
                       game.constraints === nothing ? nothing : (wi -> game.constraints(_embed(w, r, wi))))
    res = dgsqp_solve(sub, w[r]; solver = solver)
    w_br = copy(w); w_br[r] .= res.w
    # A best response must be feasible; if the inner solve did not manage, report no improvement.
    if !isempty(res.λ) && res.cond.p_feas > 1e-5
        return 0.0, w
    end
    return max(0.0, game.costs[i](w) - game.costs[i](w_br)), w_br
end

"""
    best_response_gap(game::GameProblem, sol::GNEPSolution; kwargs...) -> Vector{Float64}

Best-response gaps of every player of a solution of `game` (any solver), using the control sequences
in `sol`. See the `OpenLoopGame` method.
"""
function best_response_gap(game::GameProblem, sol::GNEPSolution; kwargs...)
    olg = OpenLoopGame(game)
    w = join_controls([Matrix{Float64}(get_trajectory(sol, i).controls) for i in 1:game.n_players])
    return [best_response_gap(olg, w, i; kwargs...)[1] for i in 1:game.n_players]
end
