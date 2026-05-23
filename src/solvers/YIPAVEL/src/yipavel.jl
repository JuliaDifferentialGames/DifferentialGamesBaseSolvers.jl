# ============================================================================
# yipavel.jl — Distributed Primal-Dual GNE Solver (Yi & Pavel 2017)
#
# Reference:
#   Yi, P. & Pavel, L. (2017). "A distributed primal-dual algorithm for
#   computation of generalized Nash equilibria via operator splitting."
#   IEEE CDC 2017, pp. 3841–3846.
#
# ── Problem class ─────────────────────────────────────────────────────────────
#   Targets ConvexGameProblem{T}: an N-player convex game with separable
#   dynamics, convex private feasible sets Ωᵢ, and a shared affine coupling
#   constraint Ax ≤ b (supply ≤ capacity) or Ax ≥ b (paper form).
#
#   Each player i controls uᵢ ∈ Ωᵢ ⊂ ℝ^{mᵢ} (vectorized over N time steps)
#   and minimises fᵢ(u₁,…,uₙ), which may couple across players.
#
# ── Algorithm 1 (Yi & Pavel 2017, fully-connected communication graph) ────────
#   Initialise: x = 0, z = 0, λ = 0.
#
#   for k = 0, 1, 2, …:
#     # Primal (projected gradient descent on Lagrangian, simultaneous):
#     x_{i,k+1} = P_{Ωᵢ}[ x_{i,k} − τᵢ(∇_{xᵢ}fᵢ(x_k) − Aᵢᵀλ_{i,k}) ]
#
#     # Auxiliary (Laplacian-based consensus on λ):
#     z_{i,k+1} = z_{i,k} + νᵢ Σ_{j≠i}(λ_{i,k} − λ_{j,k})
#
#     # Dual (projected gradient ascent on dual, using z for consensus):
#     λ_{i,k+1} = P_{ℝ₊ᵐ}[ λ_{i,k} − σᵢ(
#                   Aᵢ(2x_{i,k+1} − x_{i,k}) − bᵢ
#                   + Σ_{j≠i}[2(z_{i,k+1}−z_{j,k+1}) − (z_{i,k}−z_{j,k})]
#                   + Σ_{j≠i}(λ_{i,k} − λ_{j,k}) ) ]
#
#   Coupling is stored internally in Ax ≥ b form.  The public API accepts
#   Ax ≤ b (the natural form for capacity constraints) and converts via sign flip.
#
# ── Namespace note ────────────────────────────────────────────────────────────
#   Included by DifferentialGamesBaseSolvers. All DGB types are available via
#   `using DifferentialGamesBase` in the parent module.
# ============================================================================

# ============================================================================
# YiPavel solver struct
# ============================================================================

"""
    YiPavel <: GameSolver

Distributed primal-dual operator-splitting GNE solver for `ConvexGameProblem{T}`.

Implements Algorithm 1 of Yi & Pavel (2017) for computing the *variational GNE*
(v-NE) of an N-player convex game with a globally shared affine coupling
constraint.

Each iteration simultaneously updates all players' primal decisions, auxiliary
consensus variables (z), and local copies of the coupling multiplier (λ):

- **Primal**: projected gradient descent on the local Lagrangian
  `Lᵢ = fᵢ + λᵢᵀ(b − Aᵢxᵢ)`.
- **Auxiliary (z)**: PI-like dynamics that accumulate the Laplacian of λ to
  enforce consensus across the fully-connected communication graph.
- **Dual**: projected ascent on λᵢ, using an extrapolated constraint residual
  and the consensus correction from z.

Convergence to the v-NE is guaranteed under Assumptions 1–3 of Yi & Pavel
(2017) when step-sizes satisfy Theorem 2 (Lemma 3).

# Fields
- `max_iter` : Maximum number of iterations (default 5000).
- `tol`      : Convergence tolerance ‖Δx‖_∞ < tol (default 1e-6).
- `τ`        : Primal step size τᵢ (same for all players; default 0.05).
- `ν`        : Auxiliary (z) step size νᵢ (default 0.02).
- `σ`        : Dual step size σᵢ (default 0.05).

# Usage
Pass cost functions and coupling data as keyword arguments to `solve`:
```julia
sol = solve(game, YiPavel();
            cost_fns    = [f1, f2, ...],   # f_i(seg_1,...,seg_N) -> Real
            coupling_A  = A,               # m × Σdⱼ matrix
            coupling_b  = b,               # m-vector, Ax ≤ b
            coupling_leq = true,           # true = Ax ≤ b (default)
            lb          = [lb_1, ...],     # per-player lower bounds
            ub          = [ub_1, ...])     # per-player upper bounds
```

# References
Yi, P. & Pavel, L. (2017). A distributed primal-dual algorithm for computation
of generalized Nash equilibria via operator splitting. *IEEE CDC*, pp. 3841–3846.
"""
@kwdef struct YiPavel <: GameSolver
    max_iter::Int     = 5000
    tol     ::Float64 = 1e-6
    τ       ::Float64 = 0.05   # primal step size (same for all players)
    ν       ::Float64 = 0.02   # auxiliary (z) step size
    σ       ::Float64 = 0.05   # dual step size
end

# ── Capabilities ──────────────────────────────────────────────────────────────

function solver_capabilities(::Type{YiPavel})
    return [:ConvexGame, :GeneralizedNash, :VariationalNash,
            :SharedConstraints, :SeparableDynamics, :DiscreteTime]
end

# ============================================================================
# Internal helpers
# ============================================================================

# Box projection
_yp_project_box(x, lb, ub) = min.(ub, max.(lb, x))

# Non-negative orthant projection
_yp_project_nonneg(x) = max.(zero(eltype(x)), x)

"""
    _yp_laplacian_row(i, v) -> Vector

Row i of the graph Laplacian times `v`, for a fully-connected unit-weight graph:
  (Lv)ᵢ = (N−1)vᵢ − Σ_{j≠i} vⱼ

For a single player (N=1), returns zero (no consensus needed).
"""
function _yp_laplacian_row(i::Int, v::Vector{<:AbstractVector})
    N = length(v)
    N == 1 && return zero(v[i])
    s = zero(v[1])
    for j in 1:N
        j == i && continue
        s = s .+ v[j]
    end
    return (N - 1) .* v[i] .- s
end

"""
    _yp_rollout(f_i, x0, U, dt) -> Matrix

Euler rollout for player i's dynamics.  ForwardDiff-compatible: element type is
promoted so dual numbers propagate through the rollout and into cost functions.
"""
function _yp_rollout(f_i, x0::AbstractVector{S}, U::AbstractMatrix, dt::Real) where {S}
    n     = length(x0)
    N     = size(U, 2)
    TE    = promote_type(S, eltype(U))
    X     = Matrix{TE}(undef, n, N + 1)
    X[:, 1] .= x0
    for k in 1:N
        X[:, k+1] .= X[:, k] .+ dt .* f_i(X[:, k], U[:, k], nothing, (k-1)*dt)
    end
    return X
end

# ============================================================================
# Gradient computation
# ============================================================================

"""
    _yp_gradient_coupled(cost_fn_i, x_all, i, d_i) -> Vector

Compute ∇_{xᵢ} fᵢ(x₁,…,xₙ) via ForwardDiff over the joint decision vector.
`cost_fn_i(seg_1,…,seg_N)` is player i's objective as a function of all
players' (vectorised) control trajectories.
"""
function _yp_gradient_coupled(
    cost_fn_i::Function,
    x_all    ::Vector{<:AbstractVector{T}},
    i        ::Int,
    d_i      ::Vector{Int}
) where {T}
    offsets = [0; cumsum(d_i)]
    np      = length(d_i)
    x_joint = vcat(x_all...)
    g = ForwardDiff.gradient(x_joint) do v
        segs = [v[offsets[j]+1:offsets[j+1]] for j in 1:np]
        cost_fn_i(segs...)
    end
    return g[offsets[i]+1:offsets[i+1]]
end

"""
    _yp_gradient_lq(fwd, i, u_flat, N, dt) -> Vector

Compute ∇_{vec(Uᵢ)} Jᵢ(Xᵢ,Uᵢ) for a player with an LQ objective,
using ForwardDiff through the Euler rollout.  Works for any (smooth) dynamics.
"""
function _yp_gradient_lq(
    fwd   ::GameProblem{T},
    i     ::Int,
    u_flat::AbstractVector,
    N     ::Int,
    dt    ::T
) where {T}
    obj   = get_objective(fwd, i)
    m_i   = fwd.metadata.control_dims[i]
    n_i   = fwd.metadata.state_dims[i]
    off_s = fwd.metadata.state_offsets[i]
    x0_i  = fwd.initial_state[off_s+1:off_s+n_i]
    dyn_i = fwd.dynamics.player_dynamics[i]

    return ForwardDiff.gradient(u_flat) do uv
        Ui = reshape(uv, m_i, N)
        Xi = _yp_rollout(dyn_i, x0_i, Ui, dt)
        J  = zero(eltype(uv))
        for k in 1:N
            J += evaluate_stage_cost(obj.stage_cost, Xi[:, k], Ui[:, k], nothing, k)
        end
        J += evaluate_terminal_cost(obj.terminal_cost, Xi[:, N+1], nothing)
        J * obj.scaling
    end
end

# ============================================================================
# Public solve overload
# ============================================================================

"""
    solve(game::ConvexGameProblem{T}, solver::YiPavel; kwargs...) -> GNEPSolution{T}

Compute a variational GNE via Algorithm 1 of Yi & Pavel (2017).

# Keyword Arguments
- `cost_fns::Vector{Function}` — Per-player objective closures
  `fᵢ(seg₁,…,segₙ) -> Real` where `segⱼ = vec(Uⱼ)` is player j's
  vectorised control trajectory.  Required when objectives couple across
  players.  Falls back to the LQ objectives in the `ConvexGameProblem` when
  not provided.

- `coupling_A::Matrix{T}` — Full coupling matrix `A ∈ ℝ^{m × Σdⱼ}`, where
  columns are ordered `[A₁ | A₂ | … | Aₙ]` and `dⱼ = mⱼ × N`.
- `coupling_b::Vector{T}` — Coupling RHS `b ∈ ℝ^m`.
- `coupling_leq::Bool = true` — `true` ⟹ coupling is `Ax ≤ b` (default,
  natural for capacity constraints); `false` ⟹ `Ax ≥ b` (paper form).

- `lb::Vector{Vector{T}}` — Per-player lower bounds on `vec(Uᵢ)`.
- `ub::Vector{Vector{T}}` — Per-player upper bounds on `vec(Uᵢ)`.

- `warmstart` — Previous `GNEPSolution` or `WarmstartData` (optional).
- `verbose::Bool = false` — Print per-iteration diagnostics.
"""
function solve(
    game        ::ConvexGameProblem{T},
    solver      ::YiPavel;
    warmstart   ::Union{Nothing, GNEPSolution{T}, WarmstartData{T}} = nothing,
    verbose     ::Bool                         = false,
    cost_fns    ::Union{Nothing, Vector{Function}}  = nothing,
    coupling_A  ::Union{Nothing, Matrix{T}}    = nothing,
    coupling_b  ::Union{Nothing, Vector{T}}    = nothing,
    coupling_leq::Bool                         = true,
    lb          ::Union{Nothing, Vector{Vector{T}}} = nothing,
    ub          ::Union{Nothing, Vector{Vector{T}}} = nothing,
    kwargs...
) where {T}
    fwd = game.forward_game
    np  = n_players(game)
    N   = n_steps(game)
    dt  = fwd.time_horizon.dt
    tf  = fwd.time_horizon.tf

    m_i = [control_dim(game, i) for i in 1:np]   # control dim per player
    n_i = [state_dim(game, i)   for i in 1:np]   # state dim per player
    d_i = m_i .* N                                # vectorised control dim

    # ── Primal initialisation ─────────────────────────────────────────────────
    x = if warmstart isa WarmstartData{T} && warmstart.trajectories !== nothing
        [vec(warmstart.trajectories[i].controls) for i in 1:np]
    elseif warmstart isa GNEPSolution{T}
        [vec(warmstart.trajectories[i].controls) for i in 1:np]
    else
        [zeros(T, d_i[i]) for i in 1:np]
    end

    # ── Coupling constraint (internally: Ax ≥ b) ─────────────────────────────
    A_blk, b_eff, m_con = _yp_parse_coupling(
        coupling_A, coupling_b, coupling_leq, np, d_i, T
    )

    # ── Box constraints ───────────────────────────────────────────────────────
    lb_i = lb !== nothing ? lb : [fill(T(-Inf), d_i[i]) for i in 1:np]
    ub_i = ub !== nothing ? ub : [fill(T( Inf), d_i[i]) for i in 1:np]

    # ── Dual variables ────────────────────────────────────────────────────────
    z = [zeros(T, m_con) for _ in 1:np]
    λ = [zeros(T, m_con) for _ in 1:np]

    τ_v = T(solver.τ)
    ν_v = T(solver.ν)
    σ_v = T(solver.σ)

    # ── Working arrays (pre-allocated) ────────────────────────────────────────
    x_new = [copy(xi) for xi in x]
    z_new = [copy(zi) for zi in z]
    λ_new = [copy(λi) for λi in λ]

    converged = false
    iters     = 0
    Δ_hist    = T[]
    t_start   = time()

    for iter in 1:solver.max_iter
        iters = iter

        # ── Primal update (simultaneous across all players) ───────────────────
        for i in 1:np
            gi = if cost_fns !== nothing
                _yp_gradient_coupled(cost_fns[i], x, i, d_i)
            else
                _yp_gradient_lq(fwd, i, x[i], N, dt)
            end
            # Projected gradient descent: x_new = P_{Ωᵢ}[x − τ(∇f − Aᵢᵀλᵢ)]
            x_new[i] = _yp_project_box(
                x[i] .- τ_v .* (gi .- A_blk[i]' * λ[i]),
                lb_i[i], ub_i[i]
            )
        end

        # ── Auxiliary update (Laplacian consensus on λ) ───────────────────────
        for i in 1:np
            z_new[i] = z[i] .+ ν_v .* _yp_laplacian_row(i, λ)
        end

        # ── Dual update ────────────────────────────────────────────────────────
        for i in 1:np
            Lλ_i        = _yp_laplacian_row(i, λ)
            Lz_new_i    = _yp_laplacian_row(i, z_new)
            Lz_i        = _yp_laplacian_row(i, z)
            consensus_z = 2 .* Lz_new_i .- Lz_i

            # Extrapolated local constraint residual (Aᵢ(2x_new − x) − b)
            local_term  = A_blk[i] * (2 .* x_new[i] .- x[i]) .- b_eff

            λ_new[i] = _yp_project_nonneg(
                λ[i] .- σ_v .* (local_term .+ consensus_z .+ Lλ_i)
            )
        end

        # ── Convergence check ─────────────────────────────────────────────────
        Δ = maximum(maximum(abs.(x_new[i] .- x[i])) for i in 1:np)
        push!(Δ_hist, Δ)

        verbose && @printf("[YiPavel %4d] Δx = %.3e  max_λ = %.3e\n",
                           iter, Δ, maximum(maximum(λi) for λi in λ))

        # Swap buffers
        x, x_new = x_new, x
        z, z_new = z_new, z
        λ, λ_new = λ_new, λ

        Δ < T(solver.tol) && (converged = true; break)
    end

    return _yp_assemble(
        game, fwd, x, λ, cost_fns,
        converged, iters, time() - t_start, Δ_hist,
        N, dt, tf, m_i, n_i, d_i
    )
end

# ============================================================================
# Private helpers
# ============================================================================

"""
    _yp_parse_coupling(A, b, leq, np, d_i, T)

Parse coupling keyword arguments into per-player block matrices (Ax ≥ b form).
Returns `(A_blk, b_eff, m_con)`.
"""
function _yp_parse_coupling(A, b, leq, np, d_i, ::Type{T}) where {T}
    if A !== nothing && b !== nothing
        m_con   = length(b)
        offsets = [0; cumsum(d_i)]
        A_blk   = [A[:, offsets[i]+1:offsets[i+1]] for i in 1:np]
        b_eff   = copy(b)
        if leq
            # Ax ≤ b  →  −Ax ≥ −b
            A_blk = [.-ai for ai in A_blk]
            b_eff  = .-b_eff
        end
        # Each player owns an equal share of the global RHS (Yi & Pavel §III.A):
        #   b = Σᵢ bᵢ  where  bᵢ = b/N  (even split)
        b_eff = b_eff ./ np
        return A_blk, b_eff, m_con
    else
        # No coupling — trivial (b = −∞ means constraint is always satisfied)
        m_con  = 1
        A_blk  = [zeros(T, 1, d_i[i]) for i in 1:np]
        b_eff  = fill(T(-Inf), 1)
        return A_blk, b_eff, m_con
    end
end

"""
    _yp_assemble(game, fwd, x, λ, cost_fns, ...)

Assemble a `GNEPSolution` from the converged primal/dual iterates.
"""
function _yp_assemble(
    game      ::ConvexGameProblem{T},
    fwd       ::GameProblem{T},
    x         ::Vector{Vector{T}},
    λ         ::Vector{Vector{T}},
    cost_fns  ::Union{Nothing, Vector{Function}},
    converged ::Bool,
    iters     ::Int,
    wall_time ::Float64,
    Δ_hist    ::Vector{T},
    N         ::Int,
    dt        ::T,
    tf        ::T,
    m_i       ::Vector{Int},
    n_i       ::Vector{Int},
    d_i       ::Vector{Int}
) where {T}
    np    = n_players(game)
    times = collect(range(T(0), tf, length=N+1))

    trajectories = map(1:np) do i
        Ui    = reshape(x[i], m_i[i], N)
        dyn_i = fwd.dynamics.player_dynamics[i]
        off_s = fwd.metadata.state_offsets[i]
        x0_i  = fwd.initial_state[off_s+1:off_s+n_i[i]]
        Xi    = _yp_rollout(dyn_i, x0_i, Ui, dt)

        Ji = if cost_fns !== nothing
            Float64(cost_fns[i]([x[j] for j in 1:np]...))
        else
            obj = get_objective(fwd, i)
            J   = zero(T)
            for k in 1:N
                J += evaluate_stage_cost(obj.stage_cost, Xi[:, k], Ui[:, k], nothing, k)
            end
            J += evaluate_terminal_cost(obj.terminal_cost, Xi[:, N+1], nothing)
            J * obj.scaling
        end

        Trajectory{T}(i, Xi, Ui, times, T(Ji))
    end

    strategy = OpenLoopStrategy(
        [reshape(x[i], m_i[i], N) for i in 1:np],
        m_i,
        times
    )

    info = Dict{Symbol, Any}(
        :λ_final => λ,
        :Δ_hist  => Δ_hist,
    )

    return GNEPSolution(
        fwd, trajectories;
        strategy         = strategy,
        equilibrium_type = :GeneralizedNash,
        converged        = converged,
        iterations       = iters,
        solve_time       = wall_time,
        solver_info      = info
    )
end
